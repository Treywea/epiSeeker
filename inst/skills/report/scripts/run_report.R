## report: collect the outputs and provenance of epiSeeker skill runs into one report -> input
## checks (each folder is a skill run; its outputs are unchanged since the run) -> run table,
## every STOP / WARN row, interpret() status, cross-run consistency -> provenance. No model is
## called and no number is computed here: everything comes from the runs. See ../SKILL.md.
##
## Rscript run_report.R --runs run1,run2,... | --runs_dir <folder of runs> --out <dir> [--title "..."]
## Rscript run_report.R --example --out <dir>   (runs two small skills and one that stops, then collects them)

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

a <- skill_args(list(title = "epiSeeker analysis report"))
if (isTRUE(a$example)) {
    ## three runs without a model: base-methylation and base-singlemolecule examples, and
    ## base-methylation on the bundled table with levels in percent, which a check stops. The
    ## child R sessions get this session's library paths.
    sk <- system.file("skills", package = "epiSeeker")
    rdir <- file.path(a$out, "runs")
    dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
    pct <- utils::read.delim(system.file("extdata", "demo_bisseq.txt", package = "epiSeeker"))
    pct$Methylation <- pct$Methylation * 100
    utils::write.table(pct, file.path(rdir, "demo_bisseq_percent.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    jobs <- list(
        "base-methylation" = c(file.path(sk, "base-methylation/scripts/run_base_methylation.R"), "--example"),
        "base-singlemolecule" = c(file.path(sk, "base-singlemolecule/scripts/run_base_singlemolecule.R"), "--example"),
        "base-methylation_percent" = c(file.path(sk, "base-methylation/scripts/run_base_methylation.R"),
                                       "--input", file.path(rdir, "demo_bisseq_percent.txt"), "--format", "table",
                                       "--genome", "hg38", "--names", "percent_levels"))
    for (j in names(jobs)) {
        system2(file.path(R.home("bin"), "Rscript"), c(jobs[[j]], "--out", file.path(rdir, j)),
                stdout = file.path(rdir, paste0(j, ".log")), stderr = file.path(rdir, paste0(j, ".log")),
                env = c("EPISEEKER_MODEL=", paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))))
    }
    a$runs <- paste(file.path(rdir, names(jobs)), collapse = ",")
    a$title <- "epiSeeker report of the example runs"
}
skill_start("report", a, a$out)

## ---- 1. input checks -------------------------------------------------------------------
runs <- skill_split(a$runs)
if (!is.null(a$runs_dir)) runs <- c(runs, list.dirs(a$runs_dir, recursive = FALSE))
if (!length(runs)) skill_check("input", "runs given", "STOP", "give --runs dir1,dir2 or --runs_dir <folder>")
runs <- normalizePath(runs, mustWork = FALSE)
names(runs) <- make.unique(basename(runs))
for (r in names(runs)) {
    ok <- file.exists(file.path(runs[[r]], "provenance", c("run.txt", "checks.tsv")))
    skill_check("input", paste("is a skill run:", r), if (all(ok)) "PASS" else "STOP",
                if (all(ok)) runs[[r]] else paste(runs[[r]], "has no provenance/run.txt or checks.tsv"))
}
## run.txt: key <tab> value [<tab> md5]
read_run <- function(d) {
    f <- file.path(d, "provenance", "run.txt")
    skill_input(f, "run")
    p <- strsplit(readLines(f, warn = FALSE), "\t", fixed = TRUE)
    key <- vapply(p, `[`, "", 1); val <- vapply(p, function(x) if (length(x) >= 2) x[2] else "", "")
    md5 <- vapply(p, function(x) if (length(x) >= 3) x[3] else NA_character_, "")
    one <- function(k) { v <- val[key == k]; if (length(v)) v[1] else NA_character_ }
    out <- key == "output"
    now <- unname(tools::md5sum(file.path(d, val[out])))
    list(skill = one("skill"), status = one("status"), seconds = as.numeric(one("seconds")),
         version = one("version_epiSeeker"), commit = one("epiSeeker_commit"), genome = one("arg_genome"),
         txdb = one("arg_txdb"), interpret = one("interpret"), outputs = val[out],
         changed = val[out][is.na(now) | now != md5[out]])
}
info <- lapply(runs, read_run)
ck <- do.call(rbind, lapply(names(runs), function(r) {
    x <- utils::read.delim(file.path(runs[[r]], "provenance", "checks.tsv"), quote = "", colClasses = "character")
    if (nrow(x)) cbind(run = r, x) else NULL
}))
for (r in names(runs)) {
    ch <- info[[r]]$changed
    skill_check("input", paste("outputs unchanged since the run:", r), if (length(ch)) "WARN" else "PASS",
                if (length(ch)) paste0(length(ch), " of ", length(info[[r]]$outputs), " files changed or missing: ",
                                       paste(utils::head(ch, 3), collapse = ", ")) else
                    paste(length(info[[r]]$outputs), "files, md5 as recorded"))
}

## ---- 2. cross-run checks and tables --------------------------------------------------------
v <- unique(stats::na.omit(vapply(info, `[[`, "", "version")))
skill_check("consistency", "runs used one epiSeeker version", if (length(v) > 1) "WARN" else "PASS", paste(v, collapse = ", "))
g <- vapply(info, `[[`, "", "genome")
gu <- unique(stats::na.omit(g))
skill_check("consistency", "runs declare one genome", if (length(gu) > 1) "WARN" else "PASS",
            paste0(paste(names(g), ifelse(is.na(g), "none", g), sep = ": ", collapse = "; ")))
tx <- unique(stats::na.omit(vapply(info, `[[`, "", "txdb")))
if (length(tx) > 1) skill_check("consistency", "runs used one TxDb", "NOTE", paste(tx, collapse = ", "))
st <- vapply(info, `[[`, "", "status")
stopped <- names(st)[grepl("^STOPPED", st)]
skill_check("consistency", "runs stopped by a check", if (length(stopped)) "WARN" else "PASS",
            if (length(stopped)) paste0(paste(stopped, collapse = ", "), ": their analyses stopped; no result of theirs is in this report")
            else paste(length(st), "runs completed"))
it <- vapply(info, `[[`, "", "interpret")
inc <- names(it)[grepl("^INCOMPLETE", it)]
skill_check("consistency", "interpret() answers complete", if (length(inc)) "WARN" else "PASS",
            if (length(inc)) paste(inc, collapse = ", ") else "no incomplete answer")
fb <- if (is.null(ck)) NULL else ck[ck$check == "forbidden statements in the interpretation" & ck$status == "WARN", , drop = FALSE]
skill_check("consistency", "interpretations without flagged statements", if (!is.null(fb) && nrow(fb)) "WARN" else "PASS",
            if (!is.null(fb) && nrow(fb)) paste(fb$run, collapse = ", ") else "none flagged")

cnt <- function(r, s) if (is.null(ck)) 0L else sum(ck$run == r & ck$status == s)
rt <- data.frame(run = names(runs), skill = vapply(info, `[[`, "", "skill"), status = unname(st),
                 PASS = vapply(names(runs), cnt, 0L, s = "PASS"), NOTE = vapply(names(runs), cnt, 0L, s = "NOTE"),
                 WARN = vapply(names(runs), cnt, 0L, s = "WARN"), STOP = vapply(names(runs), cnt, 0L, s = "STOP"),
                 interpret = sub("^skipped:.*$", "skipped (no model)", sub(" \\(model.*$", "", unname(it))), epiSeeker = vapply(info, `[[`, "", "version"),
                 genome = unname(g), outputs_changed = lengths(lapply(info, `[[`, "changed")),
                 seconds = vapply(info, `[[`, 0, "seconds"), row.names = NULL)
skill_table(rt, "runs.tsv")
if (!is.null(ck)) {
    skill_table(ck, "checks_all.tsv")
    att <- ck[ck$status %in% c("STOP", "WARN"), , drop = FALSE]
    skill_table(att, "attention.tsv")
    long <- as.data.frame(table(run = ck$run, status = factor(ck$status, c("PASS", "NOTE", "WARN", "STOP"))))
    skill_figure(ggplot2::ggplot(long, ggplot2::aes(x = .data$run, y = .data$Freq, fill = .data$status)) +
                     ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_minimal() +
                     ggplot2::scale_fill_manual(values = c(PASS = "grey70", NOTE = "#8fb3d9", WARN = "#e0a030", STOP = "#c0392b")) +
                     ggplot2::labs(x = NULL, y = "checks", fill = NULL), "check_status.pdf", width = 7, height = 1.5 + 0.35 * length(runs))
}

## ---- 3. report: the runs' own words, with the checks first ---------------------------------
excerpt <- function(r) {
    f <- file.path(runs[[r]], "report", "report.md")
    if (!file.exists(f)) return("(no report.md)")
    l <- readLines(f, warn = FALSE)
    cut <- which(l == "## Checks")
    if (length(cut)) l <- l[seq_len(cut[1] - 1)]
    l <- l[-1]                                      # its title, replaced by the run heading
    sub("^(#+) ", "\\1## ", l)                       # demote headings under the run heading
}
md <- function(x) gsub("\\|", "/", x)
att_lines <- if (is.null(ck) || !nrow(att)) "none" else
    c("| run | step | check | status | detail |", "|---|---|---|---|---|",
      paste0("| ", att$run, " | ", att$step, " | ", md(att$check), " | ", att$status, " | ", md(substr(att$detail, 1, 300)), " |"))
sec <- c(list(
    "Runs" = c("| run | skill | status | PASS | NOTE | WARN | STOP | interpret() | genome |", "|---|---|---|---|---|---|---|---|---|",
               paste0("| ", rt$run, " | ", rt$skill, " | ", rt$status, " | ", rt$PASS, " | ", rt$NOTE, " | ", rt$WARN, " | ",
                      rt$STOP, " | ", rt$interpret, " | ", rt$genome, " |"),
               "", paste0("Folders: ", paste0("`", runs, "`", collapse = ", "))),
    "Needs attention (every STOP and WARN of the runs)" = att_lines),
    stats::setNames(lapply(names(runs), function(r) c(
        paste0("Folder `", runs[[r]], "`; status ", st[[r]], "; interpret(): ", it[[r]], "."),
        if (file.exists(file.path(runs[[r]], "report", "interpretation.txt")))
            paste0("Model answer: `", file.path(runs[[r]], "report", "interpretation.txt"), "` (with its digest in evidence_digest.txt)"),
        "", excerpt(r))), paste0("Run ", names(runs), " (", vapply(info, `[[`, "", "skill"), ")")))

## ---- 4. interpret(): not called here -------------------------------------------------------
.skill$interpret <- "not called: the report collects the interpret() results of the runs"
skill_check("interpret", "interpret() not called by the report", "NOTE", "each run's answer and its provenance stay in that run's folder")
skill_report(a$title, sec)
saveRDS(list(runs = rt, checks = ck), skill_path("provenance", "objects.rds"))
skill_finish()
message("report done: ", skill_path("report", "report.md"))
