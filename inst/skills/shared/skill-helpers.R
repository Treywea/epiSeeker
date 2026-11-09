## Helpers shared by the scripts of the epiSeeker agent skills (inst/skills/<skill>/scripts/).
## They hold no analysis: every number in the outputs comes from an epiSeeker function. What
## they do is the bookkeeping every skill must do the same way -- read arguments, make the
## fixed output directories, record each check (PASS / NOTE / WARN / STOP), capture warnings
## instead of letting them scroll by, call interpret() last, and write the provenance.
##
## A script sources this file with
##   source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
## so it works from the installed package, also when the skill folder was copied elsewhere.

suppressPackageStartupMessages(library(epiSeeker))

.skill <- new.env(parent = emptyenv())

## `--key value` and `--flag` arguments; a harness can pass them with
## options(epiSeeker.skill.args = c(...)) instead of the command line.
skill_args <- function(defaults = list()) {
    a <- getOption("epiSeeker.skill.args", commandArgs(trailingOnly = TRUE))
    out <- defaults
    i <- 1
    while (i <= length(a)) {
        k <- a[i]
        if (!startsWith(k, "--")) stop("unexpected argument '", k, "'; use --key value")
        k <- gsub("-", "_", sub("^--", "", k))
        if (i < length(a) && !startsWith(a[i + 1], "--")) {
            out[[k]] <- a[i + 1]; i <- i + 2
        } else {
            out[[k]] <- TRUE; i <- i + 1
        }
    }
    out
}

## comma-separated argument -> character vector; NULL stays NULL
skill_split <- function(x) if (is.null(x) || isTRUE(x)) NULL else trimws(strsplit(x, ",", fixed = TRUE)[[1]])

## Start a run: fixed output directories, empty check table, the arguments.
skill_start <- function(skill, args, out) {
    if (is.null(out)) stop("give --out <directory>")
    dirs <- file.path(out, c("figures", "tables", "report", "provenance"))
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    .skill$skill <- skill
    .skill$out <- normalizePath(out)
    .skill$args <- args
    .skill$t0 <- Sys.time()
    .skill$checks <- data.frame(step = character(0), check = character(0), status = character(0),
                                detail = character(0), stringsAsFactors = FALSE)
    .skill$warnings <- character(0)
    .skill$inputs <- character(0)
    .skill$notes <- list()
    .skill$interpret <- "not run"
    invisible(.skill$out)
}

skill_path <- function(...) file.path(.skill$out, ...)

## Record one check. STOP writes the provenance and stops: the analysis must not go on with a
## silent error in its input. The condition has class "skill_stop" so a harness can tell a
## check that fired from a crash.
skill_check <- function(step, check, status = c("PASS", "NOTE", "WARN", "STOP"), detail = "") {
    status <- match.arg(status)
    .skill$checks[nrow(.skill$checks) + 1, ] <- list(step, check, status, detail)
    message(sprintf("[%s] %s: %s%s", status, step, check, if (nzchar(detail)) paste0(" -- ", detail) else ""))
    if (status == "STOP") {
        skill_report(paste(.skill$skill, "report (stopped)"),
                     list(Stopped = paste0("A check stopped the run at ", step, ": ", check, " -- ", detail,
                                           ". Nothing after this point was computed.")))
        skill_finish(stopped = TRUE)
        stop(structure(class = c("skill_stop", "error", "condition"),
                       list(message = paste0("STOP at ", step, ": ", check, " -- ", detail), call = NULL)))
    }
    invisible(status)
}

## Run `expr`, keep its warnings in the record (and as a WARN check) instead of losing them.
skill_capture <- function(expr, step) {
    withCallingHandlers(expr, warning = function(w) {
        msg <- conditionMessage(w)
        .skill$warnings <- c(.skill$warnings, paste0(step, ": ", msg))
        skill_check(step, "warning raised", "WARN", substr(gsub("\n", " ", msg), 1, 300))
        invokeRestart("muffleWarning")
    })
}

## Record an input file with its md5.
skill_input <- function(path, role) {
    for (p in path) {
        md5 <- if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
        .skill$inputs <- c(.skill$inputs, paste(role, p, md5, sep = "\t"))
    }
    invisible(path)
}

skill_note <- function(key, value) { .skill$notes[[key]] <- value; invisible(value) }

## TxDb / OrgDb / BSgenome objects from package names, so that a run names its annotation.
skill_load <- function(pkg, what = pkg) {
    if (is.null(pkg)) return(NULL)
    if (!requireNamespace(pkg, quietly = TRUE)) stop("package '", pkg, "' is not installed")
    getExportedValue(pkg, what)
}

## ---- input checks shared by the skills ----------------------------------------------------

## The genome the user declares must be the genome of the TxDb. Nothing else catches a peak
## file of one build annotated with a TxDb of another when the peaks record no genome.
skill_check_genome <- function(declared, TxDb, step = "input") {
    if (is.null(declared) || isTRUE(declared)) {
        skill_check(step, "genome build declared", "STOP",
                    "give --genome (e.g. hg38, mm10): the build of the input coordinates, as stated by whoever produced them")
    }
    gt <- tryCatch(unique(stats::na.omit(GenomeInfoDb::genome(TxDb))), error = function(e) character(0))
    if (!length(gt)) {
        skill_check(step, "genome of TxDb", "WARN", "TxDb records no genome; cannot compare with --genome")
    } else if (!declared %in% gt) {
        skill_check(step, "declared genome matches TxDb", "STOP",
                    paste0("--genome ", declared, " but TxDb is ", paste(gt, collapse = ",")))
    } else {
        skill_check(step, "declared genome matches TxDb", "PASS", paste0(declared))
    }
}

## The genome recorded on the regions (often none) must not contradict the declared one.
skill_check_recorded_genome <- function(gr, declared, what = "regions", step = "input") {
    gp <- unique(stats::na.omit(GenomeInfoDb::genome(gr)))
    if (!length(gp)) {
        skill_check(step, paste("genome recorded on", what), "NOTE",
                    "none recorded (usual for bed files); relying on --genome and the out-of-bound check")
    } else if (!all(gp %in% declared)) {
        skill_check(step, paste("genome recorded on", what), "STOP",
                    paste0(paste(gp, collapse = ","), " recorded, --genome ", declared))
    } else {
        skill_check(step, paste("genome recorded on", what), "PASS", paste(gp, collapse = ","))
    }
}

## Sequence names: a different style ('1' vs 'chr1') drops every region silently in some tools
## and all of them here. Fixed only on request, and the fix is recorded.
skill_check_seqlevels <- function(gr, TxDb, fix = FALSE, what = "regions", step = "input") {
    st <- tryCatch(GenomeInfoDb::seqlevelsStyle(gr), error = function(e) NA_character_)
    tt <- tryCatch(GenomeInfoDb::seqlevelsStyle(TxDb), error = function(e) NA_character_)
    if (!any(is.na(st)) && !any(is.na(tt)) && !length(intersect(st, tt))) {
        if (!isTRUE(fix)) {
            skill_check(step, paste("seqlevels style of", what), "STOP",
                        paste0(paste(st, collapse = "/"), " vs TxDb ", paste(tt, collapse = "/"),
                               "; re-run with --fix_seqlevels to convert (recorded) or fix the input"))
        }
        GenomeInfoDb::seqlevelsStyle(gr) <- tt[1]
        skill_check(step, paste("seqlevels style of", what), "NOTE",
                    paste0("converted ", paste(st, collapse = "/"), " -> ", tt[1], " (--fix_seqlevels)"))
    } else {
        skill_check(step, paste("seqlevels style of", what), "PASS", paste(st, collapse = "/"))
    }
    used <- unique(as.character(GenomicRanges::seqnames(gr)))
    miss <- setdiff(used, GenomeInfoDb::seqlevels(TxDb))
    if (length(miss)) {
        n <- sum(as.character(GenomicRanges::seqnames(gr)) %in% miss)
        skill_check(step, paste("sequences of", what, "absent from TxDb"), "NOTE",
                    paste0(n, " of ", length(gr), " on ", paste(utils::head(miss, 5), collapse = ", "),
                           if (length(miss) > 5) ", ..." else "", " (annotateSeq() reports them as dropped)"))
    }
    gr
}

## What annotateSeq() found (annoChecks): a suspect build stops the run unless the user
## accepts it explicitly; dropped regions are recorded.
skill_check_anno <- function(anno, what, allow_suspect = FALSE, step = "correctness") {
    ac <- annoChecks(anno)
    if (is.null(ac)) {
        skill_check(step, paste("annoChecks of", what), "WARN", "no checks recorded (object made by an old version)")
        return(invisible(NULL))
    }
    b <- ac$build
    if (isTRUE(b$suspect)) {
        skill_check(step, paste("genome build of", what), if (isTRUE(allow_suspect)) "WARN" else "STOP",
                    paste(c(b$reasons, if (!isTRUE(allow_suspect)) "fix the build (liftOver or the right TxDb); --allow_suspect only if you know why"), collapse = "; "))
    } else if (!isTRUE(b$reliable)) {
        skill_check(step, paste("genome build of", what), "NOTE",
                    paste0("no peak beyond chromosome ends, but only ", b$n_peaks,
                           " regions: the build check is unreliable below 200"))
    } else {
        skill_check(step, paste("genome build of", what), "PASS",
                    paste0("0 of ", b$n_peaks, " regions beyond the chromosome ends of TxDb"))
    }
    d <- ac$dropped
    if (!is.null(d) && d$n > 0) {
        bs <- d$by_seqlevel
        skill_check(step, paste("regions dropped from", what), if (d$n / d$n_total > 0.05) "WARN" else "NOTE",
                    paste0(d$n, " of ", d$n_total, ": ",
                           paste0(utils::head(bs$seqlevel, 3), " (", utils::head(bs$n, 3), ", ",
                                  utils::head(bs$reason, 3), ")", collapse = "; ")))
    } else {
        skill_check(step, paste("regions dropped from", what), "PASS", "none")
    }
    invisible(ac)
}

## ---- outputs ---------------------------------------------------------------------------

skill_table <- function(df, name) {
    f <- skill_path("tables", name)
    utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(f)
}

skill_figure <- function(p, name, width = 7, height = 5) {
    f <- skill_path("figures", name)
    grDevices::pdf(f, width = width, height = height)
    on.exit(grDevices::dev.off())
    print(p)
    invisible(f)
}

## interpret() is the last step. Without a model (--model or EPISEEKER_MODEL) it is skipped and
## the report says so; it never falls back to a default provider.
skill_interpret <- function(x, model, ...) {
    if (is.null(model) || isTRUE(model) || !nzchar(model)) model <- Sys.getenv("EPISEEKER_MODEL")
    if (!nzchar(model)) {
        .skill$interpret <- "skipped: no model (give --model provider:model or set EPISEEKER_MODEL)"
        skill_check("interpret", "interpret() run", "NOTE", .skill$interpret)
        return(invisible(NULL))
    }
    res <- skill_capture(clusterProfiler::interpret(x, model = model, ...), "interpret")
    writeLines(utils::capture.output(print(res)), skill_path("report", "interpretation.txt"))
    writeLines(attr(res, "evidence")$text, skill_path("report", "evidence_digest.txt"))
    prov <- clusterProfiler::provenance(res)
    writeLines(utils::capture.output(print(prov)), skill_path("provenance", "interpret_provenance.txt"))
    saveRDS(res, skill_path("provenance", "interpretation.rds"))
    ok <- isTRUE(prov$ok)
    .skill$interpret <- paste0(if (ok) "ok" else "INCOMPLETE", " (model ", model, ")")
    skill_check("interpret", "model returned a complete answer", if (ok) "PASS" else "WARN",
                if (ok) "" else "fields are NA; see provenance/interpret_provenance.txt")
    invisible(res)
}

## Phrases the skills do not allow in an interpretation, flagged for the reader (the text of
## the model is never edited). `patterns` is a named character vector of regular expressions.
skill_flag_phrases <- function(res, patterns, step = "interpret") {
    if (is.null(res)) return(invisible(character(0)))
    txt <- unlist(res, use.names = FALSE)
    txt <- txt[!is.na(txt)]
    hits <- character(0)
    for (nm in names(patterns)) {
        m <- grep(patterns[[nm]], txt, ignore.case = TRUE, value = TRUE, perl = TRUE)
        if (length(m)) hits <- c(hits, paste0(nm, ": \"", substr(m[1], 1, 160), "\""))
    }
    if (length(hits)) {
        skill_check(step, "forbidden statements in the interpretation", "WARN", paste(hits, collapse = " | "))
    } else {
        skill_check(step, "forbidden statements in the interpretation", "PASS", "none found")
    }
    invisible(hits)
}

## Report: the sections given (markdown), then the check table, the warnings and the
## interpret() status. Every number in `sections` must come from a package function.
skill_report <- function(title, sections) {
    ck <- .skill$checks
    lines <- c(paste0("# ", title), "",
               paste0("Skill `", .skill$skill, "`, epiSeeker ", as.character(utils::packageVersion("epiSeeker")),
                      ", ", format(.skill$t0, "%Y-%m-%d %H:%M:%S")), "",
               unlist(lapply(names(sections), function(s) c(paste0("## ", s), "", sections[[s]], ""))),
               "## Checks", "", "| step | check | status | detail |", "|---|---|---|---|",
               paste0("| ", ck$step, " | ", ck$check, " | ", ck$status, " | ", gsub("\\|", "/", ck$detail), " |"), "",
               "## interpret()", "", .skill$interpret, "")
    writeLines(lines, skill_path("report", "report.md"))
}

## Provenance: what ran, on what, with which package versions and annotation, what the checks
## found and which files were written (with md5). Written also when a check stops the run.
skill_finish <- function(stopped = FALSE) {
    pkgs <- c("epiSeeker", "clusterProfiler", "aisdk", "GenomicRanges", "GenomeInfoDb")
    ver <- vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), "")
    desc <- utils::packageDescription("epiSeeker")
    sha <- desc$RemoteSha
    if (is.null(sha)) sha <- desc$GithubSHA1
    utils::write.table(.skill$checks, skill_path("provenance", "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    outs <- list.files(.skill$out, recursive = TRUE)
    outs <- outs[!grepl("^provenance/(run|checks)\\.", outs)]
    md5 <- unname(tools::md5sum(file.path(.skill$out, outs)))
    a <- .skill$args
    lines <- c(
        paste0("skill\t", .skill$skill),
        paste0("status\t", if (stopped) "STOPPED by a check" else "completed"),
        paste0("started\t", format(.skill$t0, "%Y-%m-%d %H:%M:%S %Z")),
        paste0("seconds\t", round(as.numeric(difftime(Sys.time(), .skill$t0, units = "secs")), 1)),
        paste0("R\t", R.version.string),
        paste0("version_", pkgs, "\t", ver),
        paste0("epiSeeker_commit\t", if (is.null(sha)) "not recorded by the installation" else sha),
        paste0("arg_", names(a), "\t", vapply(a, function(v) paste(format(v), collapse = ","), "")),
        paste0("input\t", .skill$inputs),
        paste0("note_", names(.skill$notes), "\t", vapply(.skill$notes, function(v) paste(format(v), collapse = "; "), "")),
        paste0("warning\t", gsub("\n", " ", .skill$warnings)),
        paste0("interpret\t", .skill$interpret),
        paste0("output\t", outs, "\t", md5)
    )
    writeLines(lines, skill_path("provenance", "run.txt"))
    writeLines(utils::capture.output(utils::sessionInfo()), skill_path("provenance", "sessionInfo.txt"))
    invisible(lines)
}
