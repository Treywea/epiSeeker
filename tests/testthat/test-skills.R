library(epiSeeker)

context("agent skills in inst/skills: format, structure and the functions they name")

skills_dir <- system.file("skills", package = "epiSeeker")
skill_names <- setdiff(list.dirs(skills_dir, full.names = FALSE, recursive = FALSE), "shared")

frontmatter <- function(path) {
    l <- readLines(path, warn = FALSE)
    end <- which(l == "---")
    if (length(end) < 2 || end[1] != 1) return(NULL)
    fm <- l[2:(end[2] - 1)]
    kv <- regmatches(fm, regexec("^([a-z_-]+):\\s*(.*)$", fm))
    stats::setNames(lapply(kv, `[`, 3), vapply(kv, `[`, "", 2))
}

test_that("the nine skills are installed with SKILL.md, a script and a worked example", {
    expect_true(all(c("range-tf", "range-histone", "range-atac", "range-singlecell", "range-spatial",
                      "base-methylation", "base-singlemolecule", "interaction-loops", "report") %in% skill_names))
    expect_true(all(file.exists(file.path(skills_dir, "shared", c("skill-helpers.R", "range-steps.R", "group-steps.R")))))
    ## the ENCODE region sets of the examples, with their source and licence
    enc <- system.file("extdata", "encode", package = "epiSeeker")
    expect_true(all(file.exists(file.path(enc, c("SOURCE.txt", "LNCaP_R1881_DNase_ENCFF570MWX_chr1-2.bed.gz",
                                                 "prostate_H3K27me3_ENCFF198QRW_chr1-2.bed.gz")))))
    for (s in skill_names) {
        expect_true(file.exists(file.path(skills_dir, s, "SKILL.md")), info = s)
        expect_gte(length(list.files(file.path(skills_dir, s, "scripts"), pattern = "\\.R$")), 1)
        expect_true(file.exists(file.path(skills_dir, s, "references", "worked-example.md")), info = s)
    }
})

test_that("SKILL.md front matter follows the Agent Skills format", {
    for (s in skill_names) {
        fm <- frontmatter(file.path(skills_dir, s, "SKILL.md"))
        expect_false(is.null(fm), info = s)
        expect_identical(fm$name, s)
        expect_match(fm$name, "^[a-z0-9]+(-[a-z0-9]+)*$")
        expect_lte(nchar(fm$name), 64)
        expect_gt(nchar(fm$description), 50)
        expect_lte(nchar(fm$description), 1024)
    }
})

test_that("every SKILL.md has the unified structure, in order", {
    heads <- c("## 1. Input checks", "## 2. Analysis", "## 3. Correctness checks", "## 4. interpret() last",
               "## 5. Provenance", "## Do not say", "## Run it")
    for (s in skill_names) {
        l <- readLines(file.path(skills_dir, s, "SKILL.md"), warn = FALSE)
        pos <- vapply(heads, function(h) { i <- which(startsWith(l, h)); if (length(i)) i[1] else NA_integer_ }, 0L)
        expect_false(anyNA(pos), info = paste(s, paste(heads[is.na(pos)], collapse = ", ")))
        expect_false(is.unsorted(pos), info = s)
    }
})

## the functions a skill names must exist, so that a renamed or removed function breaks a test
## rather than an agent's run
known <- function(fn, pkg = NULL) {
    if (!is.null(pkg)) {
        if (!requireNamespace(pkg, quietly = TRUE)) return(TRUE)   # cannot check, not a failure
        return(fn %in% getNamespaceExports(pkg) || (pkg == "epiSeeker" && exists(fn, envir = asNamespace(pkg))))
    }
    fn %in% getNamespaceExports("epiSeeker") ||
        any(vapply(c("base", "stats", "utils", "methods", "grDevices", "tools"), function(p)
            exists(fn, envir = asNamespace(p), inherits = FALSE), TRUE))
}

## terminal tokens of an R file, in reading order
tokens <- function(f) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    pd <- pd[pd$terminal, , drop = FALSE]
    pd[order(pd$line1, pd$col1), c("token", "text")]
}
## names given a function with `<-` in a file
assigned <- function(tk) tk$text[tk$token == "SYMBOL" & c(tk$token[-1], "") == "LEFT_ASSIGN" &
                                 c(tk$token[-(1:2)], "", "") == "FUNCTION"]

test_that("functions called by the skill scripts exist", {
    helpers <- list.files(file.path(skills_dir, "shared"), pattern = "\\.R$", full.names = TRUE)
    own <- unlist(lapply(helpers, function(f) assigned(tokens(f))))
    expect_true(all(c("skill_check", "skill_interpret", "skill_finish", "range_inputs", "group_inputs") %in% own))
    for (f in c(helpers, list.files(file.path(skills_dir, skill_names, "scripts"), full.names = TRUE))) {
        tk <- tokens(f)
        defined <- c(own, assigned(tk))
        for (i in which(tk$token == "SYMBOL_FUNCTION_CALL")) {
            fn <- tk$text[i]
            ## pkg::fn is SYMBOL_PACKAGE, NS_GET, SYMBOL_FUNCTION_CALL
            pkg <- if (i > 2 && tk$token[i - 1] %in% c("NS_GET", "NS_GET_INT")) tk$text[i - 2] else NULL
            if (is.null(pkg) && fn %in% defined) next
            expect_true(known(fn, pkg), info = paste0(basename(f), ": ", if (!is.null(pkg)) paste0(pkg, "::"), fn))
        }
    }
})

test_that("functions named in SKILL.md exist", {
    allow <- c("metadata", "data")   # S4Vectors::metadata(), utils::data() as written in prose
    for (s in skill_names) {
        txt <- paste(readLines(file.path(skills_dir, s, "SKILL.md"), warn = FALSE), collapse = "\n")
        code <- unlist(regmatches(txt, gregexpr("`[^`\n]+`", txt)))
        m <- regmatches(code, regexpr("([A-Za-z][A-Za-z0-9.]*::)?[A-Za-z.][A-Za-z0-9._]*(?=\\()", code, perl = TRUE))
        for (x in unique(m)) {
            parts <- strsplit(x, "::", fixed = TRUE)[[1]]
            ok <- if (length(parts) == 2) known(parts[2], parts[1]) else x %in% allow || known(x)
            expect_true(ok, info = paste0(s, "/SKILL.md: ", x))
        }
    }
})

test_that("base-singlemolecule adjusts co-occurrence with modCooccurrence(adjust = \"read\"), not its own model", {
    f <- file.path(skills_dir, "base-singlemolecule", "scripts", "run_base_singlemolecule.R")
    tk <- tokens(f)
    calls <- tk$text[tk$token == "SYMBOL_FUNCTION_CALL"]
    expect_false(any(calls %in% c("glm", "mantelhaen.test")))
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_match(txt, "modCooccurrence\\([^)]*adjust = \"read\"[^)]*covariate = cv")
    ## the table keeps its columns
    for (col in c("odds_ratio_unadjusted", "p_unadjusted", "odds_ratio_adjusted", "p_adjusted"))
        expect_match(txt, col, fixed = TRUE)
    ## and the skill still asks for distance-matched control pairs
    md <- paste(readLines(file.path(skills_dir, "base-singlemolecule", "SKILL.md"), warn = FALSE), collapse = " ")
    expect_match(md, "distance-matched control pairs", fixed = TRUE)
})

test_that("range-tf, range-histone and range-atac run seqAssayQC() as an input check", {
    for (s in c("range-tf", "range-histone", "range-atac")) {
        f <- list.files(file.path(skills_dir, s, "scripts"), pattern = "\\.R$", full.names = TRUE)
        calls <- tokens(f)$text[tokens(f)$token == "SYMBOL_FUNCTION_CALL"]
        expect_true(all(c("range_assay_qc", "range_keep_assay_qc", "range_assay_section") %in% calls), info = s)
        md <- paste(readLines(file.path(skills_dir, s, "SKILL.md"), warn = FALSE), collapse = " ")
        expect_match(md, "seqAssayQC(", fixed = TRUE)
        expect_match(md, "The check is partial", fixed = TRUE)
        expect_match(md, "59%", fixed = TRUE)
    }
})

test_that("range_assay_qc() reports a flag as WARN and no flag as a NOTE that does not verify the assay", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    env <- new.env()
    sys.source(file.path(skills_dir, "shared", "skill-helpers.R"), envir = env)
    sys.source(file.path(skills_dir, "shared", "range-steps.R"), envir = env)
    out <- tempfile("assayqc")
    on.exit(unlink(out, recursive = TRUE))
    suppressMessages(env$skill_start("test", list(), out))
    ## 600 domains of 20 kb declared as TF: flagged (median width); the ENCODE H3K27me3 example
    ## declared as TF: a broad mark the check cannot tell from a TF, so no flag
    dom <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1e6, by = 2e5, length.out = 600), width = 20000))
    k27 <- readPeakFile(system.file("extdata", "encode", "prostate_H3K27me3_ENCFF198QRW_chr1-2.bed.gz", package = "epiSeeker"))
    qc <- suppressMessages(env$range_assay_qc(list(domains = dom, k27 = k27), txdb, "TF"))
    expect_false(qc$domains$consistent)
    expect_true(qc$k27$consistent)
    ck <- env$.skill$checks
    expect_equal(ck$status[ck$check == "assay label of domains (seqAssayQC)"], "WARN")
    expect_match(ck$detail[ck$check == "assay label of domains (seqAssayQC)"], "does not look like TF: median width")
    expect_equal(ck$status[ck$check == "assay label of k27 (seqAssayQC)"], "NOTE")
    expect_match(ck$detail[ck$check == "assay label of k27 (seqAssayQC)"], "no flag does not verify the assay", fixed = TRUE)
    expect_false(any(ck$status == "PASS" & grepl("seqAssayQC", ck$check)))
    tab <- utils::read.delim(file.path(out, "tables", "assay_qc.tsv"))
    expect_equal(tab$flagged, c(TRUE, FALSE))
    ## kept with the annotation interpret() reads; no --assay: a NOTE, nothing checked
    an <- env$range_keep_assay_qc(list(k27 = suppressWarnings(annotateSeq(k27, TxDb = txdb, verbose = FALSE))), qc)
    expect_identical(annoChecks(an$k27)$assay, qc$k27)
    expect_null(suppressMessages(env$range_assay_qc(list(k27 = k27), txdb, NULL)))
    expect_match(paste(env$range_assay_section(qc), collapse = " "), "no flag does not verify the assay", fixed = TRUE)
})
