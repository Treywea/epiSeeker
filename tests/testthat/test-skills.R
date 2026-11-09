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

test_that("the four skills are installed with SKILL.md, a script and a worked example", {
    expect_true(all(c("range-tf", "base-methylation", "interaction-loops", "range-singlecell") %in% skill_names))
    expect_true(file.exists(file.path(skills_dir, "shared", "skill-helpers.R")))
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
    helpers <- file.path(skills_dir, "shared", "skill-helpers.R")
    own <- assigned(tokens(helpers))
    expect_true(all(c("skill_check", "skill_interpret", "skill_finish") %in% own))
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
