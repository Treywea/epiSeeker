## Adapter from differential methylation callers to the range result classes.
##
## A DMR is an interval, so DMR results go to the range group: readDMR() turns the output of
## DSS (callDMR), methylKit (getMethylDiff, on bases or tiled regions) or dmrseq into a GRanges
## with the same columns whatever the caller -- effect, direction, pvalue, qvalue, n_CpG, caller
## -- and records where it came from in metadata(<GRanges>)$dmrSource, like annoSource and
## modSource. No new class or slot: annotateSeq() gives a csAnno whose anno keeps the columns
## and the metadata, and interpret() reads them.

#' Read DMRs from DSS, methylKit or dmrseq
#'
#' Turns the DMR (or DMC) table of a differential methylation caller into a `GRanges` ready for
#' [annotateSeq()], with columns that mean the same whatever the caller:
#'
#' `effect`: DSS `diff.Methy` (meanMethy1 - meanMethy2); methylKit `meth.diff / 100`
#' (treatment - control); dmrseq `meanDiff` when present (from `dmrseq::meanDiff()`), else `beta`
#' (not a methylation difference).
#'
#' `direction`: `"hyper"` / `"hypo"` from the sign of `effect` (DSS multi-factor tables, which
#' have no `diff.Methy`: from the sign of `areaStat`).
#'
#' `pvalue`, `qvalue`: methylKit `pvalue`, `qvalue`; dmrseq `pval`, `qval`; `NA` for DSS, whose
#' `callDMR()` gives none.
#'
#' `n_CpG`: DSS `nCG`; dmrseq `L`; `NA` for methylKit. `caller`: the caller.
#'
#' `direction` is `"hyper"` when methylation is higher in the first group of the caller's
#' comparison (DSS group1, methylKit treatment = 1, the non-reference level of dmrseq's
#' `testCovariate`) and `"hypo"` when lower; name the groups with `groups` so the provenance says
#' which. The other columns of the caller are kept unchanged. `metadata()$dmrSource` records the
#' caller, the input, what `effect` is, the comparison, the genome when known, and the numbers of
#' hyper and hypo regions. Nothing is filtered: DSS and methylKit tables are usually already
#' thresholded by the caller, dmrseq returns every candidate region (keep, e.g.,
#' `x[x$qvalue < 0.05]`).
#'
#' @param x the caller's result: a `data.frame` (DSS `callDMR()`; methylKit `getData()` of a
#'   `methylDiff`), a methylKit `methylDiff` object, a `GRanges` (dmrseq), or the path of one of
#'   these tables written with a header (tab- or comma-separated, optionally gzipped; e.g. by
#'   `write.table()`), with the caller's 1-based coordinates -- not a BED conversion.
#' @param caller `"auto"` (from the column names), `"DSS"`, `"methylKit"` or `"dmrseq"`.
#' @param groups optional names of the two groups compared, first the one where `"hyper"` means
#'   more methylated (DSS group1, methylKit treatment, dmrseq non-reference level), then the
#'   reference.
#' @param genome optional genome build (e.g. `"hg38"`), recorded in the provenance; taken from
#'   a methylKit object's assembly or the `GRanges` seqinfo when not given.
#' @return a `GRanges` with columns `effect`, `direction`, `pvalue`, `qvalue`, `n_CpG`,
#'   `caller` and the caller's other columns; `metadata()$dmrSource` as described.
#' @seealso [annotateSeq()], [interpret()]
#' @examples
#' dss <- data.frame(chr = c("chr1", "chr2"), start = c(1001, 5001), end = c(1500, 5300),
#'                   length = c(500, 300), nCG = c(12, 8), meanMethy1 = c(0.8, 0.2),
#'                   meanMethy2 = c(0.3, 0.6), diff.Methy = c(0.5, -0.4),
#'                   areaStat = c(60, -35))
#' dmr <- readDMR(dss, groups = c("tumor", "normal"))
#' dmr
#' S4Vectors::metadata(dmr)$dmrSource$n_direction
#' @export
readDMR <- function(x, caller = c("auto", "DSS", "methylKit", "dmrseq"), groups = NULL,
                    genome = NULL) {
    caller <- match.arg(caller)
    if (!is.null(groups) && (!is.character(groups) || length(groups) != 2)) {
        stop("`groups` must be two names: the group where \"hyper\" means more methylated, then the reference.")
    }
    input <- if (is.character(x) && length(x) == 1) normalizePath(x, mustWork = TRUE) else
        paste0("<", class(x)[1], " object>")
    extra <- list()
    if (is.character(x) && length(x) == 1) {
        df <- .dmr_read_table(x)
    } else if (is(x, "GRanges")) {
        if (is.null(genome)) {
            g <- unique(stats::na.omit(GenomeInfoDb::genome(x)))
            if (length(g) == 1) genome <- g
        }
        df <- as.data.frame(x)
    } else if (is(x, "methylDiff") || is(x, "methylDiffDB")) {
        sl <- function(n) if (methods::.hasSlot(x, n)) methods::slot(x, n) else NULL
        extra <- list(resolution = sl("resolution"), sample.ids = sl("sample.ids"),
                      treatment = sl("treatment"), context = sl("context"),
                      destranded = sl("destranded"))
        if (is.null(genome) && !is.null(sl("assembly"))) genome <- sl("assembly")
        ## a tabix-backed methylDiffDB holds no rows itself; methylKit (not a dependency) reads it
        if (is(x, "methylDiffDB")) stop("give the rows of a methylDiffDB: readDMR(methylKit::getData(x)).")
        df <- methods::as(x, "data.frame")
        if (caller == "auto") caller <- "methylKit"
    } else if (is.data.frame(x)) {
        df <- as.data.frame(x)
    } else {
        stop("`x` must be a data.frame, a methylDiff, a GRanges or a file path.")
    }
    if (!nrow(df)) stop("`x` has no rows.")
    if (caller == "auto") caller <- .dmr_detect(names(df))
    res <- switch(caller, DSS = .dmr_dss(df), methylKit = .dmr_methylkit(df), dmrseq = .dmr_dmrseq(df))
    chr <- intersect(c("chr", "seqnames", "chrom", "Chrom", "Chr", "chromosome"), names(df))[1]
    if (is.na(chr) || !all(c("start", "end") %in% tolower(names(df)))) {
        stop("cannot find the chromosome, start and end columns.")
    }
    st <- df[[names(df)[tolower(names(df)) == "start"][1]]]
    en <- df[[names(df)[tolower(names(df)) == "end"][1]]]
    strand <- if ("strand" %in% names(df)) as.character(df$strand) else rep("*", nrow(df))
    strand[!strand %in% c("+", "-", "*")] <- "*"
    eff <- res$effect
    dir_from <- if (all(is.na(eff)) && !is.null(res$direction_from)) res$direction_from else eff
    direction <- ifelse(is.na(dir_from) | dir_from == 0, NA_character_, ifelse(dir_from > 0, "hyper", "hypo"))
    used <- c(chr, names(df)[tolower(names(df)) %in% c("start", "end", "strand", "width")], res$used)
    keep <- df[, setdiff(names(df), used), drop = FALSE]
    gr <- GenomicRanges::GRanges(as.character(df[[chr]]),
                                 IRanges::IRanges(as.numeric(st), as.numeric(en)), strand = strand)
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        effect = as.numeric(eff), direction = direction,
        pvalue = as.numeric(res$pvalue), qvalue = as.numeric(res$qvalue),
        n_CpG = as.numeric(res$n_CpG), caller = caller, keep, check.names = FALSE)
    if (caller == "methylKit" && identical(extra$resolution, "base")) {
        message("methylKit resolution is \"base\": these are single differentially methylated CpGs (DMCs), not regions.")
    }
    nd <- table(factor(direction, levels = c("hyper", "hypo")), useNA = "ifany")
    S4Vectors::metadata(gr)$dmrSource <- c(list(
        caller = caller, input = input, effect = res$effect_desc,
        comparison = res$comparison,
        groups = if (is.null(groups)) NULL else stats::setNames(groups, c("hyper_in", "reference")),
        genome = if (is.null(genome)) NA_character_ else as.character(genome),
        columns = res$columns, n_regions = length(gr),
        n_direction = stats::setNames(as.integer(nd), ifelse(is.na(names(nd)), "none", names(nd))),
        median_width = stats::median(GenomicRanges::width(gr))
    ), extra[!vapply(extra, is.null, TRUE)])
    gr
}

.dmr_read_table <- function(file) {
    con <- if (grepl("\\.gz$", file)) gzfile(file) else file(file)
    first <- readLines(con, n = 1)
    close(con)
    sep <- if (grepl("\t", first)) "\t" else ","
    df <- utils::read.table(file, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                            check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[1] <- sub("^#", "", names(df)[1])
    ## write.csv() leaves an unnamed first column of row names
    if (identical(names(df)[1], "")) df <- df[, -1, drop = FALSE]
    df
}

.dmr_detect <- function(cols) {
    if (all(c("beta", "stat") %in% cols) && any(c("qval", "pval") %in% cols)) return("dmrseq")
    if ("meth.diff" %in% cols) return("methylKit")
    if (all(c("nCG", "areaStat") %in% cols)) return("DSS")
    stop("cannot tell the DMR caller from the columns (", paste(utils::head(cols, 12), collapse = ", "),
         "); give `caller`. Supported: DSS callDMR(), methylKit getMethylDiff(), dmrseq().")
}

.dmr_need <- function(df, cols, caller) {
    miss <- setdiff(cols, names(df))
    if (length(miss)) stop(caller, " table is missing column(s): ", paste(miss, collapse = ", "))
}

.dmr_dss <- function(df) {
    .dmr_need(df, c("nCG", "areaStat"), "DSS callDMR()")
    has_diff <- "diff.Methy" %in% names(df)
    na <- rep(NA_real_, nrow(df))
    list(effect = if (has_diff) df$diff.Methy else na,
         direction_from = df$areaStat,
         pvalue = na, qvalue = na, n_CpG = df$nCG,
         used = c(if (has_diff) "diff.Methy", "nCG", "length"),
         columns = c(effect = if (has_diff) "diff.Methy" else NA_character_, n_CpG = "nCG"),
         effect_desc = if (has_diff) "methylation difference, meanMethy1 - meanMethy2 (fraction)" else
             "none (multi-factor callDMR); direction from the sign of areaStat",
         comparison = "DSS group1 vs group2 (hyper = higher in group1)")
}

.dmr_methylkit <- function(df) {
    .dmr_need(df, c("meth.diff", "pvalue", "qvalue"), "methylKit getMethylDiff()")
    list(effect = df$meth.diff / 100, pvalue = df$pvalue, qvalue = df$qvalue,
         n_CpG = rep(NA_real_, nrow(df)), used = c("meth.diff", "pvalue", "qvalue"),
         columns = c(effect = "meth.diff / 100", pvalue = "pvalue", qvalue = "qvalue"),
         effect_desc = "methylation difference, treatment - control (fraction; meth.diff / 100)",
         comparison = "methylKit treatment = 1 vs treatment = 0 (hyper = higher in treatment 1)")
}

.dmr_dmrseq <- function(df) {
    .dmr_need(df, c("beta", "L"), "dmrseq()")
    md <- "meanDiff" %in% names(df)
    na <- rep(NA_real_, nrow(df))
    list(effect = if (md) df$meanDiff else df$beta,
         pvalue = if ("pval" %in% names(df)) df$pval else na,
         qvalue = if ("qval" %in% names(df)) df$qval else na,
         n_CpG = df$L, used = c(if (md) "meanDiff", "pval", "qval", "L"),
         columns = c(effect = if (md) "meanDiff" else "beta", pvalue = "pval", qvalue = "qval", n_CpG = "L"),
         effect_desc = if (md) "methylation difference from dmrseq::meanDiff() (fraction)" else
             "dmrseq beta: coefficient of the condition on the arcsine-transformed scale, not a methylation difference",
         comparison = "dmrseq testCovariate: non-reference vs reference level (hyper = higher in the non-reference level)")
}

## "DMRs: 1,234 from DSS; effect = ...; hyper in tumor vs normal: 800, hypo: 434". Empty for
## regions that do not come from readDMR().
.dmr_line <- function(gr, prefix = "DMRs: ") {
    src <- S4Vectors::metadata(gr)$dmrSource
    if (is.null(src)) return("")
    nd <- table(factor(gr$direction, levels = c("hyper", "hypo")))
    grp <- if (!is.null(src$groups)) paste0(" (", src$groups[1], " vs ", src$groups[2], ")") else ""
    pq <- if (all(is.na(gr$qvalue))) paste0("q-values: not given by ", src$caller) else
        paste0("q < 0.05: ", sum(gr$qvalue < 0.05, na.rm = TRUE), " of ", length(gr))
    paste0(prefix, format(length(gr), big.mark = ","), " from ", src$caller, "; effect = ", src$effect,
           "; ", src$comparison, grp, "; hyper ", nd[["hyper"]], ", hypo ", nd[["hypo"]], "; ", pq, "\n")
}
