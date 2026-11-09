## Assay QC: does an annotated peak set look like the assay the user declares? A declared
## H3K4me3 set with 10% of its peaks at promoters, or a "TF" set of 5 kb domains, is analysed
## and interpreted as if the label were right unless something says otherwise. The check
## compares three features of the set with the ranges seen in ENCODE peak sets of the declared
## assay; it never changes the annotation. Its result is kept in metadata(<anno>)$annoChecks
## when annotateSeq() is given `assay`, and is a line of the interpret() digest.
## Ranges: inst/extdata/assayQC_ranges.tsv, derived from the training half of a stratified
## sample of ENCODE GRCh38 peak files; derivation and held-out evaluation in epiSeeker-LLM
## checks/20261109_seqassayqc (DESIGN.md, RESULTS.md).

.assayqc_min_peaks <- 200
.assayqc_features <- c("promoter_1kb", "width_median", "gene_body")
.assayqc_marks <- "H[234][AB]?K[0-9]+(ME[123]|AC)"

## the features of a csAnno that the check uses, and a few recorded alongside
.assay_metrics <- function(x) {
    gr <- x@anno
    mc <- S4Vectors::mcols(gr)
    d <- mc$distanceToTSS
    w <- GenomicRanges::width(gr)
    a <- if ("annotation" %in% names(mc)) as.character(mc$annotation) else NULL
    q <- if (length(w)) stats::quantile(w, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
    frac <- function(z) if (is.null(a) || !length(a)) NA_real_ else mean(z)
    c(n_peaks = length(gr),
      promoter_1kb = if (length(d)) mean(abs(d) <= 1000, na.rm = TRUE) else NA_real_,
      promoter_anno = frac(grepl("^Promoter", a)),
      width_median = q[2], width_q25 = q[1], width_q75 = q[3],
      wide_5kb = if (length(w)) mean(w > 5000) else NA_real_,
      gene_body = frac(grepl("^(5' UTR|3' UTR|Exon|Intron)", a)),
      distal_intergenic = frac(a == "Distal Intergenic"))
}

## the stored ranges, read once per session
.assayqc_env <- new.env(parent = emptyenv())
.assayqc_ranges <- function() {
    if (is.null(.assayqc_env$ranges)) {
        f <- system.file("extdata", "assayQC_ranges.tsv", package = "epiSeeker")
        hdr <- readLines(f)
        .assayqc_env$source <- sub("^#\\s*", "", hdr[startsWith(hdr, "#")])
        .assayqc_env$ranges <- utils::read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
    }
    .assayqc_env$ranges
}

## declared assay -> class of the range table, with a note when the class is assumed
.assay_class <- function(assay, classes) {
    s <- trimws(as.character(assay))
    none <- function(why) list(class = NA_character_, note = why)
    m <- regmatches(toupper(s), regexpr(.assayqc_marks, toupper(s)))
    if (length(m)) {
        mark <- sub("ME([123])$", "me\\1", sub("AC$", "ac", m))
        if (mark %in% classes) return(list(class = mark, note = NULL))
        return(none(paste0("no ENCODE reference range for ", mark)))
    }
    if (grepl("ATAC", s, ignore.case = TRUE)) return(list(class = "ATAC", note = NULL))
    if (grepl("DNase|\\bDHS\\b", s, ignore.case = TRUE)) return(list(class = "DNase", note = NULL))
    if (grepl("^(TF|transcription factor)\\b", s, ignore.case = TRUE) || grepl("\\bTF\\b", s))
        return(list(class = "TF", note = NULL))
    if (grepl("input|IgG|control|mock", s, ignore.case = TRUE))
        return(none("a control is not an assay with a reference range"))
    if (grepl("^(H2A|H2B|H3|H4|macroH2A|histone)", s, ignore.case = TRUE))
        return(none(paste0("no ENCODE reference range for ", s)))
    tok <- trimws(gsub("ChIP-seq|ChIP|CUT&RUN|CUT&Tag|CUTandRUN|CUTandTag|-seq", "", s, ignore.case = TRUE))
    if (grepl("^[A-Za-z][A-Za-z0-9-]{0,14}$", tok))
        return(list(class = "TF", note = paste0("'", s, "' taken as a ChIP of a transcription factor or other chromatin protein")))
    none(paste0("'", s, "' does not name a transcription factor, histone mark, ATAC or DNase"))
}

## "12.3%" / "1,250 bp"
.assayqc_fmt <- function(feature, v) {
    if (feature == "width_median") paste0(format(round(v), big.mark = ",", trim = TRUE), " bp")
    else paste0(round(100 * v, 1), "%")
}
.assayqc_label <- c(promoter_1kb = "promoter fraction (<= 1 kb from a TSS)",
                    width_median = "median width",
                    gene_body = "gene-body fraction (UTR, exon or intron, outside the promoter)")

#' Does a peak set look like the assay it is said to be?
#'
#' Compares an annotated peak set with peak sets of the declared assay in ENCODE: the fraction
#' of peaks within 1 kb of a TSS, the median width and the fraction in gene bodies (UTR, exon
#' or intron, outside the promoter). A feature outside the range of the declared assay is a
#' flag, with the reason in words, e.g. a declared H3K4me3 set with 10% of its peaks at
#' promoters. The check is a rule, not a model: `interpret()` puts its result in the evidence
#' digest, and a set flagged here cannot be called consistent with its assay.
#'
#' The ranges (`inst/extdata/assayQC_ranges.tsv`, see `source` in the result) come from the
#' training half of a stratified sample of ENCODE GRCh38 peak files (TF ChIP-seq, eight
#' histone marks, ATAC-seq and DNase-seq, many cell types and tissues), widened by a margin
#' chosen so that at most 5% of training files are flagged under their own label (leave one
#' out; m = 2 standard deviations). On the held-out half (110 files of other biosamples or
#' targets) 3.6% of correctly labelled sets were flagged, and 59% of swaps between a narrow
#' assay (TF, H3K4me3, H3K9ac, H3K27ac, ATAC, DNase) and a broad mark (H3K27me3, H3K36me3,
#' H3K9me3, H3K79me2) were flagged: almost every swap with H3K9me3, most broad marks declared
#' as H3K4me3, H3K9ac or DNase, but **no broad mark declared as "TF"**, because the ENCODE peak
#' sets of transcription factors span nearly every promoter and gene-body fraction and the
#' ENCODE peaks of broad marks (MACS2 narrowPeak) are not wider than TF peaks. So a flag is
#' informative, but no flag does not confirm the label, least of all "TF". The ranges are from
#' ENCODE pipelines (MACS2 for histone marks, IDR for TF and ATAC, hotspot2 for DNase): peak
#' sets from other callers, e.g. broad-peak callers for H3K27me3, can differ in width.
#' Details: epiSeeker-LLM checks/20261109_seqassayqc.
#'
#' Declared assays are matched to the table by name: a histone mark in the string
#' ("H3K27me3 CUT&Tag"), "ATAC", "DNase"/"DHS", "TF", or a single other symbol, which is taken
#' as a transcription factor (e.g. "CTCF"; noted in the result). Not assessable (no flag):
#' fewer than 200 peaks, an assay without a range (e.g. H3K4me2, input), or no assay.
#'
#' @param x a `csAnno` object from [annotateSeq()]
#' @param assay the declared assay or target, e.g. "H3K4me3", "TF", "CTCF", "ATAC"
#' @return a list: `assay`, `class` (the row of the range table used, `NA` if none), `note`,
#'   `assessable`, `consistent` (`FALSE` when flagged, `TRUE` when not contradicted, `NA` when
#'   not assessable), `reasons`
#'   (one per feature outside its range, or why it is not assessable), `fits` (the classes
#'   whose ranges contain all three features), `metrics` (the features and a few more:
#'   `n_peaks`, `promoter_anno`, width quartiles, `wide_5kb`, `distal_intergenic`), `ranges`
#'   (the comparison, a data frame) and `source` (where the ranges come from).
#' @examples
#' data(peakAnno)
#' seqAssayQC(peakAnno, "H3K4me3")$consistent
#' seqAssayQC(peakAnno, "H3K36me3")$reasons
#' @seealso [annoChecks()], [annotateSeq()] (argument `assay`)
#' @export
seqAssayQC <- function(x, assay) {
    if (!is(x, "csAnno")) stop("`x` must be a csAnno object from annotateSeq().")
    if (missing(assay) || is.null(assay) || length(assay) != 1 || is.na(assay) || !nzchar(trimws(assay)))
        stop("`assay` must be one declared assay, e.g. \"H3K4me3\", \"TF\" or \"ATAC\".")
    .assayqc_decide(.assay_metrics(x), assay)
}

## the decision on the metrics of one peak set (also used by the evaluation on ENCODE files)
.assayqc_decide <- function(met, assay, rng = .assayqc_ranges()) {
    cl <- .assay_class(assay, unique(rng$class))
    out <- list(assay = as.character(assay), class = cl$class, note = cl$note, assessable = FALSE,
                consistent = NA, reasons = character(0), fits = character(0), metrics = met,
                ranges = NULL, source = .assayqc_env$source)
    inside <- function(r) {
        v <- met[r$feature]
        !is.na(v) & v >= r$lower & v <= r$upper
    }
    enough <- !is.na(met[["n_peaks"]]) && met[["n_peaks"]] >= .assayqc_min_peaks
    ## which classes the set fits, whatever was declared: a hint for the reader or the model
    if (enough) {
        out$fits <- names(Filter(isTRUE, lapply(split(rng, rng$class)[unique(rng$class)],
                                                function(r) all(inside(r)))))
    }
    if (is.na(cl$class)) {
        out$reasons <- cl$note
        return(out)
    }
    if (!enough) {
        out$reasons <- paste0("only ", met[["n_peaks"]], " peaks; the check needs at least ", .assayqc_min_peaks)
        return(out)
    }
    r <- rng[rng$class == cl$class, , drop = FALSE]
    r <- r[match(.assayqc_features, r$feature), , drop = FALSE]
    r$value <- unname(met[r$feature])
    r$inside <- inside(r)
    out$ranges <- r[, c("feature", "value", "lower", "upper", "train_min", "train_max", "n_train", "inside")]
    rownames(out$ranges) <- NULL
    out$assessable <- TRUE
    bad <- r[!is.na(r$value) & !r$inside, , drop = FALSE]
    out$consistent <- nrow(bad) == 0
    out$reasons <- vapply(seq_len(nrow(bad)), function(i) {
        f <- bad$feature[i]
        sprintf("%s %s is %s %s-%s, the range of ENCODE %s peak sets", .assayqc_label[[f]],
                .assayqc_fmt(f, bad$value[i]), if (bad$value[i] < bad$lower[i]) "below" else "above",
                .assayqc_fmt(f, bad$lower[i]), .assayqc_fmt(f, bad$upper[i]), cl$class)
    }, character(1))
    if (any(is.na(r$value)))
        out$note <- c(out$note, paste0("not available: ", paste(r$feature[is.na(r$value)], collapse = ", ")))
    out
}

## one line for show(), warnings and the interpret() digest
.assayqc_line <- function(qc) {
    if (is.null(qc)) return("")
    fits <- if (length(qc$fits)) paste0("; the profile fits ENCODE ", paste(qc$fits, collapse = ", ")) else
        "; the profile fits no class of the ENCODE table"
    if (isFALSE(qc$consistent)) {
        paste0("ASSAY QC INCONSISTENT - the regions do not look like the declared assay (", qc$assay, "): ",
               paste(qc$reasons, collapse = "; "), fits)
    } else if (isTRUE(qc$consistent)) {
        paste0("assay QC: not contradicted for ", qc$class, " (promoter fraction, median width and gene-body ",
               "fraction within the ENCODE range; this check misses many mislabels, so it does not confirm the label)", if (!is.null(qc$note)) paste0("; ", paste(qc$note, collapse = "; ")))
    } else {
        paste0("assay QC: not assessable (", paste(qc$reasons, collapse = "; "), ")")
    }
}
