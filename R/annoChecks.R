## Checks for silent errors of annotateSeq(): peaks on another genome build than TxDb, and
## peaks dropped because TxDb has no feature on their sequence. The checks never change the
## annotation; what they find is reported (warning) and kept in metadata(<anno>)$annoChecks,
## next to annoSource, so that interpret() and analysis scripts can read it.
## Calibration: epiSeeker-LLM checks/20261105_correctness-checks (DESIGN.md, RESULTS.md).

## shifts (bp) of the local background of the TSS-proximity check: far enough to break any
## precise relation to TSSs, near enough to keep the gene density of the region
.tss_shifts <- c(-2e5, -1e5, -5e4, 5e4, 1e5, 2e5)
.tss_window <- 1000
.tss_min_peaks <- 200
## |log2(observed / shifted)| below this would mean "no closer to TSSs than chance". Calibrated
## on 2026-11-05 (116 files, 246 file x TxDb pairs): every cutoff that raises no false alarm on
## correctly matched data detected no mismatch that the out-of-bound test missed, while the
## distal factors and broad marks of matched data reach 0.11-0.2. So it is not used for the
## warning (NA), and annotateSeq() does not compute it; the code is kept for re-calibration.
.tss_log2_cutoff <- NA_real_
## Below this many peaks the build check is recorded as unreliable: only peaks near chromosome
## ends can fall outside them, and the four mismatched sets the calibration missed had 50-160
## peaks. Recorded in annoChecks() and shown by show(), not warned, so that small example and
## test sets do not warn on every call.
.build_min_peaks <- 200

#' What the checks of annotateSeq found
#'
#' `annotateSeq()` checks two silent errors and keeps what it finds with the annotation:
#' \itemize{
#'   \item `build`: whether the peaks look like coordinates of another genome build than
#'     `TxDb`: peaks extend beyond the chromosome ends of `TxDb` (`n_out_of_bound`, with the
#'     chromosomes in `out_of_bound_seqlevels`). `suspect` is then `TRUE` and `annotateSeq()`
#'     warns. The genomes recorded on the peaks and on `TxDb` are kept (`genome_peaks`,
#'     `genome_TxDb`); peaks read from bed files usually record none, and peaks recorded with
#'     another genome than `TxDb` already stop with an error in GenomicRanges. On 116 ENCODE
#'     and example files annotated with their own build and with another one, the
#'     out-of-bound test raised no false alarm and found 132
#'     of 136 mismatched pairs; it can miss small sets (the four it missed had 50-150 peaks),
#'     because only peaks near chromosome ends can fall outside them. `n_peaks` is the number
#'     of peaks checked and `reliable` is `FALSE` below 200 of them; `show()` then says that
#'     the build check is unreliable for this set (no warning is raised for that alone).
#'   \item `dropped`: peaks that were not annotated because `TxDb` has no feature on their
#'     sequence, per sequence and reason (`not in TxDb`, `named differently`, `no feature`).
#' }
#' Objects made before these checks existed return `NULL`.
#'
#' @param x a `csAnno` object
#' @return a list with elements `build` and `dropped`, or `NULL`
#' @examples
#' data(peakAnno)
#' annoChecks(peakAnno)
#' @export
annoChecks <- function(x) {
    if (!is(x, "csAnno")) stop("`x` must be a csAnno object.")
    S4Vectors::metadata(x@anno)$annoChecks
}

## the transcription start sites of the features annotateSeq() uses
.feature_tss <- function(features) {
    tss <- GenomicRanges::resize(features, width = 1, fix = "start")
    GenomicRanges::reduce(GenomicRanges::granges(tss), ignore.strand = TRUE)
}

## fraction of positions within each window of the nearest TSS; NA when nothing to measure
.near_tss <- function(pos, tss, windows) {
    if (!length(pos)) return(rep(NA_real_, length(windows)))
    d <- GenomicRanges::distanceToNearest(pos, tss, ignore.strand = TRUE)
    dist <- rep(Inf, length(pos))
    dist[S4Vectors::queryHits(d)] <- S4Vectors::mcols(d)$distance
    vapply(windows, function(w) mean(dist <= w), numeric(1))
}

## All the build signals, without a decision: used by annotateSeq() and by the calibration.
.anno_check_build <- function(peak.gr, features, TxDb, tss = FALSE, windows = c(500, .tss_window, 3000)) {
    gp <- unique(stats::na.omit(GenomeInfoDb::genome(peak.gr)))
    gt <- tryCatch(unique(stats::na.omit(GenomeInfoDb::genome(TxDb))), error = function(e) character(0))
    ## recorded only: peaks tagged with another genome than TxDb never get here, GenomicRanges
    ## stops the nearest-feature step first ("sequences ... have incompatible genomes")
    out <- list(genome_peaks = if (length(gp)) paste(gp, collapse = ",") else NA_character_,
                genome_TxDb = if (length(gt)) paste(gt, collapse = ",") else NA_character_)

    ## S2: coordinates beyond the chromosome ends of TxDb
    sl <- tryCatch(GenomeInfoDb::seqlengths(TxDb), error = function(e) NULL)
    chr <- as.character(GenomicRanges::seqnames(peak.gr))
    len <- if (is.null(sl)) rep(NA_real_, length(chr)) else unname(sl[chr])
    oob <- !is.na(len) & GenomicRanges::end(peak.gr) > len
    out$n_out_of_bound <- sum(oob)
    out$out_of_bound_seqlevels <- names(sort(table(chr[oob]), decreasing = TRUE))

    if (!tss) return(out)

    ## S3: TSS proximity against the same peaks shifted locally (calibration only, see above)
    tss <- .feature_tss(features)
    keep <- chr %in% as.character(GenomeInfoDb::seqlevelsInUse(tss)) & !oob
    mid <- GenomicRanges::resize(GenomicRanges::granges(peak.gr[keep]), width = 1, fix = "center")
    out$n_tss_peaks <- length(mid)
    out$tss_window <- .tss_window
    if (length(mid) < .tss_min_peaks) {
        out$tss_ratio <- NA_real_
        out$tss_fraction <- out$tss_fraction_shifted <- NA_real_
        out$tss_ratio_windows <- stats::setNames(rep(NA_real_, length(windows)), windows)
        return(out)
    }
    obs <- .near_tss(mid, tss, windows)
    ## shifted copies stay on their chromosome: within its length, or within the span of its
    ## TSSs when TxDb has no lengths
    lim <- if (!is.null(sl) && !all(is.na(sl))) sl else
        vapply(split(GenomicRanges::end(tss), as.character(GenomicRanges::seqnames(tss))), max, numeric(1))
    mchr <- as.character(GenomicRanges::seqnames(mid))
    hits <- n <- numeric(length(windows))
    for (s in .tss_shifts) {
        p <- GenomicRanges::start(mid) + s
        ok <- p >= 1 & !is.na(lim[mchr]) & p <= lim[mchr]
        if (!any(ok)) next
        sh <- GenomicRanges::GRanges(mchr[ok], IRanges::IRanges(p[ok], width = 1))
        hits <- hits + .near_tss(sh, tss, windows) * sum(ok)
        n <- n + sum(ok)
    }
    shifted <- hits / n
    ratio <- obs / shifted
    names(ratio) <- windows
    w <- match(.tss_window, windows)
    out$tss_fraction <- obs[w]
    out$tss_fraction_shifted <- shifted[w]
    out$tss_ratio <- unname(ratio[w])
    out$tss_ratio_windows <- ratio
    out
}

## the decision on the build signals
.build_suspect <- function(b, cutoff = .tss_log2_cutoff) {
    tss_flat <- !is.na(cutoff) && !is.null(b$tss_ratio) && !is.na(b$tss_ratio) && b$tss_ratio > 0 &&
        abs(log2(b$tss_ratio)) < cutoff
    reasons <- c(
        if (b$n_out_of_bound > 0) paste0(b$n_out_of_bound, " peak(s) extend beyond the chromosome ends of TxDb (",
                                         paste(utils::head(b$out_of_bound_seqlevels, 3), collapse = ", "), ")"),
        if (tss_flat) sprintf(paste("the peaks are no closer to TSSs than the same peaks shifted by 50-200 kb",
                                    "(%.1f%% vs %.1f%% within 1 kb)"),
                              100 * b$tss_fraction, 100 * b$tss_fraction_shifted)
    )
    list(suspect = length(reasons) > 0, reasons = reasons)
}

## Peaks the nearest-feature step dropped, by sequence and reason.
.anno_check_dropped <- function(before, after, TxDb, features) {
    chr <- as.character(GenomicRanges::seqnames(before))
    b <- table(chr)
    a <- table(as.character(GenomicRanges::seqnames(after)))
    tab <- b - as.integer(ifelse(is.na(a[names(b)]), 0, a[names(b)]))
    tab <- sort(tab[tab > 0], decreasing = TRUE)
    if (!length(tab)) return(list(n = 0L, n_total = length(chr), by_seqlevel = NULL))
    txdb_sl <- tryCatch(GenomeInfoDb::seqlevels(TxDb), error = function(e) character(0))
    feat_sl <- as.character(GenomeInfoDb::seqlevelsInUse(features))
    reason <- vapply(names(tab), function(s) {
        if (s %in% txdb_sl) return("no feature")
        alt <- if (startsWith(s, "chr")) sub("^chr", "", s) else paste0("chr", s)
        if (alt %in% c(txdb_sl, feat_sl)) "named differently" else "not in TxDb"
    }, character(1))
    list(n = sum(tab), n_total = length(chr),
         by_seqlevel = data.frame(seqlevel = names(tab), n = as.integer(tab), reason = unname(reason),
                                  stringsAsFactors = FALSE))
}

.dropped_message <- function(d) {
    bs <- d$by_seqlevel
    parts <- vapply(split(bs, bs$reason)[unique(bs$reason)], function(z) {
        what <- switch(z$reason[1],
            "named differently" = "named differently from TxDb (e.g. '1' vs 'chr1'; align with GenomeInfoDb::seqlevelsStyle())",
            "no feature" = "on sequences of TxDb without any gene",
            "not in TxDb" = "on sequences TxDb does not have (contigs, scaffolds, chrEBV, ...)")
        shown <- utils::head(z, 3)
        paste0(sum(z$n), " ", what, ": ", paste0(shown$seqlevel, " (", shown$n, ")", collapse = ", "),
               if (nrow(z) > 3) ", ..." else "")
    }, character(1))
    paste0(d$n, " of ", d$n_total, " peaks were dropped, since no feature of 'TxDb' can be found for them: ",
           paste(parts, collapse = "; "), ".")
}
