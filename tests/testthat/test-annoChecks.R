library(epiSeeker)

context("test the checks of annotateSeq: dropped peaks and genome build")

skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

## n peaks of width 200 centred on TSSs of chr1-chr5 (a TF-like set), optionally shifted
tss_peaks <- function(n = 400, shift = 0) {
    tx <- GenomicFeatures::transcripts(txdb)
    tx <- tx[as.character(GenomicRanges::seqnames(tx)) %in% paste0("chr", 1:5)]
    tss <- unique(GenomicRanges::resize(tx, 1))[seq(1, by = 50, length.out = n)]
    GenomicRanges::GRanges(GenomicRanges::seqnames(tss),
                           IRanges::IRanges(GenomicRanges::start(tss) - 100 + shift, width = 200))
}

test_that("dropped peaks are counted per sequence with the reason", {
    pk <- suppressWarnings(c(tss_peaks(20), GenomicRanges::GRanges(c("chrEBV", "chrEBV", "1"),
                                                                 IRanges::IRanges(c(100, 500, 1e6), width = 100))))
    expect_warning(x <- annotateSeq(pk, TxDb = txdb, verbose = FALSE), "3 of 23 peaks were dropped")
    d <- annoChecks(x)$dropped
    expect_equal(d$n, 3)
    expect_equal(d$by_seqlevel$seqlevel, c("chrEBV", "1"))
    expect_equal(d$by_seqlevel$reason, c("not in TxDb", "named differently"))
    expect_warning(annotateSeq(pk, TxDb = txdb, verbose = FALSE), "seqlevelsStyle")
})

test_that("peaks at TSSs on the right build raise no warning and are recorded as fine", {
    expect_no_warning(x <- annotateSeq(tss_peaks(), TxDb = txdb, verbose = FALSE))
    b <- annoChecks(x)$build
    expect_false(b$suspect)
    expect_equal(b$n_out_of_bound, 0)
    expect_equal(annoChecks(x)$dropped$n, 0)
    expect_true(b$reliable)
    expect_equal(b$n_peaks, 400)
    expect_false(any(grepl("unreliable", capture.output(show(x)))))
})

test_that("a small set is recorded as unreliable for the build check, without a warning", {
    expect_no_warning(x <- annotateSeq(tss_peaks(50), TxDb = txdb, verbose = FALSE))
    b <- annoChecks(x)$build
    expect_false(b$suspect)
    expect_false(b$reliable)
    expect_equal(b$n_peaks, 50)
    expect_output(show(x), "with 50 peaks the genome-build check is unreliable")
})

test_that("coordinates beyond the chromosome ends of TxDb are flagged", {
    pk <- c(tss_peaks(), GenomicRanges::GRanges("chr1", IRanges::IRanges(249e6, width = 200)))   # chr1 of hg38 ends at 248,956,422
    expect_warning(x <- annotateSeq(pk, TxDb = txdb, verbose = FALSE), "may not be on the genome build")
    b <- annoChecks(x)$build
    expect_true(b$suspect)
    expect_equal(b$n_out_of_bound, 1)
    expect_output(show(x), "may not be on this genome build")
})

test_that("a genome recorded on the peaks that differs from TxDb already stops in GenomicRanges", {
    pk <- tss_peaks()
    ## only the sequences in use: with all ~700 of TxDb the message, which lists them, is cut at
    ## R's error-message limit (getOption("warning.length")) before "incompatible genomes"
    GenomeInfoDb::seqlevels(pk) <- GenomeInfoDb::seqlevelsInUse(pk)
    GenomeInfoDb::genome(pk) <- "hg19"
    expect_error(annotateSeq(pk, TxDb = txdb, verbose = FALSE), "genomes")
})

test_that("the TSS-proximity signal kept for calibration separates TSS peaks from shifted ones", {
    f <- epiSeeker:::getGene(txdb, by = "transcript")
    on <- epiSeeker:::.anno_check_build(tss_peaks(), f, txdb, tss = TRUE)
    off <- epiSeeker:::.anno_check_build(tss_peaks(shift = 30000), f, txdb, tss = TRUE)
    expect_gt(on$tss_ratio, 4)
    expect_lt(abs(log2(off$tss_ratio)), 1)
    expect_null(epiSeeker:::.anno_check_build(tss_peaks(), f, txdb)$tss_ratio)   # not computed by default
})

test_that("the shipped example files raise no build warning", {
    for (f in c(unlist(getSampleFiles()), system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))) {
        x <- annotateSeq(f, TxDb = txdb, verbose = FALSE)
        expect_false(annoChecks(x)$build$suspect, info = basename(f))
    }
})

test_that("objects made before the checks existed still work", {
    data(peakAnno, package = "epiSeeker")
    expect_null(annoChecks(peakAnno))
    expect_output(show(peakAnno), "peaks were annotated")
})
