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

## custom regions (TE-like) on chr1-chr5 in the coordinates of another build: the chromosome
## lengths of T2T-CHM13, whose chr1 is 569 kb shorter than hg38's (the 2025 case study)
chm13_len <- c(chr1 = 248387328, chr2 = 242696752, chr3 = 201105948, chr4 = 193574945, chr5 = 182045439)
te_regions <- function(with_lengths = TRUE) {
    set.seed(3)
    chr <- sample(names(chm13_len), 3000, TRUE)
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(round(stats::runif(3000, 1e5, 1.8e8)), width = 2000),
                                 strand = sample(c("+", "-"), 3000, TRUE))
    if (with_lengths) GenomeInfoDb::seqlengths(gr) <- chm13_len[GenomeInfoDb::seqlevels(gr)]
    gr
}
## hg38 peaks: TSS peaks plus 5 near the end of hg38 chr1, beyond the end of CHM13 chr1
hg38_peaks <- function() c(tss_peaks(), GenomicRanges::GRanges("chr1", IRanges::IRanges(248.5e6 + 1e4 * (1:5), width = 300)))

test_that("a GRanges annotation with chromosome lengths gets the build check (cases 2025)", {
    expect_warning(x <- annotateSeq(hg38_peaks(), TxDb = te_regions(), verbose = FALSE),
                   "may not be on the genome build")
    b <- annoChecks(x)$build
    pk <- hg38_peaks()
    n_oob <- sum(GenomicRanges::end(pk) > chm13_len[as.character(GenomicRanges::seqnames(pk))])   # the 5 + TSS peaks near chr1 end
    expect_gte(n_oob, 5)
    expect_true(b$suspect)
    expect_equal(b$n_out_of_bound, n_oob)
    expect_equal(b$annotation, "GRanges")
    expect_true(b$seqlengths_known)
    expect_output(show(x), "may not be on this genome build")
    expect_match(epiSeeker:::.checks_line(x), "GENOME BUILD SUSPECT")
    ## the same peaks against the regions with hg38 lengths: no warning
    te38 <- te_regions(FALSE)
    GenomeInfoDb::seqlengths(te38) <- GenomeInfoDb::seqlengths(txdb)[GenomeInfoDb::seqlevels(te38)]
    expect_no_warning(y <- annotateSeq(hg38_peaks(), TxDb = te38, verbose = FALSE))
    expect_false(annoChecks(y)$build$suspect)
})

test_that("a GRanges annotation without chromosome lengths is recorded as not checked", {
    expect_message(x <- annotateSeq(hg38_peaks(), TxDb = te_regions(FALSE), verbose = FALSE),
                   "no chromosome lengths")
    b <- annoChecks(x)$build
    expect_false(b$suspect)
    expect_false(b$seqlengths_known)
    expect_false(b$reliable)
    expect_match(b$note, "could not be checked")
    expect_output(show(x), "genome build was not checked")
    expect_match(epiSeeker:::.checks_line(x), "NOT CHECKED")
})

test_that("chromosome lengths recorded on the peaks that differ from the annotation are evidence", {
    pk <- hg38_peaks()
    GenomeInfoDb::seqlevels(pk) <- names(chm13_len)
    GenomeInfoDb::seqlengths(pk) <- GenomeInfoDb::seqlengths(txdb)[names(chm13_len)]
    b <- epiSeeker:::.build_check(pk, te_regions(), ref_label = "the regions")
    expect_true(b$suspect)
    expect_setequal(b$seqlength_conflicts, names(chm13_len))
})

test_that("getTagMatrix checks the build of the peaks against the regions (cases 2025)", {
    w <- getBioRegion(TxDb = te_regions(), by = "TE", type = "body", upstream = rel(0.2), downstream = rel(0.2))
    expect_warning(m <- getTagMatrix(hg38_peaks(), windows = w, nbin = 100, verbose = FALSE),
                   "may not be on the genome build of the regions")
    expect_true(annoChecks(m)$build$suspect)
    pk <- hg38_peaks()
    expect_equal(annoChecks(m)$build$n_out_of_bound,
                 sum(GenomicRanges::end(pk) > chm13_len[as.character(GenomicRanges::seqnames(pk))]))
    ## a list of windows: one check per set, named after it
    wg <- getBioRegion(TxDb = txdb, by = "gene", type = "body", upstream = rel(0.2), downstream = rel(0.2))
    expect_warning(ml <- getTagMatrix(hg38_peaks(), windows = list(gene = wg, TE = w), nbin = 100, verbose = FALSE),
                   "the regions 'TE'")
    ck <- annoChecks(ml)
    expect_false(ck$gene$build$suspect)
    expect_true(ck$TE$build$suspect)
    ## without lengths: a note, recorded as not checked
    w0 <- getBioRegion(TxDb = te_regions(FALSE), by = "TE", type = "body", upstream = rel(0.2), downstream = rel(0.2))
    expect_message(m0 <- getTagMatrix(hg38_peaks(), windows = w0, nbin = 100, verbose = FALSE), "no chromosome lengths")
    expect_false(annoChecks(m0)$build$seqlengths_known)
    ## the TxDb path: hg38 peaks and hg38 TxDb raise nothing
    expect_no_warning(m1 <- getTagMatrix(tss_peaks(), TxDb = txdb, type = "start_site", by = "gene",
                                         upstream = 1000, downstream = 1000, verbose = FALSE))
    expect_false(annoChecks(m1)$build$suspect)
    expect_true(annoChecks(m1)$build$reliable)
})
