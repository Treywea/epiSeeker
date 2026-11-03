library(epiSeeker)

context("test result structures: annoSource, peakOverlapResult, csAnnoList")

test_that("annotateSeq records where the annotation came from", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))[1:20]
    x <- annotateSeq(peak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                     annoDb = "org.Hs.eg.db", verbose = FALSE)
    src <- annoSource(x)
    expect_equal(src$genome, "hg38")
    expect_match(src$TxDb$source, "UCSC")
    expect_equal(src$TxDb$organism, "Homo sapiens")
    expect_equal(src$annoDb$name, "org.Hs.eg.db")
    expect_equal(src$annoDb$version, as.character(utils::packageVersion("org.Hs.eg.db")))
    expect_output(show(x), "Genome: hg38")
    ## combine_csAnno keeps the source of the first object
    expect_equal(annoSource(combine_csAnno(x, x))$TxDb, src$TxDb)
    ## no annoDb
    y <- annotateSeq(peak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                     verbose = FALSE)
    expect_null(annoSource(y)$annoDb)
})

test_that("a csAnno made before annoSource existed still works", {
    data(peakAnno, package = "epiSeeker")
    expect_false(methods::.hasSlot(peakAnno, "annoSource"))
    src <- annoSource(peakAnno)
    expect_true(is.character(src$genome))
    expect_null(src$TxDb)
    expect_output(show(peakAnno), "peaks were annotated")
})

test_that("enrichPeakOverlap adds effect sizes and stays a data.frame", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1001, 5001), width = 100))
    t1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1051, width = 100))      # half of one query peak
    t2 <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1001, width = 100))      # no overlap
    res <- enrichPeakOverlap(q, list(a = t1, b = t2),
                             TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                             nShuffle = 0, mc.cores = 1, verbose = FALSE)
    expect_s3_class(res, "peakOverlapResult")
    expect_s3_class(res, "data.frame")
    expect_equal(res$jaccard, c(50 / 250, 0))     # 50 bp shared; union 200 + 50
    expect_equal(res$frac_query, c(0.5, 0))
    expect_equal(res$frac_target, c(1, 0))
    expect_equal(attr(res, "params")$nShuffle, 0)
    expect_equal(attr(res, "params")$genome, "hg38")
    ## data.frame code keeps working
    expect_equal(nrow(res[res$N_OL > 0, ]), 1)
    expect_true(is.data.frame(as.data.frame(res)))
    ## pool = FALSE goes through the per-target path and is classed as well
    res2 <- enrichPeakOverlap(q, list(a = t1, b = t2),
                              TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                              nShuffle = 0, pool = FALSE, mc.cores = 1, verbose = FALSE)
    expect_s3_class(res2, "peakOverlapResult")
    expect_equal(res2$jaccard, res$jaccard)
    expect_false(attr(res2, "params")$pool)
})

test_that("csAnnoList checks its elements and works wherever a list does", {
    data(peakAnno, package = "epiSeeker")
    x <- csAnnoList(a = peakAnno, b = peakAnno)
    expect_s4_class(x, "csAnnoList")
    expect_true(validObject(x))
    expect_equal(names(x), c("a", "b"))
    expect_true(is.list(x))
    expect_true(is(x, "list"))
    expect_s4_class(x[["a"]], "csAnno")
    expect_identical(csAnnoList(list(a = peakAnno, b = peakAnno)), x)
    expect_s3_class(plotAnnoBar(x), "ggplot")
    expect_s3_class(plotDistToTSS(x), "ggplot")
    expect_equal(length(lapply(x, annoSource)), 2)
    expect_error(csAnnoList(peakAnno, peakAnno), "unique name")
    expect_error(csAnnoList(a = peakAnno, b = 1), "csAnno")
})
