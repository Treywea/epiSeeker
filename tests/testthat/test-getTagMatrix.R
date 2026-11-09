library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for getTagMatrix")

test_that("getTagMatrix function for single peak file", {
    data(demo_peak)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    # make window through txdb object
    tagMatrix <- getTagMatrix(demo_peak,
        type = "start_site", by = "gene",
        upstream = 3000, downstream = 3000,
        TxDb = txdb, weightCol = "V7"
    )

    # test tagmatrix record
    expect_equal(attr(tagMatrix, "type"), "start_site")
    expect_equal(attr(tagMatrix, "by"), "gene")

    # test tagMatrix upstream and downstream
    expect_equal(dim(tagMatrix)[2], 6001)
})

test_that("body tag matrix works with a sub-1kb flank extension", {
    ## issue #250 (ChIPseeker): upstream/downstream = 500 made ChIPseeker's
    ## binning code divide by zero.  epiSeeker bins the flanks separately, so
    ## this guards the same scenario here
    data(demo_peak)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    tagMatrix <- getTagMatrix(demo_peak,
        type = "body", by = "gene",
        upstream = 500, downstream = 500, nbin = 50,
        TxDb = txdb, weightCol = "V7",
        verbose = FALSE
    )

    expect_true(is.matrix(tagMatrix))
    expect_true(ncol(tagMatrix) > 50)
    expect_false(anyNA(tagMatrix))
    ## the flanks have round(50 * 0.15) = 8 bins each and carry values (a fractional count of
    ## bins, 7.5, gave flank columns of zeros before, and the 4 genes with peaks only in a
    ## flank were left out)
    expect_equal(ncol(tagMatrix), 66)
    expect_gt(sum(tagMatrix[, 1:8]), 0)
    expect_gt(sum(tagMatrix[, 59:66]), 0)
    expect_equal(attr(tagMatrix, "upstream"), 500)
    expect_equal(attr(tagMatrix, "downstream"), 500)
})

# test the getPromoters can run normally or not
test_that("getPromoters runs with default parameters", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    res <- getPromoters(txdb)

    expect_s4_class(res, "GRanges")
    expect_true(length(res) > 0)
})

make_issue236_case <- function(type = c("site", "body")) {
    type <- match.arg(type)
    peak_seqlevels <- c("chr1", "chr2", "chr10")
    peak_starts <- c(101L, 201L, 301L)

    if (type == "site") {
        peak <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                peak_seqlevels,
                levels = peak_seqlevels
            )),
            ranges = IRanges::IRanges(peak_starts, width = 5L),
            strand = "*"
        )
        S4Vectors::mcols(peak)$weight <- c(1, 2, 10)

        window_index <- c(1L, 3L, 2L)
        window_seqlevels <- peak_seqlevels[window_index]
        windows <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                window_seqlevels,
                levels = window_seqlevels
            )),
            ranges = IRanges::IRanges(
                peak_starts[window_index],
                width = 5L
            ),
            strand = "+"
        )
        attr(windows, "type") <- "start_site"
        attr(windows, "upstream") <- 2
        attr(windows, "downstream") <- 2
        attr(windows, "label") <- "center"
    } else {
        peak <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                peak_seqlevels,
                levels = peak_seqlevels
            )),
            ranges = IRanges::IRanges(peak_starts, width = 20L),
            strand = "*"
        )
        S4Vectors::mcols(peak)$weight <- c(1, 2, 10)

        window_index <- rep(c(1L, 3L, 2L), each = 2L)
        window_seqlevels <- peak_seqlevels[window_index]
        windows <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                window_seqlevels,
                levels = unique(window_seqlevels)
            )),
            ranges = IRanges::IRanges(
                peak_starts[window_index] + rep(c(0L, 10L), 3L),
                width = 10L
            ),
            strand = "+"
        )
        attr(windows, "type") <- "body"
        attr(windows, "upstream") <- 0
        attr(windows, "downstream") <- 0
        attr(windows, "label") <- c("SS", "TS")
    }

    attr(windows, "by") <- "issue236_demo"
    list(peak = peak, windows = windows)
}

test_that("getTagMatrix keeps site windows aligned across seqlevel orders", {
    case <- make_issue236_case("site")

    tag_matrix <- getTagMatrix(
        peak = case$peak,
        windows = case$windows,
        weightCol = "weight",
        verbose = FALSE,
        ignore_strand = TRUE
    )

    expect_equal(unname(rowSums(tag_matrix)), c(5, 50, 10))
})

test_that("getTagMatrix keeps body windows aligned across seqlevel orders", {
    case <- make_issue236_case("body")

    tag_matrix <- getTagMatrix(
        peak = case$peak,
        windows = case$windows,
        weightCol = "weight",
        nbin = 5,
        verbose = FALSE,
        ignore_strand = TRUE
    )

    expect_equal(
        unname(rowSums(tag_matrix)),
        c(5, 5, 50, 50, 10, 10)
    )
})

test_that("regions not longer than nbin are reported, not dropped silently (cases 2025)", {
    ## 20 long and 30 short regions, every one with a peak; nbin = 500
    st <- seq(1e5, by = 2e4, length.out = 50)
    wd <- c(rep(5000L, 20), rep(300L, 30))
    reg <- GenomicRanges::GRanges("chr1", IRanges::IRanges(st, width = wd), strand = "+")
    ## (plus one peak far beyond the last region: windows that pass the end of the last peak
    ## of a chromosome are left out before binning, an older and separate behaviour)
    pk <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(st + 50, 5e6), width = 100))
    w <- getBioRegion(TxDb = reg, by = "TE", type = "body", upstream = 0, downstream = 0)
    expect_warning(m <- suppressMessages(getTagMatrix(pk, windows = w, nbin = 500, verbose = FALSE)),
                   "30 of the 50 regions with peaks \\(60.0%\\) are not in the tag matrix")
    expect_equal(nrow(m), 20)
    s <- annoChecks(m)$short_regions
    expect_equal(s$n_dropped, 30)
    expect_equal(s$n_with_peaks, 50)
    expect_equal(s$fraction_dropped, 0.6)
    ## with rel() flanks the regions must be longer than all bins (nbin plus the flanks)
    w2 <- getBioRegion(TxDb = reg[1:20], by = "TE", type = "body", upstream = rel(0.2), downstream = rel(0.2))
    expect_warning(suppressMessages(getTagMatrix(pk, windows = w2, nbin = 5000, verbose = FALSE)),
                   "not longer than the 7000 bins")
    ## bp flanks, peaks over a flank and the body: the rows stay, with an empty body
    w3 <- getBioRegion(TxDb = reg[21:50], by = "TE", type = "body", upstream = 1000, downstream = 1000)
    pk3 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(st - 500, 5e6), width = c(rep(1500, 50), 100)))
    expect_warning(m3 <- suppressMessages(getTagMatrix(pk3, windows = w3, nbin = 500, verbose = FALSE)), "empty body")
    expect_equal(annoChecks(m3)$short_regions$n_empty_body, 30)
})

test_that("the binned body matrix equals the per-base implementation", {
    ## getTagMatrix_body_internal() sums the run-length coverage over the bins instead of
    ## expanding every window to one value per base (17-19 GB on the human genes)
    set.seed(1)
    pk <- GenomicRanges::GRanges(sample(c("chr1", "chr2"), 3000, TRUE),
                                 IRanges::IRanges(sample(1:5e6, 3000), width = sample(100:2000, 3000, TRUE)))
    win <- GenomicRanges::GRanges(sample(c("chr1", "chr2"), 300, TRUE),
                                  IRanges::IRanges(sample(1:4.5e6, 300), width = sample(c(50:30000), 300, TRUE)))
    for (cv in list(GenomicRanges::coverage(pk), GenomicRanges::coverage(pk, weight = stats::runif(3000)))) {
        for (nb in c(100L, 800L, 1120L)) {
            a <- epiSeeker:::getTagMatrix_body_internal(cv, win, nb, c("chr1", "chr2"))
            b <- epiSeeker:::.getTagMatrix_body_internal_dense(cv, win, nb, c("chr1", "chr2"))
            expect_equal(a, b, tolerance = 1e-12)
        }
    }
})
