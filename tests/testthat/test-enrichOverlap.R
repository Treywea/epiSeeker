library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for enrichAnnoOverlap ")

test_that("enrichAnnoOverlap works with example peakfile", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    peakfile <- system.file("extdata", "demo_peak.txt", package = "epiSeeker")

    res <- enrichAnnoOverlap(peakfile, peakfile, txdb)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("qSample", "tSample", "qLen", "tLen", "N_OL", "pvalue", "p.adjust") %in% names(res)))
    expect_equal(nrow(res), 1)
})


test_that("enrichPeakOverlap works with GRanges and file input", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    peakfile <- system.file("extdata", "demo_peak.txt", package = "epiSeeker")
    gr <- readPeakFile(peakfile)[1:10]

    # reduce shuffle = 5 to speed up test
    res <- enrichPeakOverlap(gr, peakfile, txdb, nShuffle = 5, mc.cores = 1, verbose = FALSE)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("qSample", "tSample", "qLen", "tLen", "N_OL", "pvalue", "p.adjust") %in% names(res)))
    expect_equal(nrow(res), 1)
})


test_that("enrichPeakOverlap accepts a single GRanges as target", {
    ## issue #84 of ChIPseeker: a bare GRanges was not wrapped into a list, so
    ## the permutation test failed with
    ## "GRanges objects don't support [[, as.list(), lapply()"
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    set.seed(1)
    q <- GRanges("chr1", IRanges(sample(1000000:1400000, 5), width = 300))
    t <- GRanges("chr1", IRanges(sample(1000000:1400000, 5), width = 300))

    res <- enrichPeakOverlap(q, t, TxDb = txdb, nShuffle = 20,
                             mc.cores = 1, verbose = FALSE)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 1)
    expect_equal(res$qLen, 5)
    expect_equal(res$tLen, 5)
    expect_true(res$N_OL >= 0 && res$N_OL <= 5)
})


test_that("the overlap count is symmetric while the tested ratio is not", {
    ## issue #84: N_OL is direction free, the p-value normalises by the target
    ## size and shuffles the target, so swapping the arguments can change it
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    set.seed(11)
    hot <- sample(1000000:1400000, 3)
    A <- GRanges("chr1", IRanges(c(sample(hot, 2),
                                   sample(2000000:2400000, 18)), width = 300))
    B <- GRanges("chr1", IRanges(c(sample(hot, 3),
                                   sample(3000000:3400000, 97)), width = 300))

    set.seed(1)
    ab <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)
    set.seed(1)
    ba <- enrichPeakOverlap(B, A, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)

    ## the number of overlapping peaks is direction free ...
    expect_equal(ab$N_OL, ba$N_OL)
    expect_equal(ab$N_OL, length(intersect(A, B)))
    expect_equal(ab$qLen, ba$tLen)
    expect_equal(ab$tLen, ba$qLen)

    ## ... while the ratio the p-value is based on is target normalised, so it
    ## changes when the two peak sets differ in size
    expect_false(isTRUE(all.equal(ab$N_OL / ab$tLen, ba$N_OL / ba$tLen)))

    expect_true(all(c(ab$pvalue, ba$pvalue) > 0))
    expect_true(all(c(ab$pvalue, ba$pvalue) <= 1))
})


test_that("symmetric=TRUE gives the same p-value in both directions", {
    ## issue #84 of ChIPseeker: the opt-in symmetric p-value combines both
    ## directions as min(1, 2*min(p, p_rev))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    set.seed(11)
    hot <- sample(1000000:1400000, 3)
    A <- GRanges("chr1", IRanges(c(sample(hot, 2),
                                   sample(2000000:2400000, 18)), width = 300))
    B <- GRanges("chr1", IRanges(c(sample(hot, 3),
                                   sample(3000000:3400000, 97)), width = 300))

    set.seed(1)
    d.ab <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                              verbose = FALSE
    )
    set.seed(1)
    s.ab <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                              verbose = FALSE, symmetric = TRUE
    )
    set.seed(1)
    s.ba <- enrichPeakOverlap(B, A, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                              verbose = FALSE, symmetric = TRUE
    )

    ## direction free ...
    expect_equal(s.ab$pvalue, s.ba$pvalue)
    expect_equal(s.ab$N_OL, s.ba$N_OL)

    ## ... and consistent with the one-sided test it is built from
    set.seed(1)
    p.rev <- epiSeeker:::enrichOverlap.peak.mirrored(A, list(B), txdb, 100,
                                                     mc.cores = 1, verbose = FALSE
    )$pvalue
    expect_equal(s.ab$pvalue, min(1, 2 * min(d.ab$pvalue, p.rev)))
    expect_true(s.ab$pvalue >= d.ab$pvalue)
})


test_that("shuffle", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p <- GRanges(
        seqnames = c("chr1", "chr3"),
        ranges = IRanges(
            start = c(1, 100),
            end = c(50, 130)
        )
    )

    res <- shuffle(p, TxDb = txdb)

    expect_s4_class(res, "GRanges")
    expect_equal(length(res), 2)
    expect_true(all(width(res) == width(p)))
})


test_that("enrichOverlap.peak.internal", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    query <- GRanges("chr1", IRanges(c(100, 200), c(150, 250)))
    target <- list(
        GRanges("chr1", IRanges(120, 160)),
        GRanges("chr1", IRanges(1000, 1100))
    )

    res <- epiSeeker:::enrichOverlap.peak.internal(
        query.gr = query,
        target.gr = target,
        TxDb = txdb,
        nShuffle = 5,
        mc.cores = 1,
        verbose = FALSE
    )

    expect_true(is.list(res))
    expect_true(all(c("pvalue", "overlap") %in% names(res)))
    expect_length(res$pvalue, length(target))
})


test_that("shuffle() takes chromosome lengths by name, not by position", {
    ## seqlevels chr1, chr2, chr10 (not alphabetical): the peaks of chr10 were placed
    ## along the length of chr2 and could fall off the end of chr10
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    len <- GenomeInfoDb::seqlengths(txdb)[c("chr1", "chr2", "chr10")]
    p <- GRanges(c(rep("chr1", 5), rep("chr2", 5), rep("chr10", 200)),
                 IRanges(1000, width = 100))
    GenomeInfoDb::seqlevels(p) <- c("chr1", "chr2", "chr10")
    set.seed(3)
    res <- shuffle(p, TxDb = txdb)
    on10 <- as.character(seqnames(res)) == "chr10"
    expect_equal(sum(on10), 200)
    expect_true(all(start(res)[on10] <= len[["chr10"]]))
})

test_that("pool = FALSE names every target set and adjusts across them", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    set.seed(5)
    q <- GRanges("chr1", IRanges(sample(1e6:5e6, 50), width = 300))
    tl <- list(a = q[1:20], b = GRanges("chr1", IRanges(sample(1e6:5e6, 30), width = 300)),
               c = GRanges("chr1", IRanges(sample(2e7:3e7, 30), width = 300)))
    set.seed(1)
    res <- enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 50, pool = FALSE, mc.cores = 1, verbose = FALSE)
    expect_equal(res$tSample, c("a", "b", "c"))
    expect_equal(res$p.adjust, p.adjust(res$pvalue, method = "BH"))
    ## each p-value is that of the target set on its own
    set.seed(1)
    one <- enrichPeakOverlap(q, tl["a"], TxDb = txdb, nShuffle = 50, mc.cores = 1, verbose = FALSE)
    expect_equal(res$pvalue[1], one$pvalue)
    expect_equal(attr(res, "params")$null, "per target set, shuffled over whole chromosomes")
    expect_true(all(c("expected", "fold_enrichment") %in% names(res)))
})

test_that("background shuffles target peaks within the background, one null per target set", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    bg <- GRanges(c("chr1", "chr1", "chr2"), IRanges(c(1e6, 4e6, 1e6), width = c(2e6, 1e6, 2e6)))
    bgi <- epiSeeker:::.prepare_background(bg)
    set.seed(2)
    pk <- GRanges(c(rep("chr1", 40), rep("chr2", 10), rep("chr3", 5)), IRanges(1, width = 300))
    s <- epiSeeker:::.shuffle_in_background(pk, bgi)
    expect_true(all(IRanges::overlapsAny(s, bg)))
    expect_equal(width(s), width(pk))
    expect_true(sum(as.character(seqnames(s)) == "chr1") >= 40)   # chr3 peaks go anywhere in bg

    ## a target set built from half of the query peaks scores far above unrelated ones
    q <- epiSeeker:::.shuffle_in_background(GRanges("chr1", IRanges(rep(1, 200), width = 300)), bgi)
    tl <- list(related = c(q[1:100], epiSeeker:::.shuffle_in_background(q[1:100], bgi)),
               unrelated = epiSeeker:::.shuffle_in_background(q, bgi))
    set.seed(1)
    res <- enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 50, mc.cores = 1, verbose = FALSE,
                             background = bg)
    expect_equal(res$tSample, c("related", "unrelated"))
    expect_equal(res$pvalue[1], 1 / 51)
    expect_gt(res$fold_enrichment[1], 2 * res$fold_enrichment[2])
    expect_equal(attr(res, "params")$null, "per target set, shuffled within background")
    ## pool is ignored with a background, symmetric works
    set.seed(1)
    res2 <- enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 50, mc.cores = 1, verbose = FALSE,
                              background = bg, pool = FALSE, symmetric = TRUE)
    expect_equal(res2$tSample, c("related", "unrelated"))
    expect_true(all(res2$pvalue > 0 & res2$pvalue <= 1))
    expect_error(enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 5, mc.cores = 1, verbose = FALSE,
                                   background = GRanges()), "no regions")
})

test_that("symmetric = TRUE is not significant when the two sets do not overlap", {
    ## the mirrored direction shuffled one query peak per permutation and compared its
    ## 0/1 overlap with the observed ratio: interleaved sets with no overlap got p = 0.02
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    st <- seq(1e6, by = 2e4, length.out = 100)
    A <- GRanges("chr1", IRanges(st, width = 500))
    B <- GRanges("chr1", IRanges(st + 1e4, width = 500))
    expect_equal(length(intersect(A, B)), 0)
    set.seed(1)
    res <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1, verbose = FALSE,
                             symmetric = TRUE)
    expect_equal(res$N_OL, 0)
    expect_equal(res$pvalue, 1)
    set.seed(1)
    p.rev <- epiSeeker:::enrichOverlap.peak.mirrored(A, list(B), txdb, 100, mc.cores = 1,
                                                     verbose = FALSE)$pvalue
    expect_equal(p.rev, 1)
})

test_that("a null value equal to the observed overlap counts against it", {
    ## two small unrelated sets: observed 0 overlaps and almost every shuffle gives 0 too,
    ## so p must be 1; with ties left out it was 1/(nShuffle+1)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    ## same chromosome (shuffle() keeps each peak on its chromosome; sets on different ones made
    ## every permutation warn "no sequence levels in common"), far apart
    A <- GRanges("chr1", IRanges(c(1e6, 5e7), width = 300))
    B <- GRanges("chr1", IRanges(c(2e7, 1e8), width = 300))
    set.seed(2)
    res <- enrichPeakOverlap(A, list(b = B), TxDb = txdb, nShuffle = 50, mc.cores = 1, verbose = FALSE)
    expect_equal(res$pvalue, 1)
    set.seed(2)
    res <- enrichPeakOverlap(A, list(b = B), TxDb = txdb, nShuffle = 50, pool = FALSE, mc.cores = 1,
                             verbose = FALSE)
    expect_equal(res$pvalue, 1)
})
