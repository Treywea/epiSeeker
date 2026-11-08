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

test_that("without a background a note says once per session that the whole-genome null is anti-conservative", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1001, 5001), width = 100))
    t <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1051, width = 100))
    run <- function(...) suppressWarnings(utils::capture.output(
        r <- enrichPeakOverlap(q, list(a = t), TxDb = txdb, nShuffle = 2, mc.cores = 1, ...)))
    rlang::reset_message_verbosity("epiSeeker_whole_genome_null")
    expect_message(run(verbose = TRUE), "anti-conservative")
    expect_no_message(run(verbose = TRUE), message = "anti-conservative")       # once per session
    rlang::reset_message_verbosity("epiSeeker_whole_genome_null")
    expect_no_message(run(verbose = FALSE), message = "anti-conservative")
    old <- options(epiSeeker.null_note = FALSE); on.exit(options(old), add = TRUE)
    expect_no_message(run(verbose = TRUE), message = "anti-conservative")
})

## ---- region-level test against a universe (method = "fisher") and the odds ratio ----------

## background of unequal widths on chr1, and peaks placed uniformly by base pair within it
fisher_bg <- function(n = 2000) {
    w <- round(exp(stats::runif(n, log(200), log(20000))))
    st <- cumsum(c(1e5, head(w, -1) + 5000))
    GRanges("chr1", IRanges(st, width = w))
}
in_bg <- function(n, bg, w = 300) {
    j <- sample(length(bg), n, replace = TRUE, prob = width(bg))
    room <- pmax(width(bg)[j] - w + 1L, 1L)
    GRanges("chr1", IRanges(start(bg)[j] + floor(stats::runif(n) * room), width = w))
}

test_that("the region-level test is calibrated on a universe of unequal region sizes", {
    set.seed(20261108)
    bg <- fisher_bg()
    p10 <- p1 <- numeric(200)
    for (r in seq_along(p10)) {
        q <- in_bg(400, bg); t <- in_bg(400, bg)
        p10[r] <- epiSeeker:::.overlap_fisher(q, list(t), bg, strata = 10)$pvalue
        p1[r] <- epiSeeker:::.overlap_fisher(q, list(t), bg, strata = 1)$pvalue
    }
    expect_gte(mean(p10 < 0.05), 0.01)
    expect_lte(mean(p10 < 0.05), 0.09)
    ## without strata, width confounds the table: wide units are hit by both sets
    expect_gt(mean(p1 < 0.05), mean(p10 < 0.05))
})

test_that("strata = 1 is Fisher's exact test on the region-level table", {
    set.seed(3)
    bg <- fisher_bg(300)
    q <- in_bg(80, bg); t <- c(q[1:20], in_bg(60, bg))
    ft <- epiSeeker:::.overlap_fisher(q, list(t), bg, strata = 1)
    u <- reduce(bg)
    hq <- overlapsAny(u, q); ht <- overlapsAny(u, t)
    tab <- table(factor(hq, c(TRUE, FALSE)), factor(ht, c(TRUE, FALSE)))
    expect_equal(ft$pvalue, fisher.test(tab, alternative = "greater")$p.value, tolerance = 1e-10)
    expect_equal(ft$odds_ratio, (tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1]))
    expect_equal(ft$a, tab[1, 1])
})

test_that("the region-level test works on a universe of tens of thousands of units", {
    ## the CMH variance multiplies four counts per stratum; as integers they overflowed to NA
    set.seed(11)
    bg <- GRanges("chr1", IRanges(seq(1, by = 3000, length.out = 60000), width = sample(200:2000, 60000, replace = TRUE)))
    q <- bg[sample(length(bg), 3000)]; t <- bg[sample(length(bg), 3000)]
    for (st in c(1, 2, 10)) {
        ft <- epiSeeker:::.overlap_fisher(q, list(t), bg, strata = st)
        expect_false(is.na(ft$pvalue))
        expect_false(is.na(ft$odds_ratio))
    }
})

test_that("query regions outside the universe become units, so the table stays consistent", {
    bg <- GRanges("chr1", IRanges(c(1e5, 2e5, 3e5), width = 1000))
    q <- GRanges("chr1", IRanges(c(1e5 + 10, 5e5, 6e5), width = 100))     # two outside
    t <- GRanges("chr1", IRanges(c(5e5 + 20, 9e5), width = 100))           # one outside
    ft <- epiSeeker:::.overlap_fisher(q, list(t), bg, strata = 1)
    expect_equal(ft$units, 5)
    expect_equal(ft$query_outside, 2)
    expect_equal(ft$target_outside, 1)
    expect_equal(ft$a, 1)          # the added unit at 5e5 is hit by both
    expect_false(is.na(ft$pvalue))
})

test_that("the fast background permutation counts what intersect() counts", {
    set.seed(5)
    q <- GRanges("chr1", IRanges(sample(1e6, 300), width = sample(100:2000, 300, replace = TRUE)))
    tl <- lapply(1:4, function(i) GRanges("chr1", IRanges(sample(1e6, 200), width = sample(50:3000, 200, replace = TRUE))))
    ## adjacent and overlapping regions within a set
    tl[[1]] <- c(tl[[1]], GRanges("chr1", IRanges(c(1000, 1100), c(1099, 1200))))
    bgi <- epiSeeker:::.prepare_background(GRanges("chr1", IRanges(1, 2e6)))
    obs <- epiSeeker:::.overlap_perm_background(q, tl, bgi, nShuffle = 0, verbose = FALSE)$overlap
    expect_equal(obs, vapply(tl, function(tt) length(intersect(q, tt)), numeric(1)), ignore_attr = TRUE)
})

test_that("method = \"fisher\" needs a background and reports the odds ratio and its test", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    set.seed(7)
    bg <- fisher_bg(500)
    q <- in_bg(100, bg)
    tl <- list(related = c(q[1:40], in_bg(60, bg)), unrelated = in_bg(100, bg))
    expect_error(enrichPeakOverlap(q, tl, TxDb = txdb, method = "fisher", verbose = FALSE), "needs `background`")
    res <- enrichPeakOverlap(q, tl, TxDb = txdb, method = "fisher", background = bg, verbose = FALSE)
    expect_equal(res$tSample, c("related", "unrelated"))
    expect_lt(res$pvalue[1], 1e-6)
    expect_gt(res$odds_ratio[1], res$odds_ratio[2])
    prm <- attr(res, "params")
    expect_equal(prm$method, "fisher")
    expect_match(prm$null, "Cochran-Mantel-Haenszel")
    ## the permutation path with a background reports the same odds ratio; without one, NA
    res2 <- enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 10, mc.cores = 1, background = bg, verbose = FALSE)
    expect_equal(res2$odds_ratio, res$odds_ratio)
    res3 <- enrichPeakOverlap(q, tl, TxDb = txdb, nShuffle = 0, mc.cores = 1, verbose = FALSE)
    expect_true(all(is.na(res3$odds_ratio)))
    ## the interpret() digest ranks by the recorded statistic and says so
    txt <- epiSeeker:::.digest_overlap(res, 30)$text
    expect_match(txt, paste0("ranked by ", prm$rank_by), fixed = TRUE)
    expect_match(txt, "odds_ratio=", fixed = TRUE)
    expect_match(txt, "no permutations", fixed = TRUE)
    attr(res, "params")$rank_by <- "odds_ratio"
    expect_match(epiSeeker:::.digest_overlap(res, 30)$text, "ranked by odds_ratio (jaccard", fixed = TRUE)
})
