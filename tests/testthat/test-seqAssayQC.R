## seqAssayQC(): the declared assay against ENCODE ranges (R/seqAssayQC.R)

rng <- epiSeeker:::.assayqc_ranges()

## metrics in the middle of the training range of a class
mid_metrics <- function(cl, n = 5000) {
    r <- rng[rng$class == cl, ]
    v <- setNames((r$train_min + r$train_max) / 2, r$feature)
    c(n_peaks = n, v)
}

## 600 domains of 20 kb on chr1: nothing like a transcription factor
domain_anno <- function(assay = NULL) {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1e6, by = 2e5, length.out = 600), width = 20000))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    annotateSeq(gr, TxDb = txdb, verbose = FALSE, assay = assay)
}

test_that("the stored range table is small and complete", {
    f <- system.file("extdata", "assayQC_ranges.tsv", package = "epiSeeker")
    expect_lt(file.size(f), 10000)
    expect_true(all(c("TF", "H3K4me3", "H3K27me3", "H3K36me3", "ATAC", "DNase") %in% rng$class))
    expect_true(all(table(rng$class) == 3))
    expect_true(all(rng$lower <= rng$train_min & rng$upper >= rng$train_max))
    expect_true(all(rng$n_train >= 5))
})

test_that("declared assays are matched to the table by name", {
    cls <- unique(rng$class)
    cl <- function(s) epiSeeker:::.assay_class(s, cls)$class
    expect_equal(cl("H3K4me3"), "H3K4me3")
    expect_equal(cl("h3k27ac ChIP-seq"), "H3K27ac")
    expect_equal(cl("CUT&Tag H3K27me3"), "H3K27me3")
    expect_equal(cl("ATAC-seq"), "ATAC")
    expect_equal(cl("DNase-seq"), "DNase")
    expect_equal(cl("TF"), "TF")
    expect_equal(cl("transcription factor ChIP-seq"), "TF")
    expect_equal(cl("CTCF"), "TF")
    expect_match(epiSeeker:::.assay_class("CTCF", cls)$note, "taken as a ChIP of a transcription factor")
    expect_equal(cl("TFAP2A ChIP-seq"), "TF")
    expect_true(is.na(cl("H3K4me2")))
    expect_match(epiSeeker:::.assay_class("H3K4me2", cls)$note, "no ENCODE reference range for H3K4me2")
    expect_true(is.na(cl("Input")))
    expect_true(is.na(cl("H2A.Z")))
    expect_true(is.na(cl("a ChIP of something in liver")))
})

test_that("a set in the middle of its own range is consistent; the PLAN example is flagged", {
    d <- epiSeeker:::.assayqc_decide
    for (cl in unique(rng$class)) {
        q <- d(mid_metrics(cl), cl)
        expect_true(q$assessable, info = cl)
        expect_true(q$consistent, info = cl)
        expect_length(q$reasons, 0)
        expect_true(cl %in% q$fits, info = cl)
    }
    ## declared H3K4me3, but only 10% of the peaks at promoters
    m <- mid_metrics("H3K4me3"); m[["promoter_1kb"]] <- 0.10
    q <- d(m, "H3K4me3")
    expect_false(q$consistent)
    expect_match(q$reasons, "promoter fraction \\(<= 1 kb from a TSS\\) 10% is below .*ENCODE H3K4me3 peak sets")
    expect_equal(q$ranges$feature, c("promoter_1kb", "width_median", "gene_body"))
    expect_equal(q$ranges$inside, c(FALSE, TRUE, TRUE))
    ## too few peaks, or no range for the assay: not assessable, never flagged
    q <- d(mid_metrics("H3K4me3", n = 150), "H3K4me3")
    expect_false(q$assessable); expect_true(is.na(q$consistent))
    expect_match(q$reasons, "only 150 peaks")
    q <- d(mid_metrics("H3K4me3"), "H3K4me2")
    expect_false(q$assessable); expect_true(is.na(q$consistent))
    expect_match(epiSeeker:::.assayqc_line(q), "^assay QC: not assessable")
})

test_that("annotateSeq(assay = ) records the check, warns, and show() says so", {
    expect_warning(x <- domain_anno("TF"), "do not look like TF: median width 20,000 bp is above")
    qc <- annoChecks(x)$assay
    expect_identical(qc, seqAssayQC(x, "TF"))
    expect_false(qc$consistent)
    expect_equal(qc$class, "TF")
    expect_equal(qc$metrics[["n_peaks"]], 600)
    expect_output(show(x), "do not look like the declared assay, TF")
    expect_match(qc$source[1], "ENCODE GRCh38")
    ## without assay: nothing recorded, no warning
    expect_no_warning(y <- domain_anno())
    expect_null(annoChecks(y)$assay)
    expect_error(seqAssayQC(y, NULL), "one declared assay")
    expect_error(seqAssayQC(y@anno, "TF"), "csAnno")
})

test_that("interpret() puts the flag in the digest and the model cannot call the set consistent", {
    skip_if_not_installed("aisdk")
    x <- suppressWarnings(domain_anno())
    ev <- epiSeeker:::.seq_digest(x, assay = "TF")
    expect_match(ev$sections$overview, "Declared assay: TF\nASSAY QC INCONSISTENT - the regions do not look like the declared assay (TF): median width",
                 fixed = TRUE)
    expect_false(ev$assay_qc$consistent)
    expect_match(epiSeeker:::.seq_system_prompt(), "ASSAY QC INCONSISTENT", fixed = TRUE)
    ## no assay, no line
    expect_false(grepl("assay QC", epiSeeker:::.seq_digest(x)$sections$overview, ignore.case = TRUE))
    seen <- NULL
    local_mocked_bindings(
        generate_object = function(...) {
            seen <<- list(...)
            list(object = list(answer = "None", candidates = list(), binding_profile = "b",
                               assay_consistency = list(verdict = "Consistent", reason = "looks fine"),
                               target_programs = "t", cofactors = "None", narrative = "n"),
                 finish_reason = "stop", raw_text = "{}")
        },
        .package = "aisdk"
    )
    res <- clusterProfiler::interpret(x, assay = "TF", model = "openai:m")
    expect_match(seen$prompt, "ASSAY QC INCONSISTENT", fixed = TRUE)
    expect_equal(res$assay_consistency$verdict, "Inconsistent")
    expect_match(res$assay_consistency$reason, "^\\[set by seqAssayQC\\(\\), the model said Consistent\\] median width")
    p <- clusterProfiler::provenance(res)
    expect_false(p$assay_qc$consistent)
    expect_true(p$assay_qc$verdict_overridden)
    ## a set the check does not flag keeps the model's verdict
    res2 <- clusterProfiler::interpret(x, assay = "H3K4me2", model = "openai:m")
    expect_equal(res2$assay_consistency$verdict, "Consistent")
    expect_false(clusterProfiler::provenance(res2)$assay_qc$verdict_overridden)
    ## a csAnnoList gets one line per group
    evl <- epiSeeker:::.seq_digest_list(csAnnoList(a = x, b = x), assay = "TF", n_regions = 3)
    expect_match(evl$sections$overview, "(a) ASSAY QC INCONSISTENT", fixed = TRUE)
    expect_match(evl$sections$overview, "(b) ASSAY QC INCONSISTENT", fixed = TRUE)
})
