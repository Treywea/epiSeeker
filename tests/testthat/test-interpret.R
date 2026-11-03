library(epiSeeker)

context("test interpret() for csAnno: evidence digest and model call")

skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
skip_if_not_installed("org.Hs.eg.db")

## the 150 peaks of the bundled sample, which carry fold_enrichment
interpret_anno <- function() {
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    annotateSeq(peak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                annoDb = "org.Hs.eg.db", verbose = FALSE)
}

toy_overlap <- function() {
    data.frame(qSample = "q", tSample = c("a.bed", "b.bed", "c.bed", "d.bed"),
               N_OL = c(50, 90, 10, 70), p.adjust = c(0.001, 0.001, 0.2, 0.01),
               jaccard = c(0.30, 0.05, 0.40, 0.10), label = c("CTCF", "CREB1", "noise", "RAD21"))
}

test_that("the digest sends specific evidence first and says how much it sends", {
    anno <- interpret_anno()
    ev <- epiSeeker:::.seq_digest(anno, context = "ctx", assay = "TF", question = "Which factor?",
                                  overlap = toy_overlap(), n_regions = 20)
    expect_equal(names(ev$sections), c("overview", "top_regions", "genes", "overlap",
                                       "distribution", "tss_distance", "width"))
    expect_match(ev$sections$top_regions, "20 of 150, ranked by fold_enrichment", fixed = TRUE)
    expect_equal(nrow(ev$evidence$regions), 20)
    ## ranked by signal, highest first
    df <- as.data.frame(anno@anno)
    expect_equal(ev$evidence$regions$start[1], df$start[which.max(df$fold_enrichment)])
    expect_match(ev$sections$overview, "Declared assay: TF", fixed = TRUE)
    expect_match(ev$sections$genes, "Promoter genes \\([0-9]+ of [0-9]+")
    ## the same object gives the same digest
    expect_identical(ev$text, epiSeeker:::.seq_digest(anno, context = "ctx", assay = "TF",
        question = "Which factor?", overlap = toy_overlap(), n_regions = 20)$text)
    skip_if_not_installed("digest")
    expect_true(all(nchar(ev$sha256) == 64))
})

test_that("overlap keeps significant samples and ranks them by effect size, not p-value", {
    ov <- epiSeeker:::.digest_overlap(toy_overlap(), n_overlap = 30)
    expect_equal(ov$table$label, c("CTCF", "RAD21", "CREB1"))   # c.bed has p.adjust 0.2
    expect_match(ov$text, "3 of 3 with p.adjust < 0.05, 4 tested; ranked by jaccard", fixed = TRUE)
    ## counts whole, also for a single row
    big <- transform(toy_overlap()[1, ], N_OL = 36797)
    expect_match(epiSeeker:::.digest_overlap(big, 30)$text, "N_OL=36797, jaccard=0.3, p.adjust=0.001", fixed = TRUE)
    expect_error(epiSeeker:::.digest_overlap(data.frame(x = 1), 30), "tSample")
})

test_that("without a signal column regions are ranked by distance to the TSS", {
    anno <- interpret_anno()
    GenomicRanges::mcols(anno@anno)$fold_enrichment <- NULL
    GenomicRanges::mcols(anno@anno)$tags <- NULL
    ev <- epiSeeker:::.seq_digest(anno, n_regions = 5)
    expect_match(ev$totals$rank_by, "distanceToTSS")
    df <- as.data.frame(anno@anno)
    expect_equal(ev$evidence$regions$start, df$start[order(abs(df$distanceToTSS))][1:5])
    expect_error(epiSeeker:::.seq_digest(anno, signal = "nope"), "not found")
})

test_that("interpret() dispatches to the csAnno method and sends the digest as the prompt", {
    skip_if_not_installed("aisdk")
    anno <- interpret_anno()
    seen <- NULL
    local_mocked_bindings(
        generate_object = function(...) {
            seen <<- list(...)
            list(object = list(answer = "CTCF", candidates = list(list(label = "CTCF", support = "a.bed")),
                               binding_profile = "b", assay_consistency = list(verdict = "Consistent", reason = "r"),
                               target_programs = "t", cofactors = "c", narrative = "n"),
                 finish_reason = "stop", raw_text = "{}")
        },
        .package = "aisdk"
    )
    res <- clusterProfiler::interpret(anno, context = "ctx", assay = "TF", question = "Which factor?",
                                      overlap = toy_overlap(), model = "openai:gpt-4o-mini")
    expect_s3_class(res, "seq_interpretation")
    expect_equal(res$answer, "CTCF")
    expect_identical(seen$prompt, attr(res, "evidence")$text)
    expect_false("temperature" %in% names(seen))          # NULL is not sent
    expect_equal(seen$max_tokens, 4096)
    expect_equal(names(seen$schema$properties), epiSeeker:::.seq_fields)
    p <- clusterProfiler::provenance(res)
    expect_true(p$ok)
    expect_equal(p$entry_point, "interpret,csAnno")
    expect_equal(p$calls[[1]]$model_requested, "gpt-4o-mini")
    expect_output(print(p), "interpret provenance")
    expect_output(print(res), "CTCF")
})

test_that("an unfinished answer is reported as a failure, not as a result", {
    skip_if_not_installed("aisdk")
    anno <- interpret_anno()
    local_mocked_bindings(
        generate_object = function(...) list(object = list(answer = "half"), finish_reason = "length", raw_text = "{"),
        .package = "aisdk"
    )
    expect_warning(res <- clusterProfiler::interpret(anno, model = "openai:gpt-4o-mini"), "no complete answer")
    expect_true(is.na(res$answer))
    p <- clusterProfiler::provenance(res)
    expect_false(p$ok)
    expect_equal(p$calls[[1]]$status, "incomplete")
    expect_output(print(res), "INCOMPLETE")

    local_mocked_bindings(generate_object = function(...) stop("HTTP 401"), .package = "aisdk")
    expect_warning(res <- clusterProfiler::interpret(anno, model = "openai:gpt-4o-mini"), "no complete answer")
    expect_equal(clusterProfiler::provenance(res)$calls[[1]]$status, "call_failed")
    expect_equal(clusterProfiler::provenance(res)$calls[[1]]$error_message, "HTTP 401")
})

test_that("enrichment results still go to clusterProfiler's own method", {
    m <- methods::selectMethod(clusterProfiler::interpret, "data.frame")
    expect_identical(m@defined[[1]], "ANY")
    expect_true(methods::existsMethod(clusterProfiler::interpret, "csAnno"))
})
