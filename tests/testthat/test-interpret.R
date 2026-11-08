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
    expect_equal(names(seen$schema$properties), c("answer", "candidates", "binding_profile",
                 "assay_consistency", "target_programs", "cofactors", "narrative"))
    expect_match(seen$prompt, "Genome: hg38 (TxDb: UCSC", fixed = TRUE)
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
    expect_equal(res$candidates, list())
    expect_true(is.na(res$assay_consistency$verdict))
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
    for (cl in c("csAnno", "GRanges", "csAnnoList", "peakOverlapResult"))
        expect_true(methods::existsMethod(clusterProfiler::interpret, cl))
})

## a mock that answers with the empty object of whatever schema it is given
mock_answer <- function(env) {
    local_mocked_bindings(
        generate_object = function(...) {
            a <- list(...); env$seen <- a
            list(object = lapply(a$schema$properties, epiSeeker:::.na_of), finish_reason = "stop", raw_text = "{}")
        },
        .package = "aisdk", .env = parent.frame()
    )
}

test_that("a GRanges is annotated first, and only against a TxDb that was given", {
    skip_if_not_installed("aisdk")
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    expect_error(clusterProfiler::interpret(peak), "Give `TxDb`")
    env <- new.env(); mock_answer(env)
    expect_message(res <- clusterProfiler::interpret(peak, TxDb = txdb, annoDb = "org.Hs.eg.db", assay = "TF",
                                                     model = "openai:m"), "Annotating 150 regions")
    expect_equal(provenance <- clusterProfiler::provenance(res)$entry_point, "interpret,csAnno")
    ## the same prompt as annotating by hand
    by_hand <- epiSeeker:::.seq_digest(annotateSeq(peak, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE),
                                       assay = "TF")
    expect_identical(env$seen$prompt, by_hand$text)
})

test_that("a csAnnoList is compared group by group", {
    skip_if_not_installed("aisdk")
    anno <- interpret_anno()
    half <- anno
    half@anno <- anno@anno[1:60]
    half@annoStat <- epiSeeker:::getGenomicAnnoStat(half@anno)
    x <- csAnnoList(all = anno, top = half)
    ev <- epiSeeker:::.seq_digest_list(x, assay = "TF", n_regions = 5)
    expect_equal(names(ev$sections), c("overview", "distribution", "tss_distance", "width", "genes", "top_regions"))
    expect_match(ev$sections$overview, "all (150 regions), top (60 regions)", fixed = TRUE)
    expect_match(ev$sections$top_regions, "Top regions of all (5 of 150", fixed = TRUE)
    expect_match(ev$sections$top_regions, "Top regions of top (5 of 60", fixed = TRUE)
    ## half is a subset of all, so nothing is specific to it
    expect_length(ev$evidence$specific_promoter_genes$top, 0)
    env <- new.env(); mock_answer(env)
    res <- clusterProfiler::interpret(x, model = "openai:m")
    expect_equal(names(env$seen$schema$properties), c("answer", "candidates", "group_profiles", "differences",
                                                      "shared", "target_programs", "narrative"))
    expect_equal(clusterProfiler::provenance(res)$entry_point, "interpret,csAnnoList")
    expect_output(print(res), "epiSeeker interpretation")
})

test_that("an overlap table alone points to the csAnno method, and its test is described", {
    ov <- toy_overlap()
    class(ov) <- c("peakOverlapResult", "data.frame")
    expect_error(clusterProfiler::interpret(ov), "overlap = ")
    attr(ov, "params") <- list(nShuffle = 1000, pool = TRUE)
    txt <- epiSeeker:::.digest_overlap(ov, 30)$text
    expect_match(txt, "p-values from 1000 shuffles; null: pooled over reference sets, shuffled over whole chromosomes",
                 fixed = TRUE)
    expect_match(txt, "smallest possible p is 0.000999", fixed = TRUE)
    expect_match(txt, "rely on the ranking, not on the p-values", fixed = TRUE)
    ## only columns the table has are named as effect sizes
    expect_match(txt, "ranked by jaccard)", fixed = TRUE)
    ## the null recorded by enrichPeakOverlap() is stated, and a background null gets no caution
    attr(ov, "params") <- list(nShuffle = 100, pool = TRUE, null = "per target set, shuffled within background")
    txt <- epiSeeker:::.digest_overlap(ov, 30)$text
    expect_match(txt, "null: per target set, shuffled within background", fixed = TRUE)
    expect_false(grepl("anti-conservative", txt))
    ## ranking stays on jaccard when fold_enrichment is present
    ov$fold_enrichment <- c(1, 50, 2, 3)
    expect_equal(epiSeeker:::.digest_overlap(ov, 30)$table$label, c("CTCF", "RAD21", "CREB1"))
    expect_match(epiSeeker:::.digest_overlap(ov, 30)$text,
                 "ranked by jaccard (fold_enrichment is reported as an effect size", fixed = TRUE)
})

test_that("the digest reports what the checks of annotateSeq found", {
    anno <- interpret_anno()
    ## 150 peaks on the right build: consistent, but the check is unreliable for so few peaks
    expect_match(epiSeeker:::.seq_digest(anno)$sections$overview,
                 "Checks: genome build not contradicted, but the check is unreliable with 150 peaks", fixed = TRUE)
    ## a suspected mismatch and dropped peaks are stated in the overview
    S4Vectors::metadata(anno@anno)$annoChecks$build$suspect <- TRUE
    S4Vectors::metadata(anno@anno)$annoChecks$build$reasons <- "3 peak(s) extend beyond the chromosome ends of TxDb (chr19)"
    S4Vectors::metadata(anno@anno)$annoChecks$dropped <- list(n = 2, n_total = 152,
        by_seqlevel = data.frame(seqlevel = "chrEBV", n = 2, reason = "not in TxDb"))
    ov <- epiSeeker:::.seq_digest(anno)$sections$overview
    expect_match(ov, "GENOME BUILD SUSPECT - 3 peak(s) extend beyond", fixed = TRUE)
    expect_match(ov, "2 of 152 peaks dropped (not in TxDb: chrEBV 2)", fixed = TRUE)
    ## objects made before the checks existed get no line
    data(peakAnno, package = "epiSeeker")
    expect_false(grepl("Checks:", epiSeeker:::.seq_digest(peakAnno)$sections$overview))
    ## and the model is told what to do with it
    expect_match(epiSeeker:::.seq_system_prompt(), "GENOME BUILD SUSPECT", fixed = TRUE)
})
