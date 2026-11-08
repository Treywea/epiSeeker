library(epiSeeker)

context("test that every epiSeeker result reaches its own interpret() method")

## answers with the empty object of whatever schema it is given, and remembers the request
mock_any <- function(env) {
    local_mocked_bindings(
        generate_object = function(...) {
            a <- list(...); env$calls <- c(env$calls, list(a))
            list(object = lapply(a$schema$properties, epiSeeker:::.na_of), finish_reason = "stop",
                 raw_text = "{}")
        },
        .package = "aisdk", .env = parent.frame()
    )
}

test_that("each class with an interpret() method dispatches to its own method", {
    for (cl in c("csAnno", "GRanges", "csAnnoList", "peakOverlapResult", "bmData", "loopAnno")) {
        m <- methods::selectMethod(clusterProfiler::interpret, cl)
        expect_identical(m@defined[[1]], cl, info = cl)
    }
    ## enrichment results still go to clusterProfiler's own method
    expect_identical(methods::selectMethod(clusterProfiler::interpret, "data.frame")@defined[[1]], "ANY")
})

test_that("objects made by the adapters and the bulk functions are all interpreted", {
    skip_if_not_installed("aisdk")
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    env <- new.env(); env$calls <- list(); mock_any(env)
    entry <- function(r) clusterProfiler::provenance(r)$entry_point

    ## range: csAnno, and a GRanges annotated on the way
    data(peakAnno, package = "epiSeeker")
    expect_equal(entry(clusterProfiler::interpret(peakAnno, model = "openai:m")), "interpret,csAnno")
    r_gr <- suppressMessages(clusterProfiler::interpret(peak, TxDb = txdb, model = "openai:m"))
    expect_equal(entry(r_gr), "interpret,csAnno")

    ## range, grouped: a csAnnoList from annotateGroups() states its resolution
    g <- GenomicRanges::GRangesList(c1 = peak[1:60], c2 = peak[61:150])
    grp <- annotateGroups(g, TxDb = txdb, resolution = "single-cell", members = c(c1 = 120, c2 = 2000))
    r_grp <- clusterProfiler::interpret(grp, model = "openai:m")
    expect_equal(entry(r_grp), "interpret,csAnnoList")
    expect_match(attr(r_grp, "evidence")$sections$overview, "Resolution: single-cell", fixed = TRUE)

    ## base: a bmData summarised from single-molecule calls
    f <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste("#ct", "st", "en", "fiber", "score", "strand", "ref_5mC", "5mC_qual", sep = "\t"),
        paste("chr1", 1000, 2000, "r1", 0, "+", "1500,1600", "240,20", sep = "\t"),
        paste("chr1", 1000, 2000, "r2", 0, "-", "1500,1600", "10,250", sep = "\t")
    ), f)
    bm <- summarizeModCalls(readModCalls(f), type = "5mC", sampleName = "s")
    r_bm <- clusterProfiler::interpret(bm, model = "openai:m")
    expect_equal(entry(r_bm), "interpret,bmData")

    ## pair: a loopAnno restricted to some genes says so in the digest
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10000, 100000, 300000), c(20000, 110000, 310000)),
                                    strand = c("+", "+", "-"), gene_id = c("G1", "G2", "G3"))
    bedpe <- data.frame(c1 = "chr1", s1 = c(40000, 8000), e1 = c(45000, 12000),
                        c2 = "chr1", s2 = c(98000, 308000), e2 = c(103000, 312000))
    la <- suppressWarnings(suppressMessages(annotateLoops(bedpe, TxDb = genes, gene_ids = c("G1", "G2"),
                                                          verbose = FALSE)))
    r_la <- clusterProfiler::interpret(la, model = "openai:m")
    expect_equal(entry(r_la), "interpret,loopAnno")
    expect_match(attr(r_la, "evidence")$sections$overview, "of one of 2 genes given (gene_ids)", fixed = TRUE)

    ## one model call each, and each prompt is the digest the result carries
    expect_length(env$calls, 5)
    expect_identical(env$calls[[5]]$prompt, attr(r_la, "evidence")$text)

    ## a lone overlap table points to the csAnno method
    ov <- data.frame(tSample = "a", N_OL = 1)
    class(ov) <- c("peakOverlapResult", "data.frame")
    expect_error(clusterProfiler::interpret(ov), "overlap = ")
})

test_that("an overlap table from method = \"fisher\" reaches the csAnno digest with its odds ratio and ranking", {
    skip_if_not_installed("aisdk")
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    env <- new.env(); env$calls <- list(); mock_any(env)
    set.seed(11)
    bg <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1e5, by = 2e4, length.out = 300), width = 5000))
    pick <- function(n) {
        j <- sample(length(bg), n, replace = TRUE)
        GenomicRanges::GRanges("chr1", IRanges::IRanges(GenomicRanges::start(bg)[j] + sample(0:4000, n, TRUE),
                                                        width = 300))
    }
    q <- pick(100)
    tl <- list(related = c(q[1:40], pick(60)), unrelated = pick(100))
    res <- enrichPeakOverlap(q, tl, TxDb = txdb, method = "fisher", background = bg, verbose = FALSE)
    expect_s3_class(res, "peakOverlapResult")
    expect_true("odds_ratio" %in% names(res))
    rank_by <- attr(res, "params")$rank_by
    expect_true(rank_by %in% names(res))

    data(peakAnno, package = "epiSeeker")
    r <- clusterProfiler::interpret(peakAnno, overlap = res, model = "openai:m")
    expect_equal(clusterProfiler::provenance(r)$entry_point, "interpret,csAnno")
    sec <- attr(r, "evidence")$sections$overlap
    expect_match(sec, paste0("ranked by ", rank_by), fixed = TRUE)
    expect_match(sec, "odds_ratio=", fixed = TRUE)
    expect_match(sec, "no permutations", fixed = TRUE)
    expect_match(sec, "related", fixed = TRUE)
    ## the model is asked once, with the digest that names the overlap
    expect_length(env$calls, 1)
    expect_match(env$calls[[1]]$prompt, "odds_ratio=", fixed = TRUE)
})
