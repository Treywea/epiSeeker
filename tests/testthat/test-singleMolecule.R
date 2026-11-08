library(epiSeeker)

context("test the single-molecule adapter")

## two Fiber-seq reads on chr1: r1 opens 1101-1160 (FIRE) with m6A in it, r2 is closed there
ft_table <- function() {
    f <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste("#ct", "st", "en", "fiber", "score", "strand", "HP", "ref_m6a", "m6a_qual", "ref_5mC", "5mC_qual",
              "ref_nuc_starts", "ref_nuc_lengths", "ref_msp_starts", "ref_msp_lengths", "fire", sep = "\t"),
        paste("chr1", 1000, 2000, "r1", 0, "+", 1, "1110,1120,1130,-1", "250,250,100,250", "1500", "240",
              "1300", "147", "1100", "60", "240", sep = "\t"),
        paste("chr1", 1000, 2000, "r2", 0, "-", 2, "1700", "250", "1500,1600", "10,250",
              "1080", "147", "1700", "50", "100", sep = "\t")
    ), f)
    f
}

test_that("fibertools tables become 1-based per-call ranges, with -1 dropped and counted", {
    x <- readModCalls(ft_table())
    expect_s4_class(x, "GRanges")
    expect_equal(sort(unique(x$type)), c("5mC", "m6A", "msp", "nuc"))
    m6a <- x[x$type == "m6A"]
    expect_equal(GenomicRanges::start(m6a[m6a$read_id == "r1"]), c(1111, 1121, 1131))
    expect_equal(m6a$modified[m6a$read_id == "r1"], c(TRUE, TRUE, FALSE))     # 100/255 < 0.5
    expect_equal(S4Vectors::metadata(x)$modSource$n_unplaced, 1)
    msp <- x[x$type == "msp"]
    expect_equal(GenomicRanges::start(msp), c(1101, 1701))
    expect_equal(GenomicRanges::end(msp), c(1160, 1750))
    expect_equal(msp$modified, c(TRUE, FALSE))                               # FIRE score > 230
    reads <- S4Vectors::metadata(x)$reads
    expect_equal(reads$start, c(1001, 1001))
    expect_equal(reads$HP, c(1L, 2L))
    ## region filter keeps reads overlapping it and only their calls inside it
    y <- readModCalls(ft_table(), region = GenomicRanges::GRanges("chr1:1100-1200"))
    expect_true(all(GenomicRanges::start(y) <= 1200 & GenomicRanges::end(y) >= 1100))
})

test_that("modkit extract full and extract calls are read", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste("read_id", "forward_read_position", "ref_position", "chrom", "mod_strand", "ref_strand",
              "ref_mod_strand", "fw_soft_clipped_start", "fw_soft_clipped_end", "alignment_start",
              "alignment_end", "read_length", "mod_qual", "mod_code", "base_qual", "ref_kmer", "query_kmer",
              "canonical_base", "modified_primary_base", "inferred", "flag", sep = "\t"),
        paste("a", 5, 99, "chr2", "+", "+", "+", 0, 0, 90, 300, 210, 0.9, "m", 30, "ACGTA", "ACGTA", "C", "C", "false", 0, sep = "\t"),
        paste("a", 9, 109, "chr2", "+", "+", "+", 0, 0, 90, 300, 210, 0.1, "m", 30, "ACGTA", "ACGTA", "C", "C", "false", 0, sep = "\t"),
        paste("b", 3, -1, "chr2", "+", "+", "+", 0, 0, 90, 300, 210, 0.9, "m", 30, "ACGTA", "ACGTA", "C", "C", "false", 0, sep = "\t")
    ), f)
    x <- readModCalls(f)
    expect_equal(S4Vectors::metadata(x)$modSource$format, "modkit")
    expect_equal(GenomicRanges::start(x), c(100, 110))
    expect_equal(x$type, c("5mC", "5mC"))
    expect_equal(x$modified, c(TRUE, FALSE))
    expect_equal(S4Vectors::metadata(x)$modSource$n_unplaced, 1)
    expect_equal(S4Vectors::metadata(x)$reads$start, 91)

    g <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste("read_id", "ref_position", "chrom", "ref_strand", "call_prob", "call_code", "canonical_base", "fail", sep = "\t"),
        paste("a", 99, "chr2", "+", 0.95, "m", "C", "false", sep = "\t"),
        paste("a", 109, "chr2", "+", 0.90, "-", "C", "false", sep = "\t"),
        paste("a", 119, "chr2", "+", 0.55, "m", "C", "true", sep = "\t")
    ), g)
    y <- readModCalls(g, format = "modkit")
    expect_equal(length(y), 2)                                               # the failed call is dropped
    expect_equal(y$type, c("5mC", "5mC"))
    expect_equal(y$modified, c(TRUE, FALSE))
})

test_that("calls summarise to a bmData that the existing functions take", {
    x <- readModCalls(ft_table())
    b <- summarizeModCalls(x, type = "5mC", sampleName = "s")
    expect_s4_class(b, "bmData")
    expect_equal(SummarizedExperiment::assayNames(b), c("Cov", "Methylation"))
    ## 5mC lists unmodified calls, so coverage = calls: 1501 has 2 calls, 1 modified (240 >= 128)
    v <- SummarizedExperiment::assay(b, "Methylation")[GenomicRanges::start(b) == 1501, 1]
    expect_equal(unname(v), 0.5)
    ## m6A of Fiber-seq lists only modified bases -> coverage = reads spanning the position
    expect_equal(epiSeeker:::.sm_lists_negatives(x[x$type == "m6A" & x$modified]), FALSE)
    a <- summarizeModCalls(x[x$type != "m6A" | x$modified], type = "m6A")
    expect_true(all(SummarizedExperiment::assay(a, "Cov") == 2))
    f <- tempfile(fileext = ".rds"); saveRDS(x, f); expect_identical(readRDS(f), x)
})

test_that("bedMethyl reads as bmData", {
    f <- tempfile(fileext = ".bed")
    writeLines(c(
        "chr1\t99\t100\tm\t10\t+\t99\t100\t255,0,0\t10\t80.00\t8\t2\t0\t0\t0\t0\t0",
        "chr1\t149\t150\tm\t5\t+\t149\t150\t255,0,0\t5\t20.00\t1\t4\t0\t0\t0\t0\t0",
        "chr1\t149\t150\th\t5\t+\t149\t150\t255,0,0\t5\t0.00\t0\t5\t0\t0\t0\t0\t0"
    ), f)
    expect_error(readBedMethyl(f), "several modification codes")
    b <- readBedMethyl(f, mod = "m", sampleName = "s1")
    expect_s4_class(b, "bmData")
    expect_equal(GenomicRanges::start(b), c(100, 150))
    expect_equal(unname(SummarizedExperiment::assay(b, "Methylation")[, 1]), c(0.8, 0.2))
    expect_equal(unname(SummarizedExperiment::assay(b, "Cov")[, 1]), c(10, 5))
})

test_that("per-read matrix and co-occurrence go in the right direction", {
    ## 20 reads: in 10 both regions are open, in 10 both closed -> positive association
    rows <- vapply(seq_len(20), function(i) {
        open <- i <= 10
        paste("chr1", 0, 3000, paste0("r", i), 0, "+", ".",
              if (open) "100,1100" else ".", if (open) "250,250" else ".",
              if (open) "100,1100" else "2500", if (open) "50,50" else "50", if (open) "240,240" else "100", sep = "\t")
    }, "")
    f <- tempfile(fileext = ".tsv")
    writeLines(c(paste("#ct", "st", "en", "fiber", "score", "strand", "HP", "ref_m6a", "m6a_qual",
                       "ref_msp_starts", "ref_msp_lengths", "fire", sep = "\t"), rows), f)
    x <- readModCalls(f)
    r1 <- GenomicRanges::GRanges("chr1:101-150"); r2 <- GenomicRanges::GRanges("chr1:1101-1150")
    rr <- c(r1, r2); names(rr) <- c("a", "b")
    m <- modReadMatrix(x, rr, type = "msp")
    expect_equal(dim(m), c(20, 2))
    expect_equal(unname(m[1:10, 1]), rep(1, 10))
    expect_equal(unname(m[11:20, 1]), rep(0, 10))
    co <- modCooccurrence(x, r1, r2, type = "msp")
    expect_equal(co$n_reads, 20)
    expect_equal(as.vector(co$table), c(10, 0, 0, 10))
    expect_lt(co$p_value, 0.001)
    ## independence: region 2 open in half of each group
    rows2 <- vapply(seq_len(20), function(i) {
        o1 <- i <= 10; o2 <- i %% 2 == 0
        st <- c(if (o1) 100, if (o2) 1100)
        paste("chr1", 0, 3000, paste0("q", i), 0, "+", ".", ".", ".",
              if (length(st)) paste(st, collapse = ",") else ".",
              if (length(st)) paste(rep(50, length(st)), collapse = ",") else ".",
              if (length(st)) paste(rep(240, length(st)), collapse = ",") else ".", sep = "\t")
    }, "")
    g <- tempfile(fileext = ".tsv")
    writeLines(c(paste("#ct", "st", "en", "fiber", "score", "strand", "HP", "ref_m6a", "m6a_qual",
                       "ref_msp_starts", "ref_msp_lengths", "fire", sep = "\t"), rows2), g)
    co2 <- modCooccurrence(readModCalls(g), r1, r2, type = "msp")
    expect_gt(co2$p_value, 0.5)
    expect_s3_class(plotModReads(x, GenomicRanges::GRanges("chr1:1-3000"), peaks = c(r1, r2)), "ggplot")
})
