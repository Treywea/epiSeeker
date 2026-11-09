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
    ## fibertools lists only calls above its ML filter -> coverage = reads spanning, of the strands
    ## that gave 5mC calls: 1501 has calls from r1 (+, modified) and r2 (-), 1601 only from r2 (-)
    v <- SummarizedExperiment::assay(b, "Methylation")[GenomicRanges::start(b) == 1501, 1]
    expect_equal(unname(v), 0.5)
    cv <- SummarizedExperiment::assay(b, "Cov")[, 1]
    expect_equal(unname(cv[GenomicRanges::start(b) == 1501]), 2)
    expect_equal(unname(cv[GenomicRanges::start(b) == 1601]), 1)
    ## the same calls from a tool that lists unmodified calls -> coverage = calls
    xm <- x; S4Vectors::metadata(xm)$modSource$format <- "modkit"
    expect_true(epiSeeker:::.sm_lists_negatives(xm[xm$type == "5mC"]))
    expect_false(epiSeeker:::.sm_lists_negatives(x[x$type == "5mC"]))
    bm <- summarizeModCalls(xm, type = "5mC")
    expect_equal(unname(SummarizedExperiment::assay(bm, "Cov")[GenomicRanges::start(bm) == 1501, 1]), 2)
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

test_that("modReadMatrix values by type, NA for reads not spanning the region", {
    x <- readModCalls(ft_table())
    rr <- GenomicRanges::GRanges(c("chr1:1101-1160", "chr1:1701-1800", "chr1:900-1100"))
    m <- modReadMatrix(x, rr, type = "msp")
    expect_equal(unname(m["r1", ]), c(1, 0, NA))
    expect_equal(unname(m["r2", ]), c(0, 0.5, NA))
    expect_equal(unname(modReadMatrix(x, rr[1], type = "nuc")[, 1]), c(0, 60 / 60))   # r2 nucleosome 1081-1227
    ## fibertools (only calls above its ML filter): modified calls per kb
    m6 <- GenomicRanges::GRanges("chr1:1101-1200")
    expect_equal(unname(modReadMatrix(x, m6, type = "m6A")[, 1]), c(20, 0))
    ## a table that lists unmodified calls: fraction modified, NA where a read has no call
    xm <- x; S4Vectors::metadata(xm)$modSource$format <- "modkit"
    expect_equal(unname(modReadMatrix(xm, m6, type = "m6A")[, 1]), c(2 / 3, NA))
    m5 <- modReadMatrix(xm, GenomicRanges::GRanges(c("chr1:1450-1650", "chr1:1700-1800")), type = "5mC")
    expect_equal(unname(m5[, 1]), c(1, 0.5))
    expect_true(all(is.na(m5[, 2])))
    ## plotModReads keeps the reads covering most of the region
    expect_s3_class(plotModReads(x, GenomicRanges::GRanges("chr1:1001-2000"), max_reads = 1), "ggplot")
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

## fibertools-like table: read i spans 0-6000, MSPs given as 0-based starts and lengths
sm_reads_table <- function(msps) {
    rows <- vapply(seq_along(msps), function(i) {
        s <- msps[[i]]
        paste("chr1", 0, 6000, paste0("r", i), 0, "+", ".", ".", ".",
              if (length(s)) paste(s[, 1], collapse = ",") else ".",
              if (length(s)) paste(s[, 2], collapse = ",") else ".",
              if (length(s)) paste(rep(100, nrow(s)), collapse = ",") else ".", sep = "\t")
    }, "")
    f <- tempfile(fileext = ".tsv")
    writeLines(c(paste("#ct", "st", "en", "fiber", "score", "strand", "HP", "ref_m6a", "m6a_qual",
                       "ref_msp_starts", "ref_msp_lengths", "fire", sep = "\t"), rows), f)
    readModCalls(f)
}

test_that("modCooccurrence(adjust = \"read\") removes an association made by per-read openness", {
    r1 <- GenomicRanges::GRanges("chr1:1001-1100"); r2 <- GenomicRanges::GRanges("chr1:3001-3100")
    ## per-read confound: "high" reads carry MSPs over much of the read and open each region in
    ## 80% of reads, "low" reads in 10%; within each kind the two regions are exactly independent
    n <- 400
    high <- rep(c(TRUE, FALSE), each = n / 2)
    k <- seq_len(n / 2)
    o1 <- c(k <= 160, k <= 20)
    o2 <- c(k <= 128 | (k > 160 & k <= 192), k <= 2 | (k > 20 & k <= 38))
    msps <- lapply(seq_len(n), function(i) {
        bg <- if (high[i]) seq(4000, 5800, by = 200) else c(4000, 5000)  # background MSPs outside the regions
        rbind(if (o1[i]) c(1000, 100), if (o2[i]) c(3000, 100), cbind(bg, if (high[i]) 150 else 60))
    })
    x <- sm_reads_table(msps)
    old <- modCooccurrence(x, r1, r2)
    expect_equal(old$adjust, "none")
    expect_null(old$adj_p_value)
    expect_lt(old$p_value, 0.01)
    expect_gt(old$odds_ratio, 1)
    lg <- modCooccurrence(x, r1, r2, adjust = "read")
    expect_identical(lg$table, old$table)
    expect_equal(lg$p_value, old$p_value)
    expect_gt(lg$adj_p_value, 0.05)
    expect_equal(lg$adjust_method, "logistic")
    expect_match(lg$covariate$how, "outside the two regions")
    cm <- modCooccurrence(x, r1, r2, adjust = "read", method = "cmh")
    expect_gt(cm$adj_p_value, 0.05)
    ## the "outside" covariate leaves out both regions: read 201 ("low", both regions open) has
    ## background MSPs of 2 x 60 bp over 6000 - 200 bp; the whole read adds the two 100-bp MSPs
    cv <- epiSeeker:::.sm_read_openness(x, "msp", exclude = c(r1, r2))
    expect_equal(unname(cv["r201"]), 120 / 5800)
    cv_all <- epiSeeker:::.sm_read_openness(x, "msp")
    expect_equal(unname(cv_all["r201"]), 320 / 6000)
    expect_gt(modCooccurrence(x, r1, r2, adjust = "read", covariate = "all")$adj_p_value, 0.05)
    ## a numeric covariate is used as given
    nm <- modCooccurrence(x, r1, r2, adjust = "read", covariate = stats::setNames(as.numeric(high), paste0("r", seq_len(n))))
    expect_gt(nm$adj_p_value, 0.05)
    expect_match(nm$covariate$how, "given by the user")
    expect_error(modCooccurrence(x, r1, r2, adjust = "read", covariate = 1:3), "named by read id")
})

test_that("modCooccurrence(adjust = \"read\") keeps a real same-molecule association", {
    r1 <- GenomicRanges::GRanges("chr1:1001-1100"); r2 <- GenomicRanges::GRanges("chr1:3001-3100")
    set.seed(7)
    n <- 200
    msps <- lapply(seq_len(n), function(i) {
        o1 <- runif(1) < 0.5
        o2 <- if (runif(1) < 0.85) o1 else !o1          # region 2 follows region 1 on the molecule
        nbg <- sample(1:8, 1)                             # openness elsewhere, unrelated to the pair
        rbind(if (o1) c(1000, 100), if (o2) c(3000, 100), cbind(4000 + 200 * seq_len(nbg), 100))
    })
    x <- sm_reads_table(msps)
    for (m in c("logistic", "cmh")) {
        co <- modCooccurrence(x, r1, r2, adjust = "read", method = m)
        expect_lt(co$adj_p_value, 0.01)
        expect_gt(co$adj_odds_ratio, 1)
    }
})
