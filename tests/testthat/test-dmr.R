library(epiSeeker)

context("readDMR(): DMR tables of DSS, methylKit and dmrseq to GRanges")

dss_table <- function() {
    data.frame(chr = c("chr1", "chr1", "chr2"), start = c(1001, 20001, 5001), end = c(1500, 20400, 5300),
               length = c(500, 400, 300), nCG = c(12, 9, 8), meanMethy1 = c(0.8, 0.7, 0.2),
               meanMethy2 = c(0.3, 0.4, 0.6), diff.Methy = c(0.5, 0.3, -0.4), areaStat = c(60, 30, -35))
}

test_that("DSS callDMR tables: two groups and multi-factor", {
    x <- readDMR(dss_table(), groups = c("tumor", "normal"), genome = "hg38")
    expect_s4_class(x, "GRanges")
    expect_equal(start(x), c(1001, 20001, 5001))
    expect_equal(end(x), c(1500, 20400, 5300))
    expect_equal(x$effect, c(0.5, 0.3, -0.4))
    expect_equal(x$direction, c("hyper", "hyper", "hypo"))
    expect_true(all(is.na(x$pvalue)) && all(is.na(x$qvalue)))
    expect_equal(x$n_CpG, c(12, 9, 8))
    expect_equal(unique(x$caller), "DSS")
    ## the caller's other columns stay, the mapped ones are not duplicated
    expect_true(all(c("meanMethy1", "meanMethy2", "areaStat") %in% names(mcols(x))))
    expect_false(any(c("diff.Methy", "nCG", "length") %in% names(mcols(x))))
    src <- S4Vectors::metadata(x)$dmrSource
    expect_equal(src$caller, "DSS")
    expect_equal(src$genome, "hg38")
    expect_equal(unname(src$groups), c("tumor", "normal"))
    expect_equal(src$n_direction[["hyper"]], 2L)
    expect_equal(src$n_direction[["hypo"]], 1L)
    expect_match(src$effect, "meanMethy1 - meanMethy2")
    ## multi-factor callDMR: no diff.Methy, direction from areaStat
    mf <- dss_table()[, c("chr", "start", "end", "length", "nCG", "areaStat")]
    y <- readDMR(mf)
    expect_true(all(is.na(y$effect)))
    expect_equal(y$direction, c("hyper", "hyper", "hypo"))
    expect_match(S4Vectors::metadata(y)$dmrSource$effect, "areaStat")
    expect_error(readDMR(dss_table(), groups = "tumor"), "two names")
})

test_that("methylKit: getData() table, written csv, methylDiff-like object", {
    mk <- data.frame(chr = c("chr1", "chr3"), start = c(352001, 9001), end = c(353000, 10000),
                     strand = c("*", "+"), pvalue = c(1e-6, 0.01), qvalue = c(1e-5, 0.04),
                     meth.diff = c(-16, 30))
    x <- readDMR(mk)
    expect_equal(unique(x$caller), "methylKit")
    expect_equal(x$effect, c(-0.16, 0.30))
    expect_equal(x$direction, c("hypo", "hyper"))
    expect_equal(x$pvalue, c(1e-6, 0.01))
    expect_equal(x$qvalue, c(1e-5, 0.04))
    expect_equal(as.character(strand(x)), c("*", "+"))
    f <- tempfile(fileext = ".csv")
    utils::write.csv(mk, f)  # with row names, quoted
    y <- readDMR(f)
    expect_equal(as.data.frame(mcols(y)), as.data.frame(mcols(x)))
    expect_equal(S4Vectors::metadata(y)$dmrSource$input, normalizePath(f))
    ## an S4 object extending data.frame, with the slots of methylKit's methylDiff
    setClass("methylDiff", representation(sample.ids = "character", assembly = "character",
                                          context = "character", treatment = "numeric",
                                          destranded = "logical", resolution = "character"),
             contains = "data.frame", where = environment())
    obj <- new("methylDiff", mk, sample.ids = c("c1", "t1"), assembly = "hg19", context = "CpG",
               treatment = c(0, 1), destranded = TRUE, resolution = "region")
    z <- readDMR(obj)
    expect_equal(z$effect, x$effect)
    src <- S4Vectors::metadata(z)$dmrSource
    expect_equal(src$genome, "hg19")
    expect_equal(src$resolution, "region")
    expect_equal(src$treatment, c(0, 1))
    obj@resolution <- "base"
    expect_message(readDMR(obj), "DMCs")
})

test_that("dmrseq: GRanges with and without meanDiff; a gzipped table with a #chr header", {
    g <- GenomicRanges::GRanges(c("chr21", "chr21"), IRanges::IRanges(c(100, 5000), c(900, 5600)),
                                L = c(10L, 7L), area = c(5, 3), beta = c(-0.4, 0.3), stat = c(-4, 3),
                                pval = c(0.001, 0.2), qval = c(0.01, 0.5),
                                index = IRanges::IRanges(c(1, 50), c(10, 56)))
    GenomeInfoDb::genome(g) <- "hg38"
    x <- readDMR(g)
    expect_equal(unique(x$caller), "dmrseq")
    expect_equal(x$effect, c(-0.4, 0.3))
    expect_equal(x$direction, c("hypo", "hyper"))
    expect_equal(x$qvalue, c(0.01, 0.5))
    expect_equal(x$n_CpG, c(10, 7))
    src <- S4Vectors::metadata(x)$dmrSource
    expect_equal(src$genome, "hg38")
    expect_match(src$effect, "not a methylation difference")
    g$meanDiff <- c(-0.2, 0.15)
    y <- readDMR(g)
    expect_equal(y$effect, c(-0.2, 0.15))
    expect_match(S4Vectors::metadata(y)$dmrSource$effect, "meanDiff")
    ## a written DSS table with "#chr", gzipped
    f <- tempfile(fileext = ".txt.gz")
    con <- gzfile(f, "w")
    d <- dss_table(); names(d)[1] <- "#chr"
    utils::write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE)
    close(con)
    expect_equal(readDMR(f)$effect, c(0.5, 0.3, -0.4))
    expect_error(readDMR(data.frame(chr = "chr1", start = 1, end = 2, foo = 1)), "cannot tell the DMR caller")
})

test_that("annotateSeq() keeps the DMR columns and provenance; the digest uses them", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    set.seed(1)
    eff <- round(runif(length(peak), -0.6, 0.6), 3)
    d <- data.frame(chr = as.character(seqnames(peak)), start = start(peak), end = end(peak),
                    length = width(peak), nCG = 5, meanMethy1 = 0.5, meanMethy2 = 0.5 - eff,
                    diff.Methy = eff, areaStat = eff * 100)
    x <- readDMR(d, groups = c("treated", "control"))
    anno <- annotateSeq(x, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
    expect_true(all(c("effect", "direction", "pvalue", "qvalue", "n_CpG", "caller") %in% names(mcols(anno@anno))))
    expect_equal(S4Vectors::metadata(anno@anno)$dmrSource$caller, "DSS")
    ev <- epiSeeker:::.seq_digest(anno, n_regions = 10)
    nh <- sum(anno@anno$direction == "hyper")
    expect_match(ev$sections$overview, paste0("DMRs: ", length(anno@anno), " from DSS"), fixed = TRUE)
    expect_match(ev$sections$overview, paste0("hyper ", nh, ", hypo ", length(anno@anno) - nh), fixed = TRUE)
    expect_match(ev$sections$overview, "(treated vs control)", fixed = TRUE)
    expect_match(ev$sections$top_regions, "ranked by |effect|", fixed = TRUE)
    expect_equal(ev$totals$rank_by, "|effect|")
    top <- ev$evidence$regions[1, ]
    i <- which(start(anno@anno) == top$start & as.character(seqnames(anno@anno)) == top$seqnames)
    expect_equal(abs(anno@anno$effect[i]), max(abs(anno@anno$effect)))
    ## an explicit signal still wins
    expect_match(epiSeeker:::.seq_digest(anno, signal = "areaStat")$sections$top_regions, "ranked by areaStat", fixed = TRUE)
    ## split by direction -> csAnnoList digest, one DMR line per group
    al <- csAnnoList(lapply(split(x, x$direction), annotateSeq, TxDb = txdb, verbose = FALSE))
    evl <- epiSeeker:::.seq_digest_list(al, n_regions = 3)
    expect_match(evl$sections$overview, "DMRs (hyper): ", fixed = TRUE)
    expect_match(evl$sections$top_regions, "ranked by |effect|", fixed = TRUE)
})
