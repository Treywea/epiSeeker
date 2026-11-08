library(epiSeeker)

context("test chromatin loops: readLoops, annotateLoops, linkPeaksToGenes")

## three genes on chr1 (TSS 10,000 +; 100,000 +; 310,000 -), as a self-defined TxDb
mock_genes <- function() {
    GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10000, 100000, 300000), c(20000, 110000, 310000)),
                           strand = c("+", "+", "-"), gene_id = c("G1", "G2", "G3"))
}

## L1 A-B (B at the G2 promoter), L2 C-D (no gene), L3 E-F (G1 and G3 promoters),
## L4 A-chr2 (between chromosomes), L5 G-B
mock_bedpe <- function() {
    data.frame(c1 = c("chr1", "chr1", "chr1", "chr1", "chr1"),
               s1 = c(40000, 150000, 8000, 40000, 115000), e1 = c(45000, 155000, 12000, 45000, 120000),
               c2 = c("chr1", "chr1", "chr1", "chr2", "chr1"),
               s2 = c(98000, 200000, 308000, 1000, 98000), e2 = c(103000, 205000, 312000, 5000, 103000),
               score = c(5, 4, 3, 2, 1))
}

test_that("readLoops reads BEDPE with a Juicer header, a plain header or none, and data frames", {
    df <- mock_bedpe()
    f1 <- tempfile(fileext = ".bedpe")
    writeLines(c("#chr1\tx1\tx2\tchr2\ty1\ty2\tscore", "# juicer_tools version 2",
                 apply(df, 1, paste, collapse = "\t")), f1)
    g1 <- readLoops(f1)
    expect_s4_class(g1, "GInteractions")
    expect_equal(length(g1), 5)
    a1 <- InteractionSet::anchors(g1, type = "first")
    expect_equal(GenomicRanges::start(a1)[1], 40001)          # 0-based BEDPE start + 1
    expect_equal(GenomicRanges::end(a1)[1], 45000)
    expect_equal(as.numeric(S4Vectors::mcols(g1)$score), df$score)
    expect_equal(S4Vectors::metadata(g1)$loopSource$source, basename(f1))

    f2 <- tempfile(fileext = ".bedpe.gz")
    con <- gzfile(f2, "w")
    writeLines(c("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcc", apply(df, 1, paste, collapse = "\t")), con)
    close(con)
    g2 <- readLoops(f2)
    expect_equal(length(g2), 5)
    expect_true("cc" %in% names(S4Vectors::mcols(g2)))

    f3 <- tempfile(fileext = ".bedpe")
    writeLines(apply(df[, 1:6], 1, paste, collapse = "\t"), f3)
    expect_equal(length(readLoops(f3)), 5)

    df2 <- df; df2$s2[2] <- NA
    g4 <- readLoops(df2)
    expect_equal(length(g4), 4)
    expect_equal(S4Vectors::metadata(g4)$loopSource$dropped, 1)
    expect_identical(readLoops(g4), g4)
    expect_error(readLoops(df[, 1:5]), "six columns")
    expect_error(readLoops("no/such/file.bedpe"), "does not exist")
})

test_that("annotateLoops marks promoter anchors with the tssRegion window and classes the loops", {
    la <- suppressWarnings(suppressMessages(annotateLoops(mock_bedpe(), TxDb = mock_genes(), verbose = FALSE)))
    expect_s4_class(la, "loopAnno")
    mc <- S4Vectors::mcols(la@loops)
    expect_equal(as.character(mc$loop_class),
                 c("promoter-distal", "distal-distal", "promoter-promoter", "distal-distal", "promoter-distal"))
    expect_equal(as.list(mc$anchor2_genes)[[1]], "G2")
    expect_setequal(c(unlist(mc$anchor1_genes[3]), unlist(mc$anchor2_genes[3])), c("G1", "G3"))  # minus-strand TSS
    expect_true(is.na(mc$span[4]))                             # between chromosomes
    expect_equal(as.numeric(mc$span[1]), 58000)
    expect_equal(sum(la@loopStat$n), 5)
    expect_equal(la@loopStat$n, c(1L, 2L, 2L))
    ## anchors are annotated by annotateSeq(): a csAnno, aligned by region_id
    expect_s4_class(la@anchors, "csAnno")
    expect_true(all(la@anchors@anno$region_id %in% seq_along(InteractionSet::regions(la@loops))))
    expect_output(show(la), "5 loops between")
    df <- as.data.frame(la)
    expect_equal(nrow(df), 5)
    expect_true(is.character(df$anchor2_genes))
    expect_s3_class(plotLoopAnno(la), "ggplot")
    expect_s3_class(plotLoopAnno(la, type = "span"), "ggplot")
    f <- tempfile(fileext = ".rds"); saveRDS(la, f)
    expect_identical(readRDS(f), la)
})

test_that("annotateLoops with gene_ids records only those genes and classes anchors by them", {
    la <- suppressWarnings(suppressMessages(annotateLoops(mock_bedpe(), TxDb = mock_genes(),
                                                          gene_ids = c("G1", "G2"), verbose = FALSE)))
    mc <- S4Vectors::mcols(la@loops)
    ## L3's G3 anchor is no longer a promoter anchor; L1 and L5 keep G2
    expect_equal(as.character(mc$loop_class),
                 c("promoter-distal", "distal-distal", "promoter-distal", "distal-distal", "promoter-distal"))
    expect_false(any(c(unlist(mc$anchor1_genes), unlist(mc$anchor2_genes)) == "G3"))
    expect_equal(as.list(mc$anchor2_genes)[[1]], "G2")
    expect_equal(S4Vectors::metadata(la@loops)$geneFilter$n_gene_ids, 2L)
    ## the default is unchanged
    la0 <- suppressWarnings(suppressMessages(annotateLoops(mock_bedpe(), TxDb = mock_genes(), verbose = FALSE)))
    expect_null(S4Vectors::metadata(la0@loops)$geneFilter)
    expect_equal(la0@loopStat$n, c(1L, 2L, 2L))
})

test_that("linkPeaksToGenes compares loop genes with the nearest gene", {
    g <- mock_genes()
    la <- suppressWarnings(suppressMessages(annotateLoops(mock_bedpe(), TxDb = g, verbose = FALSE)))
    pk <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(42000, 10500, 60000, 152000, 117000),
                                                         c(42500, 10800, 60500, 152500, 117300)))
    pa <- suppressMessages(annotateSeq(pk, TxDb = g, verbose = FALSE))
    links <- linkPeaksToGenes(pa, la)
    expect_equal(nrow(links), 5)
    expect_equal(as.character(links$link),
                 c("loop, other genes", "promoter", "no loop", "loop, no gene", "loop includes nearest"))
    expect_equal(links$nearest_gene[1], "G1")
    expect_equal(links$loop_genes[1], "G2")
    expect_equal(links$n_loops[1], 2L)                         # L1 and L4 share anchor A
    expect_equal(links$loop_genes[5], "G2")
    expect_s3_class(plotLoopLinks(links), "ggplot")
    expect_error(linkPeaksToGenes(pk, la), "must be a csAnno")
    expect_error(linkPeaksToGenes(pa, mock_bedpe()), "loopAnno")
})

test_that("annotateLoops on a TxDb records the genome and the gene symbols", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    ## one anchor at the TP53 promoter (chr17:7,687,490, minus strand)
    df <- data.frame("chr17", 7685000, 7690000, "chr17", 7740000, 7745000)
    la <- annotateLoops(df, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
    mc <- S4Vectors::mcols(la@loops)
    expect_true("7157" %in% unlist(mc$anchor1_genes))
    expect_true("TP53" %in% unlist(mc$anchor1_symbols))
    expect_equal(annoSource(la@anchors)$genome, "hg38")
    pa <- annotateSeq(GenomicRanges::GRanges("chr17", IRanges::IRanges(7741000, 7741500)), TxDb = txdb,
                      annoDb = "org.Hs.eg.db", verbose = FALSE)
    links <- linkPeaksToGenes(pa, la)
    expect_true(grepl("TP53", links$loop_symbols))
    expect_true("nearest_symbol" %in% names(links))
})
