library(epiSeeker)

context("test adapters from single-cell and spatial data to csAnnoList")

## fragments of 3 barcodes in two groups; barcodes carry the -1 suffix of fragments files
toy_fragments <- function() {
    GenomicRanges::GRanges(
        c(rep("chr1", 8), rep("chr2", 4)),
        IRanges::IRanges(c(101, 151, 201, 251, 5101, 5151, 5201, 5251, 101, 151, 201, 251), width = 50),
        barcode = c("a-1", "a-1", "b-1", "b-1", "c-1", "c-1", "c-1", "c-1", "a-1", "b-1", "c-1", "c-1"),
        count = 1L
    )
}

test_that("readFragments reads 0-based fragments and keeps barcode and count", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("# comment", "chr1\t100\t300\tAAAC-1\t1", "chr2\t500\t700\tAAAG-1\t2"), f)
    fr <- readFragments(f)
    expect_s4_class(fr, "GRanges")
    expect_equal(BiocGenerics::start(fr), c(101, 501))
    expect_equal(fr$barcode, c("AAAC-1", "AAAG-1"))
    expect_equal(fr$count, c(1L, 2L))
    expect_equal(length(readFragments(f, barcodes = "AAAG-1")), 1)
    expect_equal(length(readFragments(f, which = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000)))), 1)
})

test_that("fragmentsToGroups matches barcodes without the suffix and keeps group-specific tiles", {
    groups <- c(a = "A", b = "A", c = "B")                   # table without the -1 suffix
    expect_message(res <- fragmentsToGroups(toy_fragments(), groups, tile = 1000, min_fragments = 1),
                   "Matched 3 of 3 barcodes of `groups` to fragments \\(without the -N suffix\\)")
    expect_s4_class(res, "GRangesList")
    expect_equal(names(res), c("A", "B"))
    expect_equal(S4Vectors::metadata(res)$members, c(A = 2L, B = 1L))
    ## group A: chr1:1-1000 (4 fragments) is A-specific; chr1:5001-6000 belongs to B
    expect_true("chr1" %in% as.character(GenomicRanges::seqnames(res$A)))
    expect_equal(BiocGenerics::start(res$B)[as.character(GenomicRanges::seqnames(res$B)) == "chr1"], 5001)
    ## same input as a data.frame
    df <- data.frame(barcode = names(groups), cluster = unname(groups))
    expect_equal(suppressMessages(fragmentsToGroups(toy_fragments(), df, tile = 1000, min_fragments = 1)), res)
    ## without the specificity filter every tile with enough fragments is kept
    all_tiles <- suppressMessages(fragmentsToGroups(toy_fragments(), groups, tile = 1000, min_fragments = 1,
                                                    specific = FALSE))
    expect_true(all(lengths(all_tiles) >= lengths(res)))
    expect_error(suppressMessages(fragmentsToGroups(toy_fragments(), c(zz = "A"))), "no fragment barcode")
})

test_that("peaksToGroups reads ArchR-, Signac- and 10x-style marker tables", {
    archr <- data.frame(seqnames = "chr1", start = c(100, 500), end = c(200, 600), group_name = c("C1", "C2"))
    expect_equal(names(peaksToGroups(archr)), c("C1", "C2"))
    signac <- data.frame(avg_log2FC = c(2, 3), cluster = c("1", "1"), row.names = c("chr1-100-200", "chr2-5-60"))
    g <- peaksToGroups(signac)
    expect_equal(as.character(GenomicRanges::seqnames(g[["1"]])), c("chr1", "chr2"))
    expect_equal(g[["1"]]$avg_log2FC, c(2, 3))
    tenx <- data.frame(feature = c("chr1:100-200", "chr3:7-90"), Cluster = c("2", "3"))
    expect_equal(lengths(peaksToGroups(tenx)), c(`2` = 1L, `3` = 1L))
    lst <- list(x = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10)))
    expect_s4_class(peaksToGroups(lst), "GRangesList")
    expect_error(peaksToGroups(data.frame(a = 1, cluster = 1)), "coordinates")
    ## the wide differential table of Cell Ranger ATAC
    wide <- data.frame(`Feature ID` = c("chr1:100-200", "chr1:300-400", "chr2:5-50"),
                       `Cluster 1 Log2 fold change` = c(2, 0.5, 3), `Cluster 1 Adjusted p value` = c(0.01, 0.01, 0.2),
                       `Cluster 2 Log2 fold change` = c(1.5, 1.2, -1), `Cluster 2 Adjusted p value` = c(0.01, 0.04, 0.01),
                       check.names = FALSE)
    cr <- peaksToGroups(wide)
    expect_equal(names(cr), c("1", "2"))
    expect_equal(lengths(cr), c(`1` = 1L, `2` = 2L))            # chr2 fails padj in 1; chr1:300 fails log2fc in 1
    expect_equal(cr[["2"]]$log2fc, c(1.5, 1.2))
})

test_that("annotateGroups returns a csAnnoList that records what the groups are", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    x <- GenomicRanges::GRangesList(d1 = peak[1:60], d2 = peak[61:150])
    S4Vectors::metadata(x)$members <- c(d1 = 812, d2 = 640)
    coords <- data.frame(barcode = c("s1", "s2", "s3"), x = 1:3, y = 3:1, group = c("d1", "d2", "d2"))
    res <- annotateGroups(x, TxDb = txdb, resolution = "spatial", coords = coords)
    expect_s4_class(res, "csAnnoList")
    gi <- groupInfo(res)
    expect_equal(gi$group, c("d1", "d2"))
    expect_equal(gi$resolution, c("spatial", "spatial"))
    expect_equal(gi$n_members, c(812, 640))
    expect_equal(gi$n_regions, c(60L, 90L))
    expect_equal(nrow(S4Vectors::metadata(res$d2@anno)$group$coords), 2)
    ## works wherever a csAnnoList does, and survives saving
    expect_s3_class(plotAnnoBar(res), "ggplot")
    f <- tempfile(fileext = ".rds"); saveRDS(res, f)
    expect_identical(groupInfo(readRDS(f)), gi)
    ## bulk objects report nothing
    data(peakAnno, package = "epiSeeker")
    expect_true(is.na(groupInfo(peakAnno)$resolution))
})

test_that("interpret() states the resolution and group sizes", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
    x <- GenomicRanges::GRangesList(c1 = peak[1:60], c2 = peak[61:150])
    res <- annotateGroups(x, TxDb = txdb, resolution = "single-cell", members = c(c1 = 120, c2 = 2000))
    ev <- epiSeeker:::.seq_digest_list(res, n_regions = 5)
    expect_match(ev$sections$overview, "Resolution: single-cell; pseudo-bulk of c1 120 cells, c2 2,000 cells",
                 fixed = TRUE)
    one <- epiSeeker:::.seq_digest(res$c1, n_regions = 5)
    expect_match(one$sections$overview, "Resolution: single-cell; pseudo-bulk of 120 cells", fixed = TRUE)
    data(peakAnno, package = "epiSeeker")
    expect_false(grepl("Resolution:", epiSeeker:::.seq_digest(peakAnno, n_regions = 5)$sections$overview))
})

test_that("spot positions are read and plotted", {
    f <- tempfile(fileext = ".csv")
    writeLines(c("AAAC,1,0,0,100,200", "AAAG,1,0,1,100,210", "AAAT,0,1,0,110,200"), f)
    pos <- readSpotPositions(f)
    expect_equal(names(pos), c("barcode", "in_tissue", "row", "col", "x", "y"))
    expect_equal(pos$x, c(200, 210, 200))
    g <- c("AAAC-1" = "A", "AAAG-1" = "B", "AAAT-1" = "B")    # groups named like fragments barcodes
    p <- plotSpatialGroups(pos, groups = g)
    expect_s3_class(p, "ggplot")
    expect_equal(ggplot2::ggplot_build(p)$data[[1]]$colour[1] != ggplot2::ggplot_build(p)$data[[1]]$colour[2], TRUE)
    expect_s3_class(plotSpatialGroups(pos, value = c(AAAC = 1, AAAG = 5)), "ggplot")
    hdr <- tempfile(fileext = ".csv")
    writeLines(c("barcode,in_tissue,array_row,array_col", "AAAC,1,0,0"), hdr)
    expect_equal(readSpotPositions(hdr)$x, 0)
})
