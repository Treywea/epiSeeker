library(epiSeeker)

context("test function for loading BSgenome")

## in its own file: test-motif.R needs BSgenome.Hsapiens.UCSC.hg38 for every test
test_that("loadBSgenome returns the genome without the package being attached", {
    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    g <- getExportedValue("BSgenome.Scerevisiae.UCSC.sacCer3", "Scerevisiae")
    res <- epiSeeker:::loadBSgenome(g)
    expect_s4_class(res, "BSgenome")
    expect_equal(GenomeInfoDb::providerVersion(res), "sacCer3")
})
