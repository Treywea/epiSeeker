library(epiSeeker)
library(BSgenome.Hsapiens.UCSC.hg38)
library(universalmotif)

context("test function for motif analysis")

test_that("getMotifMatrix works on example region", {
    data(pwm_obj, package = "epiSeeker")

    region <- GRanges(
        seqnames = "chr22",
        ranges   = IRanges(start = 10525891, end = 10525991)
    )

    motif_df <- getMotifMatrix(
        region  = region,
        pwm     = pwm_obj[c(45, 120, 170)],
        ref_obj = BSgenome.Hsapiens.UCSC.hg38,
        by      = "name"
    )

    expect_true(is.null(motif_df) || is.data.frame(motif_df))

    expect_true(all(c("chr", "coordinate", "score", "strand", "motif") %in% names(motif_df)))

    attr_range <- attr(motif_df, "range")
    expect_equal(attr_range, c(10525891, 10525991))

    expect_type(motif_df$score, "double")


    expect_true(all(motif_df$strand %in% c("+", "-")))
})


test_that("plotMotifProf works with motifMatrix example", {
    data(pwm_obj, package = "epiSeeker")

    region <- GRanges(
        seqnames = "chr22",
        ranges   = IRanges(start = 10525891, end = 10525991)
    )

    motif_df <- getMotifMatrix(
        region  = region,
        pwm     = pwm_obj[c(45, 120, 170)],
        ref_obj = BSgenome.Hsapiens.UCSC.hg38
    )

    expect_true(is.null(motif_df) || is.data.frame(motif_df))

    p1 <- plotMotifProf(motif_df)

    expect_s3_class(p1, "ggplot")

    expect_true("coordinate" %in% names(motif_df))
    expect_true("score" %in% names(motif_df))
    expect_true("motif" %in% names(motif_df))

    expect_equal(attr(motif_df, "range"), c(10525891, 10525991))
})
