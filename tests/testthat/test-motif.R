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

test_that("getMotifMatrix places the hits of several regions at their own coordinates (cases 2025)", {
    data(pwm_obj, package = "epiSeeker")
    two <- GRanges("chr22", IRanges(c(23963111, 10525891), width = c(201, 101)))
    pwm <- pwm_obj[c(1:50, 120, 170)]
    key <- function(d) {
        d <- as.data.frame(d)
        if (!nrow(d)) return(character(0))
        sort(paste(d$chr, d$coordinate, d$score, d$strand, d$motif))
    }
    one <- lapply(1:2, function(i) getMotifMatrix(two[i], pwm, BSgenome.Hsapiens.UCSC.hg38))
    expect_gt(nrow(one[[1]]), 0)
    expect_gt(nrow(one[[2]]), 0)
    both <- getMotifMatrix(two, pwm, BSgenome.Hsapiens.UCSC.hg38)
    expect_equal(key(both), sort(c(key(one[[1]]), key(one[[2]]))))
    for (i in 1:2) {
        d <- both[both$region == i, ]
        expect_gt(nrow(d), 0)
        expect_true(all(d$coordinate >= start(two[i]) & d$coordinate <= end(two[i])))
    }
    expect_equal(attr(both, "range"), c(10525891L, 23963311L))
    expect_equal(attr(one[[1]], "range"), c(23963111L, 23963311L))
    expect_message(plotMotifProf(both), "2 regions")
    other <- both
    other$chr[other$region == 2] <- "chr21"
    expect_error(plotMotifProf(other), "several chromosomes")
})

test_that("getMotifMatrix maps the hits of a minus-strand region back to the genome", {
    data(pwm_obj, package = "epiSeeker")
    ## motifs scanned on both strands, so that the two orientations of a region hold the same sites
    pwm <- lapply(universalmotif::convert_motifs(pwm_obj[1:20]), function(m) { m@strand <- "+-"; m })
    r <- GRanges("chr22", IRanges(23963111, width = 201))
    rm <- r
    strand(rm) <- "-"
    span <- function(d) {
        d <- as.data.frame(d)
        a <- stats::aggregate(coordinate ~ seg_id + motif + strand + score, d, range)
        sort(paste(a$motif, a$coordinate[, 1], a$coordinate[, 2], a$strand, round(abs(a$score), 6)))
    }
    p <- getMotifMatrix(r, pwm, BSgenome.Hsapiens.UCSC.hg38)
    m <- getMotifMatrix(rm, pwm, BSgenome.Hsapiens.UCSC.hg38)
    expect_gt(nrow(p), 0)
    expect_equal(span(m), span(p))
    expect_true(all(m$coordinate >= start(r) & m$coordinate <= end(r)))
})
