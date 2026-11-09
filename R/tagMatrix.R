#' @title Get promoter region in GRanges format
#'
#' @param TxDb TxDb object
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param by one of 'gene', 'transcript'.
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' promoters <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000)
#' @return GRanges object
#' @author Guangchuang Yu
#' @export
getPromoters <- function(TxDb = NULL,
                         upstream = 1000,
                         downstream = 1000,
                         by = "gene") {
    getBioRegion(
        TxDb = TxDb,
        upstream = upstream,
        downstream = downstream,
        by = by,
        type = "start_site"
    )
}

#' @title Extend regions functions
#'
#' @param regions GRanges object
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR'.
#' @param type one of "start_site", "end_site", "body".
#' @return GRanges object
extend_gr <- function(regions, upstream, downstream, by, type) {
    check_extension(upstream = upstream, downstream = downstream, type = type)

    label <- switch(type,
        body = c("start_site", "end_site"),
        start_site = c("start_site"),
        end_site = c("end_site")
    )

    if (type == "body") {
        start_site <- start(regions)
        end_site <- end(regions)

        if (inherits(upstream, "rel")) {
            width_range <- abs(start_site - end_site)
            upstream_bp <- width_range * as.numeric(upstream)
            downstream_bp <- width_range * as.numeric(downstream)
        } else {
            upstream_bp <- upstream
            downstream_bp <- downstream
        }

        start_site <- ifelse(strand(regions) == "+", start_site - upstream_bp, start_site - downstream_bp)
        end_site <- ifelse(strand(regions) == "+", end_site + downstream_bp, end_site + upstream_bp)
    } else {
        if (type == "start_site") {
            coordinate <- ifelse(strand(regions) == "+", start(regions), end(regions))
        } else {
            coordinate <- ifelse(strand(regions) == "+", end(regions), start(regions))
        }

        # issue and code obtained from Chen Ting(NIH/NCI)
        start_site <- ifelse(strand(regions) == "+", coordinate - upstream, coordinate - downstream)
        end_site <- ifelse(strand(regions) == "+", coordinate + downstream, coordinate + upstream)
    }

    bioRegion <- GRanges(
        seqnames = seqnames(regions),
        ranges = IRanges(round(start_site), round(end_site)),
        strand = strand(regions)
    )
    bioRegion <- unique(bioRegion)

    # assign attribute
    attr(bioRegion, "type") <- type
    attr(bioRegion, "by") <- by
    attr(bioRegion, "label") <- label
    attr(bioRegion, "upstream") <- upstream
    attr(bioRegion, "downstream") <- downstream
    ## the chromosome lengths of the annotation the regions come from (the extended regions
    ## may pass the chromosome ends, so they are not set on bioRegion itself): getTagMatrix()
    ## checks the genome build of the peaks against them
    attr(bioRegion, "ref_seqinfo") <- GenomeInfoDb::seqinfo(regions)

    return(bioRegion)
}


#' @title Prepare a bioregion of selected feature
#'
#' @details this function combined previous functions getPromoters(), getBioRegion() and getGeneBody()
#' in order to solve the following issues.
#'
#' (1) <https://github.com/GuangchuangYu/ChIPseeker/issues/16>
#'
#' (2) <https://github.com/GuangchuangYu/ChIPseeker/issues/87>
#'
#' 1. function can provide a region of interest from txdb object.
#' 2. function can make region from granges object. txdb object do not contain insulator or enhancer regions.
#' Users can provide these regions through self-made granges object \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}.
#'
#' There are three kinds of way to extend regions: start_site, end_site and body.
#' We take transcript region to expain the differences of these three regions (tx: chr1 1000 1400).
#'
#' (1) body region refers to the 1000 ~ 1400 bp.
#'
#' (2) start_site region with (upstream = upstream = 100) refers to 900-1100bp.
#'
#' (3) end_site region with (upstream = upstream = 100) refers to 1300-1500bp.
#'
#' The regions keep the chromosome lengths (`seqinfo()`) of `TxDb` in `attr(, "ref_seqinfo")`;
#' `getTagMatrix()` checks the genome build of the peaks against them (see [annoChecks()]).
#' A self-made `GRanges` has lengths only when they were set (`seqlengths()`); without them
#' the build cannot be checked.
#'
#' @param TxDb TxDb object or self-made granges object.
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param by one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR', 'UTR'.
#' @param type one of "start_site", "end_site", "body".
#' @return GRanges object
#' @import IRanges GenomicRanges
#' @importFrom yulab.utils get_cache_item
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' getBioRegion(txdb)
#' @author Guangchuang Yu
#' @export
getBioRegion <- function(TxDb = NULL,
                         upstream = 1000,
                         downstream = 1000,
                         by = "gene",
                         type = "start_site") {
    type <- match.arg(type, c("start_site", "end_site", "body"))

    if (is(TxDb, "GRanges")) {
        message("#\n#.. 'TxDb' is a self-defined 'GRanges' object...\n#")

        windows <- extend_gr(
            regions = TxDb,
            upstream = upstream,
            downstream = downstream,
            by = by,
            type = type
        )
    } else {
        by <- match.arg(by, c("gene", "transcript", "exon", "intron", "3UTR", "5UTR", "UTR"))

        TxDb <- loadTxDb(TxDb)
        .epiSeekerEnv(TxDb, item = epiSeekerCache)
        # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)

        if (by == "gene" || by == "transcript") {
            regions <- getGene(TxDb, by)
        }

        if (by == "exon") {
            # exonList <- get_exonList(epiSeekerEnv)
            exonList <- get_exonList(item = epiSeekerCache)
            regions <- unlist(exonList)
        }

        if (by == "intron") {
            # intronList <- get_intronList(epiSeekerEnv)
            intronList <- get_intronList(item = epiSeekerCache)
            regions <- unlist(intronList)
        }

        if (by == "3UTR") {
            threeUTRList <- threeUTRsByTranscript(TxDb)
            regions <- unlist(threeUTRList)
        }

        if (by == "5UTR") {
            fiveUTRList <- fiveUTRsByTranscript(TxDb)
            regions <- unlist(fiveUTRList)
        }

        if (by == "UTR") {
            three_URT <- threeUTRsByTranscript(TxDb)
            three_UTR_regions <- unlist(three_URT)
            five_UTR <- fiveUTRsByTranscript(TxDb)
            five_UTR_regions <- unlist(five_UTR)
            regions <- c(three_UTR_regions, five_UTR_regions)
        }

        windows <- extend_gr(
            regions = regions,
            upstream = upstream,
            downstream = downstream,
            by = by,
            type = type
        )
    }

    return(windows)
}


#' @title getTagMatrix
#' @details getTagMatrix() function can produce the matrix for visualization.
#' Matrix represents the peak count in a windows and there are two ways to specify the 'windows':
#'
#' (1) use \code{\link{getPromoters}} and \code{\link{getBioRegion}} to get 'windows' and
#' put it into windows parameter in getTagMatrix().
#'
#' (2) use getTagMatrix() to call getPromoters()/getBioRegion().
#' In this way users do not need to input 'windows' parameter but need to input 'TxDb' parameter.
#' 'TxDb' can accept a set of packages contained annotation of regions of
#' different genomes(e.g. TxDb.Hsapiens.UCSC.hg38.knownGene).
#' Users can get the regions of interest through specific functions.
#' These specific functions are built in getPromoters()/getBioRegion().
#'
#' However, many regions can not be gain through txdb(e.g. insulator and enhancer regions),
#' Users can provide these regions in the form of granges object.
#' These self-made granges object will be passed to 'TxDb' and they will
#' be passed to makeBioRegionFromGranges() to produce the 'windows'.
#'
#' In a word, 'TxDb' parameter getTagMatrix() is a reference information.
#' Users can pass txdb object or self-made granges into it.
#'
#' @param peak (1) a peak file or GRanges object. (2) a list of peak file or GRanges object.
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param windows a collection of region
#' @param type one of "start_site", "end_site", "body"
#' @param by one of 'gene', 'transcript', 'exon', 'intron', '3UTR' , '5UTR', or specified by users
#' @param TxDb TxDb or self-made granges object, served as txdb
#' @param weightCol column name of weight, default is NULL.This column acts as a weight vaule.
#' Details see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}
#' @param nbin the amount of nbines. Calculate the tagMatrix by binning method.
#' Idea is derived from the function of deeptools(https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)
#' @param verbose print message or not
#' @param ignore_strand ignore the strand information or not
#' @examples
#' if (interactive()) {
#'     require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'     data(demo_peak)
#'     tagMatrix <- getTagMatrix(demo_peak,
#'         type = "start_site", by = "gene",
#'         upstream = 500, downstream = 500,
#'         TxDb = txdb, weightCol = "V7"
#'     )
#' }
#' @return tagMatrix
#' @author G Yu
#' @export
getTagMatrix <- function(peak,
                         upstream = 0,
                         downstream = 0,
                         windows = NULL,
                         type = NULL,
                         by = NULL,
                         TxDb = NULL,
                         weightCol = NULL,
                         nbin = NULL,
                         verbose = TRUE,
                         ignore_strand = FALSE) {
    if (is(peak, "list")) {
        if (is.null(names(peak))) {
            message("Automatically assign peaks name...\n")
            peak_names <- paste0("peak ", seq_len(length(peak)))
        } else {
            peak_names <- names(peak)
        }

        result <- lapply(peak, function(peak_item) {
            return(getTagMatrix.internal(
                peak = peak_item,
                upstream = upstream,
                downstream = downstream,
                windows = windows,
                type = type,
                by = by,
                TxDb = TxDb,
                weightCol = weightCol,
                nbin = nbin,
                verbose = verbose,
                ignore_strand = ignore_strand
            ))
        })


        names(result) <- peak_names
        attr(result[[1]], "multiple_samples") <- TRUE
    } else {
        result <- getTagMatrix.internal(
            peak = peak,
            upstream = upstream,
            downstream = downstream,
            windows = windows,
            type = type,
            by = by,
            TxDb = TxDb,
            weightCol = weightCol,
            nbin = nbin,
            verbose = verbose,
            ignore_strand = ignore_strand
        )
    }

    return(result)
}

#' @title getTagMatrix internal function
#'
#' @param peak peak file or GRanges object
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param windows a collection of region
#' @param type one of "start_site", "end_site", "body"
#' @param by one of 'gene', 'transcript', 'exon', 'intron', '3UTR' , '5UTR', or specified by users
#' @param TxDb TxDb or self-made granges object, served as txdb
#' @param weightCol column name of weight, default is NULL.
#' @param nbin the amount of nbines.
#' @param verbose print message or not
#' @param ignore_strand ignore the strand information or not
#' @return matrix
getTagMatrix.internal <- function(peak,
                                  upstream = 0,
                                  downstream = 0,
                                  windows = NULL,
                                  type = NULL,
                                  by = NULL,
                                  TxDb = NULL,
                                  weightCol = NULL,
                                  nbin = NULL,
                                  verbose = TRUE,
                                  ignore_strand = FALSE) {
    peak.gr <- loadPeak(peak)

    if (is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight = weight)
    }

    cov.len <- elementNROWS(peak.cov)
    cov.width <- GRanges(
        seqnames = names(cov.len),
        IRanges(
            start = rep(1, length(cov.len)),
            end = cov.len
        )
    )

    if (is.null(windows)) {
        if (length(by) == 1) {
            windows <- getBioRegion(
                TxDb = TxDb,
                upstream = upstream,
                downstream = downstream,
                by = by,
                type = type
            )

            windows <- subsetByOverlaps(windows, cov.width,
                type = "within", ignore.strand = FALSE
            )

            chr.idx <- intersect(
                names(peak.cov),
                unique(as.character(seqnames(windows)))
            )
        } else {
            intersect_chr.idx <- list()
            windows <- list()
            chr.idx <- list()
            if (length(type) != length(by)) {
                type <- rep(type, length(by))
            }

            for (i in seq_along(by)) {
                windows[[i]] <- getBioRegion(
                    TxDb = TxDb,
                    upstream = upstream,
                    downstream = downstream,
                    by = by[i],
                    type = type[i]
                )

                windows[[i]] <- subsetByOverlaps(windows[[i]], cov.width,
                    type = "within", ignore.strand = FALSE
                )

                intersect_chr.idx <- intersect(
                    names(peak.cov),
                    unique(as.character(seqnames(windows[[i]])))
                )

                chr.idx[[i]] <- unique(unlist(intersect_chr.idx))
            }
            names(windows) <- by
        }
    } else {
        if (is(windows, "list")) {
            raw_windows <- windows
            type <- vector(mode = "character")
            by <- vector(mode = "character")
            for (i in seq_len(length(windows))) {
                type[i] <- attr(raw_windows[[i]], "type")
                by[i] <- attr(raw_windows[[i]], "by")
            }


            intersect_chr.idx <- list()
            windows <- list()
            chr.idx <- list()

            for (i in seq_along(by)) {
                windows[[i]] <- subsetByOverlaps(raw_windows[[i]], cov.width,
                    type = "within", ignore.strand = FALSE
                )

                intersect_chr.idx <- intersect(
                    names(peak.cov),
                    unique(as.character(seqnames(raw_windows[[i]])))
                )

                chr.idx[[i]] <- unique(unlist(intersect_chr.idx))
            }

            names(windows) <- by
        } else {
            type <- attr(windows, "type")
            by <- attr(windows, "by")
            windows <- subsetByOverlaps(windows, cov.width,
                type = "within", ignore.strand = FALSE
            )

            chr.idx <- intersect(
                names(peak.cov),
                unique(as.character(seqnames(windows)))
            )
        }
    }


    ## genome build of the peaks against the annotation the regions come from (the TxDb
    ## path of annotateSeq() has the same check); recorded with every matrix
    if (is(windows, "list")) {
        build <- lapply(seq_along(windows), function(i) {
            .tagmatrix_build_check(peak.gr, windows[[i]], paste0("the regions '", by[i], "'"))
        })
    } else {
        build <- .tagmatrix_build_check(peak.gr, windows, "the regions")
    }

    if (verbose) {
        message(">> preparing tag matrix for ", type, " regions by ", by, "... ",
            format(Sys.time(), "%Y-%m-%d %X"), "\n",
            sep = " "
        )
    }


    if (is(windows, "list")) {
        result <- lapply(seq_len(length(windows)), function(idx) {
            check_windows(windows[[idx]])
            if (type[idx] == "body") {
                tmp <- getTagMatrix_body(
                    peak.cov = peak.cov,
                    windows = windows[[idx]],
                    nbin = nbin,
                    verbose = verbose,
                    ignore_strand = ignore_strand,
                    what = names(windows)[idx]
                )
            } else {
                tmp <- getTagMatrix_site(
                    peak.cov = peak.cov,
                    windows = windows[[idx]],
                    chr.idx = chr.idx[[idx]],
                    nbin = nbin,
                    verbose = verbose,
                    ignore_strand = ignore_strand
                )
            }

            attr(tmp, "checks") <- c(list(build = build[[idx]]), attr(tmp, "checks"))
            return(tmp)
        })

        names(result) <- names(windows)
    } else {
        check_windows(windows)

        if (type == "body") {
            result <- getTagMatrix_body(
                peak.cov = peak.cov,
                windows = windows,
                nbin = nbin,
                verbose = verbose,
                ignore_strand = ignore_strand
            )
        } else {
            result <- getTagMatrix_site(
                peak.cov = peak.cov,
                windows = windows,
                chr.idx = chr.idx,
                nbin = nbin,
                verbose = verbose,
                ignore_strand = ignore_strand
            )
        }
        attr(result, "checks") <- c(list(build = build), attr(result, "checks"))
    }

    message(">> done... ", format(Sys.time(), "%Y-%m-%d %X"), "\n", sep = "")

    return(result)
}

#' @title getTagMatrix function for region of site
#'
#' @param peak.cov peak coverage.
#' @param windows a collection of region.
#' @param chr.idx idx of chr.
#' @param nbin the amount of nbines
#' @param verbose print message or not
#' @param ignore_strand ignore the strand information or not
#' @return tagMatrix
#' @importFrom ggplot2 rel
#' @importFrom methods as
getTagMatrix_site <- function(peak.cov,
                              windows,
                              chr.idx,
                              nbin = NULL,
                              verbose = TRUE,
                              ignore_strand = FALSE) {
    # check nbin parameters
    is.binning <- check_bin(nbin, windows, verbose)

    peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
    tagMatrixList <- lapply(peakView, function(x) t(viewApply(x, as.vector)))

    tagMatrix <- do.call("rbind", tagMatrixList)

    # get the index of windows, that are reorganized by as(windows, "IntegerRangesList")
    idx.list <- split(seq_len(length(windows)), as.factor(seqnames(windows)))
    idx <- unlist(idx.list[chr.idx], use.names = FALSE)

    rownames(tagMatrix) <- idx
    tagMatrix <- tagMatrix[order(idx), ]

    # minus strand
    if (!ignore_strand) {
        minus.idx <- which(as.character(strand(windows)) == "-")
        tagMatrix[minus.idx, ] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
    }

    tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0, ]


    if (is.binning) {
        tagMatrix <- t(apply(tagMatrix, 1, bin_vector, nbin = nbin))
    }

    # assign attribute
    upstream <- attr(windows, "upstream")
    downstream <- attr(windows, "downstream")

    attr(tagMatrix, "upstream") <- upstream
    attr(tagMatrix, "downstream") <- downstream
    attr(tagMatrix, "label") <- c(
        paste0((-1) * upstream, "bp"),
        paste0((-0.5) * upstream, "bp"),
        attr(windows, "label"),
        paste0(0.5 * downstream, "bp"),
        paste0(downstream, "bp")
    )

    if (is.null(nbin)) {
        attr(tagMatrix, "label_position") <- c(
            1,
            1 + 0.5 * upstream,
            1 + upstream,
            1 + upstream + 0.5 * downstream,
            1 + upstream + downstream
        )
        attr(tagMatrix, "dash_pos") <- c(1 + upstream)
    } else {
        attr(tagMatrix, "label_position") <- round(seq(1, nbin, length.out = 5))
        attr(tagMatrix, "dash_pos") <- c(1 + nbin * 0.5)
    }

    attr(tagMatrix, "type") <- attr(windows, "type")
    attr(tagMatrix, "by") <- attr(windows, "by")


    return(tagMatrix)
}


#' @title getTagMatrix function for region of body
#'
#' @param peak.cov peak coverage.
#' @param windows a collection of region.
#' @param nbin the amount of nbines
#' @param verbose print message or not
#' @param ignore_strand ignore the strand information or not
#' @param what name of the set of regions, for the warnings
#' @return tagMatrix
#' @importFrom ggplot2 rel
getTagMatrix_body <- function(peak.cov,
                              windows,
                              nbin,
                              verbose = TRUE,
                              ignore_strand = FALSE,
                              what = NULL) {
    what <- if (is.null(what)) "" else paste0(" (regions '", what, "')")
    is.binning <- check_bin(nbin, windows, verbose)

    upstream <- attr(windows, "upstream")
    downstream <- attr(windows, "downstream")
    short_body <- NULL

    ## regions with any peak: only they can give a row (rows of zeros are removed below)
    n_regions <- length(windows)
    n_signal <- sum(.window_has_signal(peak.cov, windows))

    windows <- windows[width(windows) > nbin, ]
    ## windows without any peak give rows of zeros, which are removed below: they are left
    ## out here so that no intermediate matrix holds them. The row names keep indexing
    ## `windows` (after the width filter), as before.
    has_signal <- .window_has_signal(peak.cov, windows)
    row_ids <- which(has_signal)
    windows <- windows[has_signal]
    chr.idx <- intersect(
        names(peak.cov),
        unique(as.character(seqnames(windows)))
    )

    if (upstream == 0 && downstream == 0) {
        tagMatrix <- getTagMatrix_body_internal(
            peak.cov = peak.cov,
            windows = windows,
            nbin = nbin,
            chr.idx = chr.idx
        )

        label_position <- c(
            1,
            round(nbin * 0.25),
            round(nbin * 0.5),
            round(nbin * 0.75),
            nbin
        )
        label <- c(
            attr(windows, "label")[1],
            "25%",
            "50%",
            "75%",
            attr(windows, "label")[2]
        )
        dash_pos <- NULL
    } else if (inherits(upstream, "rel")) {
        upstream_bin <- nbin * as.numeric(upstream)
        downstream_bin <- nbin * as.numeric(downstream)

        tagMatrix <- getTagMatrix_body_internal(
            peak.cov = peak.cov,
            windows = windows,
            nbin = nbin + upstream_bin + downstream_bin,
            chr.idx = chr.idx
        )

        upstream_per <- as.numeric(upstream)
        downstream_per <- as.numeric(downstream)

        label_position <- c(
            1,
            round(nbin * upstream_per),
            round(nbin * upstream_per) + round(nbin * 0.25),
            round(nbin * upstream_per) + round(nbin * 0.5),
            round(nbin * upstream_per) + round(nbin * 0.75),
            round(nbin * upstream_per) + round(nbin),
            round(nbin * (upstream_per + downstream_per)) + round(nbin)
        )
        label <- c(
            paste0("-", upstream_per * 100, "%"),
            attr(windows, "label")[1],
            "25%",
            "50%",
            "75%",
            attr(windows, "label")[2],
            paste0(downstream_per * 100, "%")
        )

        dash_pos <- c(round(nbin * upstream_per), round(nbin * upstream_per) + round(nbin))
    } else {
        upstream_per <- downstream_per <- 0.15
        upstream_start <- start(windows)
        body_start <- upstream_end <- start(windows) + upstream
        downstream_end <- end(windows)
        body_end <- downstream_start <- end(windows) - downstream

        label_position <- c(
            1,
            round(nbin * upstream_per),
            round(nbin * upstream_per) + round(nbin * 0.25),
            round(nbin * upstream_per) + round(nbin * 0.5),
            round(nbin * upstream_per) + round(nbin * 0.75),
            round(nbin * upstream_per) + round(nbin),
            round(nbin * (upstream_per + downstream_per)) + round(nbin)
        )
        label <- c(
            paste0((-1) * upstream, "bp"),
            attr(windows, "label")[1],
            "25%",
            "50%",
            "75%",
            attr(windows, "label")[2],
            paste0(downstream, "bp")
        )
        dash_pos <- c(round(nbin * upstream_per), round(nbin * upstream_per) + round(nbin))

        mt_list <- list()

        body_windows <- GRanges(
            seqnames = seqnames(windows),
            ranges = IRanges(body_start, body_end),
            strand = strand(windows)
        )
        ## a body not longer than nbin gets zeros, while its flanks may not: such rows stay
        ## in the matrix with an empty body (counted below)
        short_body <- width(body_windows) <= nbin
        mt_list[["body_tagMatrix"]] <- getTagMatrix_body_internal(
            peak.cov = peak.cov,
            windows = body_windows,
            nbin = nbin,
            chr.idx = chr.idx
        )

        if (upstream == 0) {
            upstream_per <- 0
            mt_list[["upstream_tagMatrix"]] <- NULL
            label_position <- c(
                1,
                round(nbin * 0.25),
                round(nbin * 0.5),
                round(nbin * 0.75),
                round(nbin),
                round(nbin * downstream_per) + round(nbin)
            )
            label <- c(
                attr(windows, "label")[1],
                "25%",
                "50%",
                "75%",
                attr(windows, "label")[2],
                paste0(downstream, "bp")
            )
            dash_pos <- c(round(nbin))
        } else {
            upstream_windows <- GRanges(
                seqnames = seqnames(windows),
                ranges = IRanges(upstream_start, upstream_end),
                strand = strand(windows)
            )
            mt_list[["upstream_tagMatrix"]] <- getTagMatrix_body_internal(
                peak.cov = peak.cov,
                windows = upstream_windows,
                nbin = round(nbin * upstream_per),   # a fractional count gave flanks of zeros
                chr.idx = chr.idx
            )
        }

        if (downstream == 0) {
            downstream_per <- 0
            mt_list[["downstream_tagMatrix"]] <- NULL

            label_position <- c(
                1,
                round(nbin * upstream_per),
                round(nbin * upstream_per) + round(nbin * 0.25),
                round(nbin * upstream_per) + round(nbin * 0.5),
                round(nbin * upstream_per) + round(nbin * 0.75),
                round(nbin * upstream_per) + round(nbin)
            )
            label <- c(
                paste0((-1) * upstream, "bp"),
                attr(windows, "label")[1],
                "25%",
                "50%",
                "75%",
                attr(windows, "label")[2]
            )
            dash_pos <- c(round(nbin * upstream_per))
        } else {
            downstream_windows <- GRanges(
                seqnames = seqnames(windows),
                ranges = IRanges(downstream_start, downstream_end),
                strand = strand(windows)
            )

            mt_list[["downstream_tagMatrix"]] <- getTagMatrix_body_internal(
                peak.cov = peak.cov,
                windows = downstream_windows,
                nbin = round(nbin * downstream_per),   # a fractional count gave flanks of zeros
                chr.idx = chr.idx
            )
        }

        tagMatrix <- do.call("cbind", list(
            mt_list[["upstream_tagMatrix"]],
            mt_list[["body_tagMatrix"]],
            mt_list[["downstream_tagMatrix"]]
        ))
    }

    # minus strand
    if (!ignore_strand) {
        minus.idx <- which(as.character(strand(windows)) == "-")
        tagMatrix[minus.idx, ] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
    }

    keep <- rowSums(tagMatrix) != 0
    n_empty_body <- if (is.null(short_body)) 0L else sum(short_body & keep)
    tagMatrix <- tagMatrix[keep, ]
    if (!is.null(rownames(tagMatrix))) {
        rownames(tagMatrix) <- row_ids[as.integer(rownames(tagMatrix))]
    }

    ## Regions not longer than the number of bins cannot be split into them and are left
    ## out (most Alu elements with nbin = 800). That is not changed, but it is no longer
    ## silent: the number of regions with peaks that are not in the matrix is recorded and
    ## warned about.
    n_rows <- if (is.matrix(tagMatrix)) nrow(tagMatrix) else as.integer(length(tagMatrix) > 0)
    n_bins <- if (is.matrix(tagMatrix)) ncol(tagMatrix) else length(tagMatrix)
    short <- list(n_regions = n_regions, n_with_peaks = n_signal, n_rows = n_rows,
                  n_dropped = n_signal - n_rows, nbin = nbin, n_bins = n_bins,
                  fraction_dropped = if (n_signal > 0) (n_signal - n_rows) / n_signal else NA_real_,
                  n_empty_body = n_empty_body)
    if (short$n_dropped > 0) {
        why <- if (is.null(short_body)) {
            sprintf("they are not longer than the %d bins of type = \"body\" (nbin = %s%s)",
                    n_bins, nbin, if (n_bins != nbin) " plus the flanks" else "")
        } else {
            sprintf("their body is not longer than nbin = %s and their flanks have no peak", nbin)
        }
        warning(sprintf(paste0("%d of the %d regions with peaks (%.1f%%)%s are not in the tag matrix: %s, ",
                               "and they cannot be split into the bins. Use a smaller nbin, or a profile ",
                               "around a site (type = \"start_site\") for short regions such as Alu ",
                               "elements. Details: annoChecks(<tag matrix>)."),
                        short$n_dropped, n_signal, 100 * short$fraction_dropped, what, why), call. = FALSE)
    }
    if (n_empty_body > 0) {
        warning(sprintf(paste0("%d of the %d rows of the tag matrix%s have an empty body: the region body is ",
                               "not longer than nbin = %s and cannot be split into its bins, so only the ",
                               "flanks have values. Use a smaller nbin. Details: annoChecks(<tag matrix>)."),
                        n_empty_body, n_rows, what, nbin), call. = FALSE)
    }

    # assign attribute
    attr(tagMatrix, "upstream") <- upstream
    attr(tagMatrix, "downstream") <- downstream
    attr(tagMatrix, "type") <- attr(windows, "type")
    attr(tagMatrix, "by") <- attr(windows, "by")
    attr(tagMatrix, "label_position") <- label_position
    attr(tagMatrix, "label") <- label
    attr(tagMatrix, "dash_pos") <- dash_pos
    attr(tagMatrix, "checks") <- list(short_regions = short)

    return(tagMatrix)
}


#' @title get tagmatrix internal function
#'
#' @details The mean coverage in each of `nbin` bins of every window, with the bins of `bin_vector()`
#' (position i of a window of length L is in bin j when b_j <= i < b_(j+1), b = seq(0, L,
#' length.out = nbin + 1), the last bin includes L), computed as sums over the bins of the
#' run-length encoded coverage. The windows are never expanded to one value per base: that
#' took one integer per base pair of all windows (17-19 GB of memory for the genes of a human
#' genome with 20\% flanks). Windows not longer than `nbin` cannot be split into `nbin` bins
#' and get a row of zeros, as before.
#'
#' @param peak.cov peak coverage.
#' @param windows a collection of region.
#' @param nbin the amount of nbines.
#' @param chr.idx idx of chr.
#' @importFrom methods as
#' @return matrix
getTagMatrix_body_internal <- function(peak.cov, windows, nbin, chr.idx) {
    ## a fractional number of bins (e.g. nbin * 0.15) goes the old way: bin_vector() then
    ## returns ceiling(nbin) bins, which never equals nbin, and every row is zeros
    if (nbin != round(nbin)) {
        return(.getTagMatrix_body_internal_dense(peak.cov, windows, nbin, chr.idx))
    }
    nbin <- as.integer(round(nbin))
    idx.list <- split(seq_len(length(windows)), as.factor(seqnames(windows)))
    idx <- sort(unlist(idx.list[chr.idx], use.names = FALSE))
    tagMatrix <- matrix(0, nrow = length(idx), ncol = nbin,
                        dimnames = list(idx, seq_len(nbin)))   # the names bin_vector() gave
    w <- windows[idx]
    len <- width(w)
    ok <- which(len > nbin)
    chr <- as.character(seqnames(w))
    ## bins of up to ~5e6 at a time, so that the bin coordinates stay small too
    chunk <- max(1L, floor(5e6 / nbin))
    for (cc in unique(chr[ok])) {
        cv <- peak.cov[[cc]]
        if (is.integer(S4Vectors::runValue(cv))) S4Vectors::runValue(cv) <- as.numeric(S4Vectors::runValue(cv))
        rows <- ok[chr[ok] == cc]
        for (k in split(rows, ceiling(seq_along(rows) / chunk))) {
            L <- len[k]
            ## first position of bins 2..nbin, by window (rows) and bin (columns)
            b <- ceiling(outer(L / nbin, seq_len(nbin - 1L)))
            bs <- cbind(1, b)
            be <- cbind(b - 1, L)
            st <- t(bs + (start(w)[k] - 1))   # row m of bs is shifted by the start of window m
            wd <- t(be - bs + 1)
            v <- viewSums(Views(cv, IRanges(as.vector(st), width = as.vector(wd))))
            tagMatrix[k, ] <- matrix(v / as.vector(wd), ncol = nbin, byrow = TRUE)
        }
    }
    tagMatrix
}

## The earlier implementation, kept for a fractional number of bins: every window is
## expanded to one value per base (viewApply(as.vector)) and binned with bin_vector().
.getTagMatrix_body_internal_dense <- function(peak.cov, windows, nbin, chr.idx) {
    peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
    tagMatrixList <- lapply(peakView, function(x) {
        if (length((x)) == 1) {
            return(list(t(viewApply(x, as.vector))))
        } else {
            return(t(viewApply(x, as.vector)))
        }
    })

    if (length(unique(width(windows))) == 1) {
        tagMatrix <- do.call("rbind", tagMatrixList)

        tagMatrixList1 <- list()
        for (i in seq_len(nrow(tagMatrix))) {
            if (sum(tagMatrix[i, ]) != 0) {
                tmp <- bin_vector(tagMatrix[i, ], nbin = nbin)
                # extension out of genomic range
                if (length(tmp) != nbin) {
                    tagMatrixList1[[i]] <- vector(mode = "numeric", length = nbin)
                } else {
                    tagMatrixList1[[i]] <- tmp
                }
            } else {
                tagMatrixList1[[i]] <- vector(mode = "numeric", length = nbin)
            }
        }

        tagMatrix <- do.call("rbind", tagMatrixList1)
    } else {
        tagMatrixList1 <- list()
        for (i in seq_along(tagMatrixList)) {
            tagMatrixList1[[i]] <- lapply(tagMatrixList[[i]], function(x) {
                if (sum(x) != 0) {
                    tmp <- bin_vector(x, nbin = nbin)
                    # extension out of genomic range
                    if (length(tmp) != nbin) {
                        return(vector(mode = "numeric", length = nbin))
                    } else {
                        return(tmp)
                    }
                } else {
                    return(vector(mode = "numeric", length = nbin))
                }
            })
        }

        tagMatrix_list <- lapply(tagMatrixList1, function(x) {
            return(do.call("rbind", x))
        })
        tagMatrix <- do.call("rbind", tagMatrix_list)
    }

    # get the index of windows, that are reorganized by as(windows, "IntegerRangesList")
    idx.list <- split(seq_len(length(windows)), as.factor(seqnames(windows)))
    idx <- unlist(idx.list[chr.idx], use.names = FALSE)

    rownames(tagMatrix) <- idx

    tagMatrix <- tagMatrix[order(idx), ]

    return(tagMatrix)
}


## The genome-build check of getTagMatrix(): the peaks against the chromosome lengths of the
## annotation the windows were made from (attr "ref_seqinfo", set by getBioRegion()), or of
## the windows themselves for windows made otherwise. Warns or notes as annotateSeq() does.
.tagmatrix_build_check <- function(peak.gr, windows, ref_label) {
    if (!length(peak.gr)) return(NULL)
    si <- attr(windows, "ref_seqinfo")
    if (is.null(si)) si <- GenomeInfoDb::seqinfo(windows)
    b <- .build_check(peak.gr, windows, ref_label = ref_label, ref_seqinfo = si)
    .report_build(b, ref_label)
    b
}

## TRUE for every window that overlaps any coverage (a peak), FALSE otherwise (also for
## windows on chromosomes without coverage)
.window_has_signal <- function(peak.cov, windows) {
    out <- logical(length(windows))
    if (!length(windows)) return(out)
    chr <- as.character(seqnames(windows))
    for (cc in intersect(unique(chr), names(peak.cov))) {
        i <- which(chr == cc)
        cv <- peak.cov[[cc]]
        r <- ranges(windows[i])
        ## the part of each window inside the coverage vector (it ends at the last peak)
        e <- pmin(end(r), length(cv))
        s0 <- pmax(start(r), 1L)
        inside <- s0 <= e
        if (any(inside)) {
            out[i[inside]] <- viewSums(Views(cv, IRanges(s0[inside], e[inside]))) != 0
        }
    }
    out
}
