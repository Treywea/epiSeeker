#' @title Env function for epiSeeker
#' @param TxDb TxDb object
#' @param item item name
#' @param force force to update txdb item in cache or not.
#' @importFrom yulab.utils get_cache_item
#' @importFrom yulab.utils update_cache_item
#' @importFrom yulab.utils rm_cache_item
#' @importFrom yulab.utils initial_cache_item
#' @importFrom S4Vectors metadata
#' @return Returns `invisible(NULL)` invisibly. The primary purpose of this function is to manage
#'   the TXDB cache through side effects (creating, updating, or removing cached objects),
#'   rather than returning a value.
.epiSeekerEnv <- function(TxDb, item = "epiSeekerEnv", force = FALSE) {
    # get cache item
    # it will create a list if there is no a cache item
    cache_item <- get_cache_item(item)

    # if there is no TXDB cached, write in cache
    if (is.null(cache_item$TXDB)) {
        update_cache_item(item = item, list(TXDB = TxDb))
        message(">> Using Genome:", get_env_genome(), "...\n")
        return(invisible(NULL))
    }

    # force to update item
    if (force) {
        message(">> Force to update txdb in cache...\n")
        rm_cache_item(item)
        initial_cache_item(item)
        update_cache_item(item, list(TXDB = TxDb))
        message(">> Using Genome:", get_env_genome(), "...\n")
        return(invisible(NULL))
    }

    # if exist TXDB
    TXDB <- cache_item$TXDB
    m1 <- tryCatch(unlist(metadata(TXDB)), error = function(e) NULL)
    m2 <- tryCatch(unlist(metadata(TxDb)), error = function(e) NULL)
    if (!is.null(m1)) {
        m1 <- m1[!is.na(m1)]
    }
    if (!is.null(m2)) {
        m2 <- m2[!is.na(m2)]
    }
    txdb_flag <- is.character(all.equal(TXDB, TxDb))

    if (
        is.null(m1) ||
            is.null(m2) ||
            length(m1) != length(m2) ||
            any(m1 != m2) ||
            txdb_flag
    ) {
        message(">> Update txdb in cache...\n")
        rm_cache_item(item)
        initial_cache_item(item)
        update_cache_item(item, list(TXDB = TxDb))
    }

    message(">> Using Genome:", get_env_genome(), "...\n")

    invisible(NULL)

    # pos <- 1
    # envir <- as.environment(pos)
    # if (!exists("epiSeekerEnv", envir=.GlobalEnv)) {
    #     assign("epiSeekerEnv", new.env(), envir = envir)
    # }

    # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)
    # if (!exists("TXDB", envir=epiSeekerEnv, inherits=FALSE)) {
    #     ## first run
    #     assign("TXDB", TxDb, envir=epiSeekerEnv)
    # } else {
    #     TXDB <- get("TXDB", envir=epiSeekerEnv)
    #     m1 <- tryCatch(unlist(metadata(TXDB)), error=function(e) NULL)

    #     m2 <- unlist(metadata(TxDb))

    #     if (!is.null(m1)) {
    #         m1 <- m1[!is.na(m1)]
    #     }
    #     m2 <- m2[!is.na(m2)]

    #     if ( is.null(m1) || length(m1) != length(m2) || any(m1 != m2) ) {
    #         rm(epiSeekerEnv)
    #         assign("epiSeekerEnv", new.env(), envir = envir)
    #         epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)
    #         assign("TXDB", TxDb, envir=epiSeekerEnv)
    #     }
    # }
}


#' @importFrom GenomicFeatures exonsBy
#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
get_exonList <- function(item = "epiSeekerEnv") {
    TxDb <- get_cache_element(item = item, elements = "TXDB")

    exonList <- get_cache_element(item = item, elements = "exonList")

    if (is.null(exonList)) {
        exonList <- exonsBy(TxDb)
        update_cache_item(item = item, list("exonList" = exonList))
    }

    # TxDb <- get("TXDB", envir=epiSeekerEnv)
    # if ( exists("exonList", envir=epiSeekerEnv, inherits=FALSE) ) {
    #     exonList <- get("exonList", envir=epiSeekerEnv)
    # } else {
    #     exonList <- exonsBy(TxDb)
    #     assign("exonList", exonList, envir=epiSeekerEnv)
    # }
    return(exonList)
}

#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
get_intronList <- function(item = "epiSeekerEnv") {
    TxDb <- get_cache_element(item = item, elements = "TXDB")

    intronList <- get_cache_element(item = item, elements = "intronList")

    if (is.null(intronList)) {
        intronList <- intronsByTranscript(TxDb)
        update_cache_item(item = item, list("intronList" = intronList))
    }

    # TxDb <- get("TXDB", envir=epiSeekerEnv)
    # if ( exists("intronList", envir=epiSeekerEnv, inherits=FALSE) ) {
    #     intronList <- get("intronList", envir=epiSeekerEnv)
    # } else {
    #     intronList <- intronsByTranscript(TxDb)
    #     assign("intronList", intronList, envir=epiSeekerEnv)
    # }
    return(intronList)
}


getCols <- function(n) {
    col <- c(
        "#8dd3c7",
        "#ffffb3",
        "#bebada",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
        "#d9d9d9",
        "#bc80bd",
        "#ccebc5",
        "#ffed6f"
    )

    col2 <- c(
        "#1f78b4",
        "#ffff33",
        "#c2a5cf",
        "#ff7f00",
        "#810f7c",
        "#a6cee3",
        "#006d2c",
        "#4d4d4d",
        "#8c510a",
        "#d73027",
        "#78c679",
        "#7f0000",
        "#41b6c4",
        "#e7298a",
        "#54278f"
    )

    col3 <- c(
        "#a6cee3",
        "#1f78b4",
        "#b2df8a",
        "#33a02c",
        "#fb9a99",
        "#e31a1c",
        "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#ffff99",
        "#b15928"
    )

    ## colorRampPalette(brewer.pal(12, "Set3"))(n)
    col3[seq_len(n)]
}

getPalette <- function(n) {
    palette <- c(
        "RdBu",
        "RdYlGn",
        "Spectral",
        "RdYlBu",
        "PiYG",
        "PRGn",
        "PuOr",
        "BrBG",
        "RdGy"
    )

    palette[seq_len(n)]
}

getSgn <- function(data, idx, statistic, missingDataAsZero) {
    d <- data[idx, ]
    # if(statistic == "mean"){
    #     ss <- colMeans(d)
    # } else if(statistic == "median"){
    #     ss <- apply(d, 2, median)
    # } else if(statistic == "min"){
    #     ss <- apply(d, 2, min)
    # } else if(statistic == "max"){
    #     ss <- apply(d, 2, max)
    # } else if(statistic == "sum"){
    #     ss <- colSums(d)
    # } else if(statistic == "std"){
    #     ss <- apply(d, 2, sd)
    # }

    if (statistic == "mean") {
        if (missingDataAsZero) {
            ss <- colMeans(d)
        } else {
            ss <- apply(d, 2, function(x) mean(x[x != 0]))
        }
    } else if (statistic == "median") {
        if (missingDataAsZero) {
            ss <- apply(d, 2, median)
        } else {
            ss <- apply(d, 2, function(x) median(x[x != 0]))
        }
    } else if (statistic == "min") {
        if (missingDataAsZero) {
            ss <- apply(d, 2, min)
        } else {
            ss <- apply(d, 2, function(x) min(x[x != 0]))
        }
    } else if (statistic == "max") {
        if (missingDataAsZero) {
            ss <- apply(d, 2, max)
        } else {
            ss <- apply(d, 2, function(x) max(x[x != 0]))
        }
    } else if (statistic == "sum") {
        ss <- colSums(d)
    } else if (statistic == "std") {
        if (missingDataAsZero) {
            ss <- apply(d, 2, sd)
        } else {
            ss <- apply(d, 2, function(x) sd(x[x != 0]))
        }
    }

    ss <- ss / sum(ss)
    ss[is.na(ss)] <- 0

    return(ss)
}
parseBootCiPerc <- function(bootCiPerc) {
    bootCiPerc <- bootCiPerc$percent
    tmp <- length(bootCiPerc)
    ciLo <- bootCiPerc[tmp - 1]
    ciUp <- bootCiPerc[tmp]
    return(c(ciLo, ciUp))
}

## estimate CI using bootstraping
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom parallel detectCores
getTagCiMatrix <- function(
  tagMatrix,
  conf = 0.95,
  resample = 500,
  ncpus = detectCores() - 1,
  statistic_method,
  missingDataAsZero = TRUE
) {
    RESAMPLE_TIME <- resample
    trackLen <- ncol(tagMatrix)

    if (Sys.info()[1] == "Windows") {
        tagMxBoot <- boot(
            data = tagMatrix,
            statistic = function(data, idx) {
                getSgn(data, idx, statistic_method, missingDataAsZero)
            },
            R = RESAMPLE_TIME
        )
    } else {
        tagMxBoot <- boot(
            data = tagMatrix,
            statistic = function(data, idx) {
                getSgn(data, idx, statistic_method, missingDataAsZero)
            },
            R = RESAMPLE_TIME,
            parallel = "multicore",
            ncpus = ncpus
        )
    }
    message(
        ">> Running bootstrapping for tag matrix...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"),
        "\n"
    )
    # tagMxBootCi <- sapply(seq_len(trackLen), function(i) {
    #                     bootCiToken <- boot.ci(tagMxBoot, type = "perc", index = i)
    #                     ## parse boot.ci results
    #                     return(parseBootCiPerc(bootCiToken))
    #                     }
    #                 )
    tagMxBootCi <- vapply(
        X = seq_len(trackLen),
        FUN = function(i) {
            bootCiToken <- boot.ci(tagMxBoot, type = "perc", index = i)
            parseBootCiPerc(bootCiToken)
        },
        FUN.VALUE = numeric(2)
    )
    row.names(tagMxBootCi) <- c("Lower", "Upper")
    return(tagMxBootCi)
}

#' @importFrom methods missingArg
getTagCount <- function(
  tagMatrix,
  xlim,
  statistic_method,
  conf,
  missingDataAsZero = TRUE,
  ...
) {
    statistic_method <- match.arg(
        statistic_method,
        c("mean", "median", "min", "max", "sum", "std")
    )

    if (statistic_method == "mean") {
        if (missingDataAsZero) {
            ss <- colMeans(tagMatrix)
        } else {
            ss <- apply(tagMatrix, 2, function(x) mean(x[x != 0]))
        }
    } else if (statistic_method == "median") {
        if (missingDataAsZero) {
            ss <- apply(tagMatrix, 2, median)
        } else {
            ss <- apply(tagMatrix, 2, function(x) median(x[x != 0]))
        }
    } else if (statistic_method == "min") {
        if (missingDataAsZero) {
            ss <- apply(tagMatrix, 2, min)
        } else {
            ss <- apply(tagMatrix, 2, function(x) min(x[x != 0]))
        }
    } else if (statistic_method == "max") {
        if (missingDataAsZero) {
            ss <- apply(tagMatrix, 2, max)
        } else {
            ss <- apply(tagMatrix, 2, function(x) max(x[x != 0]))
        }
    } else if (statistic_method == "sum") {
        ss <- colSums(tagMatrix)
    } else if (statistic_method == "std") {
        if (missingDataAsZero) {
            ss <- apply(tagMatrix, 2, sd)
        } else {
            ss <- apply(tagMatrix, 2, function(x) sd(x[x != 0]))
        }
    }

    ss <- ss / sum(ss)
    ss[is.na(ss)] <- 0

    ## plot(1:length(ss), ss, type="l", xlab=xlab, ylab=ylab)
    pos <- value <- NULL
    dd <- data.frame(pos = c(xlim[1]:xlim[2]), value = ss)
    if (!(missingArg(conf) || is.na(conf))) {
        tagCiMx <- getTagCiMatrix(
            tagMatrix,
            conf = conf,
            statistic_method = statistic_method,
            missingDataAsZero = missingDataAsZero,
            ...
        )
        dd$Lower <- tagCiMx["Lower", ]
        dd$Upper <- tagCiMx["Upper", ]
    }
    return(dd)
}


TXID2EG <- function(txid, geneIdOnly = FALSE) {
    txid <- as.character(txid)
    if (geneIdOnly == TRUE) {
        res <- TXID2EGID(txid)
    } else {
        res <- TXID2TXEG(txid)
    }
    return(res)
}

#' @importFrom GenomicFeatures transcripts
#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
TXID2TXEG <- function(txid) {
    txid2geneid <- get_cache_element(
        item = epiSeekerCache,
        elements = "txid2geneid"
    )

    if (is.null(txid2geneid)) {
        txdb <- get_cache_element(item = epiSeekerCache, elements = "TXDB")
        txidinfo <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
        idx <- which(vapply(txidinfo$gene_id, length, integer(1)) == 0)
        txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
        txid2geneid <- paste(
            mcols(txidinfo)[["tx_name"]],
            mcols(txidinfo)[["gene_id"]],
            sep = "/"
        )
        txid2geneid <- sub("/NA", "", txid2geneid)

        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        update_cache_item(item = epiSeekerCache, list("txid2geneid" = txid2geneid))
    }

    # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)

    # if (exists("txid2geneid", envir=epiSeekerEnv, inherits=FALSE)) {
    #     txid2geneid <- get("txid2geneid", envir=epiSeekerEnv)
    # } else {
    #     txdb <- get("TXDB", envir=epiSeekerEnv)
    #     txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    #     # idx <- which(sapply(txidinfo$gene_id, length) == 0)
    #     idx <- which(vapply(txidinfo$gene_id, FUN = length, FUN.VALUE = integer(1)) == 0)
    #     txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    #     txid2geneid <- paste(mcols(txidinfo)[["tx_name"]],
    #                          mcols(txidinfo)[["gene_id"]],
    #                          sep="/")
    #     txid2geneid <- sub("/NA", "", txid2geneid)

    #     names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    #     assign("txid2geneid", txid2geneid, envir=epiSeekerEnv)
    # }
    return(as.character(txid2geneid[txid]))
}

#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
TXID2EGID <- function(txid) {
    txid2geneid <- get_cache_element(item = epiSeekerCache, elements = "txid2eg")

    if (is.null(txid2geneid)) {
        txdb <- get_cache_element(item = epiSeekerCache, elements = "TXDB")
        txidinfo <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
        idx <- which(vapply(txidinfo$gene_id, length, integer(1)) == 0)
        txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
        txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])

        names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
        update_cache_item(item = epiSeekerCache, list("txid2eg" = txid2geneid))
    }

    # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)

    # if (exists("txid2eg", envir=epiSeekerEnv, inherits=FALSE)) {
    #     txid2geneid <- get("txid2eg", envir=epiSeekerEnv)
    # } else {
    #     txdb <- get("TXDB", envir=epiSeekerEnv)
    #     txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    #     # idx <- which(sapply(txidinfo$gene_id, length) == 0)
    #     idx <- which(vapply(txidinfo$gene_id, FUN = length, FUN.VALUE = integer(1)) == 0)
    #     txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    #     txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])

    #     names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    #     assign("txid2eg", txid2geneid, envir=epiSeekerEnv)
    # }
    return(as.character(txid2geneid[txid]))
}

## according to: https://support.bioconductor.org/p/70432/#70545
## contributed by Hervé Pagès
getFirstHitIndex <- function(x) {
    ## sapply(unique(x), function(i) which(x == i)[1])
    which(!duplicated(x))
}

#' calculate the overlap matrix, which is useful for vennplot
#'
#'
#' @title overlap
#' @param Sets a list of objects
#' @return data.frame
#' @author G Yu
overlap <- function(Sets) {
    ## this function is very generic.
    ## it call the getIntersectLength function to calculate
    ## the number of the intersection.
    ## if it fail, take a look at the object type were supported by getIntersectLength function.

    nn <- names(Sets)
    w <- t(apply(
        gtools::permutations(2, length(Sets), 0:1, repeats.allowed = TRUE),
        1,
        rev
    ))
    rs <- rowSums(w)
    wd <- as.data.frame(w)
    wd$n <- NA
    for (i in length(nn):0) {
        idx <- which(rs == i)
        if (i == length(nn)) {
            len <- getIntersectLength(Sets, as.logical(w[idx, ]))
            wd$n[idx] <- len
        } else if (i == 0) {
            wd$n[idx] <- 0
        } else {
            for (ii in idx) {
                ## print(ii)
                len <- getIntersectLength(Sets, as.logical(w[ii, ]))
                ww <- w[ii, ]
                jj <- which(ww == 0)
                pp <- gtools::permutations(2, length(jj), 0:1, repeats.allowed = TRUE)

                for (aa in 2:nrow(pp)) {
                    ## 1st row is all 0, abondoned
                    xx <- jj[as.logical(pp[aa, ])]
                    ww[xx] <- ww[xx] + 1
                    bb <- t(apply(w, 1, function(i) i == ww))
                    wd$n[rowSums(bb) == length(ww)]
                    ww <- w[ii, ]
                    len <- len - wd$n[rowSums(bb) == length(ww)]
                    ww <- w[ii, ]
                }
                wd$n[ii] <- len
            }
        }
    }
    colnames(wd) <- c(names(Sets), "Weight")
    return(wd)
}


getIntersectLength <- function(Sets, idx) {
    ## only use intersect and length methods in this function
    ## works fine with GRanges object
    ## and easy to extend to other objects.
    ss <- Sets[idx]
    ol <- ss[[1]]

    if (sum(idx) == 1) {
        return(length(ol))
    }

    for (j in 2:length(ss)) {
        ol <- intersect(ol, ss[[j]])
    }
    return(length(ol))
}

loadPeak <- function(peak, verbose = FALSE) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        if (verbose) {
            message(
                ">> loading peak file...\t\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"),
                "\n"
            )
        }
        peak.gr <- readPeakFile(peak, as = "GRanges")
    } else {
        stop("peak should be GRanges object or a peak file...")
    }
    return(peak.gr)
}

#' @title load defaulst txdb
#' @param TxDb txdb.
#' @return txdb object
loadTxDb <- function(TxDb) {
    rlang::check_installed(
        "TxDb.Hsapiens.UCSC.hg38.knownGene",
        reason = "Default txdb..."
    )

    if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        if (is.null(TxDb)) {
            warning(
                ">> TxDb is not specified, use 'TxDb.Hsapiens.UCSC.hg38.knownGene' by default..."
            )
            TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        }
        return(TxDb)
    }
}

#' @importFrom AnnotationDbi get
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicFeatures transcriptsBy
getGene <- function(TxDb, by = "gene") {
    .epiSeekerEnv(TxDb, item = epiSeekerCache)
    # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)

    by <- match.arg(by, c("gene", "transcript"))

    if (by == "gene") {
        features <- get_cache_element(item = epiSeekerCache, elements = "Genes")

        if (is.null(features)) {
            features <- genes(TxDb)
            update_cache_item(item = epiSeekerCache, list("Genes" = features))
        }

        # if ( exists("Genes", envir=epiSeekerEnv, inherits=FALSE) ) {
        #     features <- get("Genes", envir=epiSeekerEnv)
        # } else {
        #     features <- suppressMessages(genes(TxDb))
        #     assign("Genes", features, envir=epiSeekerEnv)
        # }
    } else {
        features <- get_cache_element(
            item = epiSeekerCache,
            elements = "Transcripts"
        )

        if (is.null(features)) {
            features <- transcriptsBy(TxDb)
            features <- unlist(features)
            update_cache_item(item = epiSeekerCache, list("Transcripts" = features))
        }

        # if ( exists("Transcripts", envir=epiSeekerEnv, inherits=FALSE) ) {
        #     features <- get("Transcripts", envir=epiSeekerEnv)
        # } else {
        #     features <- transcriptsBy(TxDb)
        #     features <- unlist(features)
        #     assign("Transcripts", features, envir=epiSeekerEnv)
        # }
    }

    return(features)
}


#' get filenames of sample files
#'
#'
#' @title getSampleFiles
#' @return list of file names
#' @examples
#' files <- getSampleFiles()
#' @export
#' @author G Yu
getSampleFiles <- function() {
    dir <- system.file("extdata", "GEO_sample_data", package = "epiSeeker")
    files <- list.files(dir)
    ## protein <- sub("GSM\\d+_", "", files)
    ## protein <- sub("_.+", "", protein)
    protein <- gsub(
        pattern = "GSM\\d+_(\\w+_\\w+)_.*",
        replacement = "\\1",
        files
    )
    protein <- sub("_Chip.+", "", protein)
    res <- paste(dir, files, sep = "/")
    res <- as.list(res)
    names(res) <- protein
    return(res)
}
## @importFrom RCurl getURL
## getDirListing <- function (url) {
##     ## from GEOquery
##     print(url)
##     a <- getURL(url)
##     b <- textConnection(a)
##     d <- read.table(b, header = FALSE)
##     close(b)
##     return(d)
## }

is.dir <- function(dir) {
    if (file.exists(dir) == FALSE) {
        return(FALSE)
    }
    return(file.info(dir)$isdir)
}


parse_targetPeak_Param <- function(targetPeak) {
    if (length(targetPeak) == 1) {
        if (is.dir(targetPeak)) {
            files <- list.files(path = targetPeak)
            idx <- unlist(lapply(c("bed", "bedGraph", "Peak"), grep, x = files))
            idx <- sort(unique(idx))
            files <- files[idx]
            targetPeak <- sub("/$", "", targetPeak)
            res <- paste(targetPeak, files, sep = "/")
        } else {
            if (!file.exists(targetPeak)) {
                stop("bed file is not exists...")
            } else {
                res <- targetPeak
            }
        }
    } else {
        if (is.dir(targetPeak[1])) {
            stop(
                "targetPeak should be a vector of bed file names or a folder containing bed files..."
            )
        } else {
            res <- targetPeak[file.exists(targetPeak)]
            if (length(res) == 0) {
                stop("targetPeak file not exists...")
            }
        }
    }
    return(res)
}

#' @importFrom S4Vectors metadata
IDType <- function(TxDb) {
    ##
    ## IDType <- metadata(TxDb)[8,2]
    ##
    ## update: 2015-10-27
    ## now IDType change from metadata(TxDb)[8,2] to metadata(TxDb)[9,2]
    ## it may change in future too
    ##
    ## it's safe to extract via grep

    md <- metadata(TxDb)
    md[grep("Type of Gene ID", md[, 1]), 2]
}

list_to_dataframe <- function(dataList) {
    if (is.null(names(dataList))) {
        return(do.call("rbind", dataList))
    }

    cn <- lapply(dataList, colnames) %>%
        unlist() %>%
        unique()
    cn <- c(".id", cn)
    dataList2 <- lapply(seq_along(dataList), function(i) {
        data <- dataList[[i]]
        data$.id <- names(dataList)[i]
        idx <- !cn %in% colnames(data)
        if (sum(idx) > 0) {
            for (i in cn[idx]) {
                data[, i] <- NA
            }
        }
        return(data[, cn])
    })
    res <- do.call("rbind", dataList2)
    res$.id <- factor(res$.id, levels = rev(names(dataList)))
    return(res)
}


## . function was from plyr package
#' capture name of variable
#'
#' @rdname dotFun
#' @export
#' @title .
#' @param ... expression
#' @param .env environment
#' @return expression
#' @examples
#' x <- 1
#' eval(.(x)[[1]])
. <- function(..., .env = parent.frame()) {
    structure(as.list(match.call()[-1]), env = .env, class = "quoted")
}


#' @title check upstream and downstream extension
#'
#' @param upstream upstream extension. One of actual number or rel() object.
#' @param downstream downstream extension. One of actual number or rel() object.
#' @param type one of "start_site", "end_site", "body".
#' @return message or null
#' @importFrom ggplot2 rel
check_extension <- function(upstream, downstream, type) {
    if (class(upstream) != class(downstream)) {
        stop("the type of upstream and downstream should be the same...")
    }

    # value of rel object should be in (0,1)
    if (inherits(upstream, "rel")) {
        if (type %in% c("start_site", "end_site")) {
            stop("Extension for start_site and end_site can not be rel object")
        }

        if (
            as.numeric(upstream) < 0 ||
                as.numeric(upstream) > 1 ||
                as.numeric(downstream) < 0 ||
                as.numeric(downstream) > 1
        ) {
            stop("the value of rel object should be in (0,1)...")
        }
    }

    # check actual number
    if (is.numeric(upstream)) {
        if (upstream < 0 || downstream < 0) {
            stop(
                "if upstream or downstream is integer, the value of it should be greater than 0..."
            )
        }
    }
}


#' @importFrom ggplot2 rel
#' @return function
#' @examples
#' rel(0.2)
#' @export
ggplot2::rel


#' @title check windows function
#' @param windows windows
#' @importFrom methods is
#' @return message or null
check_windows <- function(windows) {
    if (!is(windows, "GRanges")) {
        stop("windows should be a GRanges object...")
    }

    if (is.null(attr(windows, "type"))) {
        stop("windows should be made from getPromoters()/getBioRegion()")
    }

    # if (length(unique(width(windows))) != 1) {stop("width of windows should be equal...")}
}

#' @title  check bin parameter method
#'
#' @param nbin numbers of bin.
#' @param windows a list of region in granges.
#' @param verbose show details or not
#' @return message or nothing
check_bin <- function(nbin, windows, verbose) {
    type <- attr(windows, "type")

    if (type == "body" && is.null(nbin)) {
        stop("plotting body region should set the nbin parameter...")
    }

    if (!is.null(nbin)) {
        if (verbose) {
            message(
                ">> binning method is used...",
                format(Sys.time(), "%Y-%m-%d %X"),
                "\n",
                sep = ""
            )
        }

        if (!is.numeric(nbin)) {
            stop("nbin should be NULL or numeric...")
        }

        is.binning <- TRUE
    } else {
        is.binning <- FALSE
    }

    return(is.binning)
}

#' @title bin vector function
#'
#' @param vec vector.
#' @param nbin number of bin.
#' @return bin list
bin_vector <- function(vec, nbin = 800) {
    breaks <- seq(0, length(vec), length.out = nbin + 1)
    group_index <- findInterval(seq_along(vec), breaks, rightmost.closed = TRUE)
    bin_means <- tapply(vec, group_index, mean)

    return(bin_means)
}

#' @title change a list grange object to matrix
#' @param gr_list grange list object
#' @param weightCol weight column of peak.
#' @importFrom stats reshape
#' @return matrix
#' @examples
#' data(demo_peak)
#' grange2mt(list(a = demo_peak, b = demo_peak), "V5")
#' @export
grange2mt <- function(gr_list, weightCol = NULL) {
    df_list <- list()
    for (i in names(gr_list)) {
        tmp <- as.data.frame(gr_list[[i]])
        tmp$type <- i
        df_list[[i]] <- tmp
    }

    df_all <- do.call(rbind, df_list)
    df_all$name <- paste0(df_all$seqnames, ":", df_all$start, "-", df_all$end)

    if (is.null(weightCol)) {
        df_all <- df_all[, c("name", "type")]
        df_all$V5 <- 1
    } else {
        df_all <- df_all[, c("name", "type", weightCol)]
        colnames(df_all) <- c("name", "type", "V5")
    }

    df_wide <- reshape(
        df_all,
        timevar = "type",
        idvar = "name",
        direction = "wide"
    )
    rownames(df_wide) <- df_wide$name
    df_wide$name <- NULL
    colnames(df_wide) <- sub("V5\\.", "", colnames(df_wide))

    mat <- as.matrix(df_wide)
    mat[is.na(mat)] <- 0
    return(mat)
}

#' @title parse peak str
#' @param peak_str peak str
#' @return data frame
#' @examples
#' parse_peak("chr1:150235946-150236624")
#' @export
parse_peak <- function(peak_str) {
    parts <- strsplit(peak_str, ":|-")[[1]]
    data.frame(
        chr = parts[1],
        start = as.numeric(parts[2]),
        end = as.numeric(parts[3])
    )
}


make_Methylation_reference <- function(input, cover_depth) {
    ## make the methylation reference from bsseq object
    methylation_reference <- input@rowRanges

    ## add the methylation situation to methylation_reference
    for (i in input@colData@rownames) {
        methylation_M <- input@assays@data@listData[["M"]][, i]

        ## fill the 0 in Cov with 1
        methylation_cov <- input@assays@data@listData[["Cov"]][, i]

        if (cover_depth) {
            command <- paste0(
                "mcols(methylation_reference)$",
                i,
                "_depth <- methylation_cov"
            )

            eval(parse(text = command))
        }

        methylation_cov[input@assays@data@listData[["Cov"]][, i] == 0] <- 1

        ## the situation of methylation is calculated by M/cov
        methylation <- round(methylation_M / methylation_cov, 3)

        command <- paste0(
            "mcols(methylation_reference)$",
            i,
            "_methylation <- methylation"
        )

        eval(parse(text = command))
    }

    return(methylation_reference)
}


#' @importFrom SummarizedExperiment assayNames
make_reference <- function(input) {
    ## make the  reference from bmData object
    reference <- input@rowRanges

    ## extract the valueNames
    aName <- assayNames(input)

    n0 <- length(aName)

    if (n0 == 1) {
        for (i in input@colData@rownames) {
            value <- input@assays@data@listData[[aName]][, i]

            command <- paste0("mcols(reference)$", i, "_", aName, "<- value")

            eval(parse(text = command))
        }
    } else {
        for (i in input@colData@rownames) {
            aName1 <- aName[1]
            aName2 <- aName[2]

            value1 <- input@assays@data@listData[[aName1]][, i]
            value2 <- input@assays@data@listData[[aName2]][, i]

            command <- paste0("mcols(reference)$", i, "_", aName1, "<- value1")

            eval(parse(text = command))

            command <- paste0("mcols(reference)$", i, "_", aName2, "<- value2")

            eval(parse(text = command))
        }
    }

    return(reference)
}

#' @importFrom methods is
loadBSgenome <- function(BSgenome) {
    if (!is(BSgenome, "BSgenome")) {
        stop(">> input must be BSgenome object...")
    }

    if (is.null(BSgenome)) {
        stop(">> please specify BSgenome object...")
    } else {
        ## get the object from BSgenome
        BSgenome_name <- attr(BSgenome, "pkgname")

        ## no `character.only` (see getGeneAnno()), and read the exports of the namespace
        ## rather than ls("package:..."), which needs the package to be attached. The
        ## genome is the first export that is a BSgenome; the earlier pattern "^[^BSgenome]"
        ## is a character class and dropped organisms starting with B or S (Scerevisiae)
        if (requireNamespace(BSgenome_name, quietly = TRUE)) {
            for (obj in getNamespaceExports(BSgenome_name)) {
                x <- getExportedValue(BSgenome_name, obj)
                if (is(x, "BSgenome")) {
                    BSgenome <- x
                    break
                }
            }
        }
    }

    return(BSgenome)
}


#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges strand<-
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges mcols<-
#' @importFrom GenomicRanges start
#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges IRanges
detect_strand_and_motif <- function(
  region,
  motif,
  BSgenome,
  methylation_reference,
  input,
  base,
  position_bias
) {
    ## detect the strand and motif information
    dmRregion <- GRanges(
        seqnames = region$chr,
        ranges = IRanges(start = region$start, end = region$end)
    )

    hits <- findOverlaps(methylation_reference, dmRregion, type = "within")

    dmR_melth <- methylation_reference[hits@from]

    ## make the chromosome
    region$chr <- gsub("chr", "", region$chr, ignore.case = TRUE)
    region$chr <- as.numeric(region$chr)

    ## detect the positive and negative strand
    positive_index <- rep(FALSE, length(dmR_melth@ranges@start))

    for (i in motif) {
        melth_sequence <- DNAStringSet(
            BSgenome[[region$chr]],
            start = dmR_melth@ranges@start + position_bias[i] - 1,
            width = 1
        )

        tmp_positive_index <- vapply(
            melth_sequence,
            function(x) grepl(base, as.character(x)) == 1,
            FUN.VALUE = logical(1)
        )

        positive_index <- positive_index | tmp_positive_index
    }

    negetive_index <- !positive_index

    ## deal with the positive strand
    dmR_melth_positive <- dmR_melth[positive_index]
    strand(dmR_melth_positive) <- "+"

    for (i in motif) {
        ## create the regex for the motif
        regex <- create_regex_patterns_positive(i)

        ## count the length of the motif
        width <- length(strsplit(i, split = "")[[1]])

        ## create the mapping sequence
        sequence <- DNAStringSet(
            BSgenome[[region$chr]],
            start = dmR_melth@ranges@start[positive_index] - position_bias[i] + 1,
            width = width
        )

        ## get the index of the specific motif
        index <- vapply(
            sequence,
            function(x) grepl(regex, as.character(x)) == 1,
            FUN.VALUE = logical(1)
        )

        ## fill in the motif columns
        mcols(dmR_melth_positive)[["motif"]][index] <- i
    }

    ## deal with the negetive strand
    dmR_melth_negetive <- dmR_melth[negetive_index]
    strand(dmR_melth_negetive) <- "-"

    for (i in motif) {
        ## create the regex for the motif
        regex <- create_regex_patterns_negative(i)

        ## count the length of the motif
        width <- length(strsplit(i, split = "")[[1]])

        ## create the mapping sequence
        sequence <- DNAStringSet(
            BSgenome[[region$chr]],
            end = dmR_melth@ranges@start[negetive_index] + position_bias[i] - 1,
            width = width
        )

        ## get the index of the specific motif
        index <- vapply(
            sequence,
            function(x) grepl(regex, as.character(x)) == 1,
            FUN.VALUE = logical(1)
        )

        ## fill in the motif columns
        mcols(dmR_melth_negetive)[["motif"]][index] <- i
    }

    ## combine and reorder the results from positive and negative strand
    dmR_melth <- c(dmR_melth_positive, dmR_melth_negetive)
    dmR_melth <- dmR_melth[order(start(dmR_melth)), ]

    ## filtered the unidentified melthylation sites
    dmR_melth <- dmR_melth[which(!is.na(mcols(dmR_melth)[["motif"]]))]

    return(dmR_melth)
}


#' create regex patterns in positive strand
#'
#' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
#' @return regex pattern
create_regex_patterns_positive <- function(motif) {
    ## split the motif
    motif_list <- as.list(strsplit(motif, split = "")[[1]])

    ## check the character in the motif
    legal_character <- c(
        "R",
        "Y",
        "M",
        "K",
        "S",
        "W",
        "H",
        "B",
        "V",
        "D",
        "N",
        "A",
        "T",
        "G",
        "C"
    )

    if (!all(unlist(motif_list) %in% legal_character)) {
        stop("please input legal motif...")
    }

    ## make the connection between degenerate bases and normal bases
    converted_motif <- lapply(motif_list, function(x) {
        if (x == "R") {
            return("[AG]")
        }

        if (x == "Y") {
            return("[CT]")
        }

        if (x == "M") {
            return("[AC]")
        }

        if (x == "K") {
            return("[GT]")
        }

        if (x == "S") {
            return("[GC]")
        }

        if (x == "W") {
            return("[AT]")
        }

        if (x == "H") {
            return("[ATC]")
        }

        if (x == "B") {
            return("[GTC]")
        }

        if (x == "V") {
            return("[GAC]")
        }

        if (x == "D") {
            return("[GAT]")
        }

        if (x == "N") {
            return("[ATGC]")
        }

        if (x == "A") {
            return("A")
        }

        if (x == "T") {
            return("T")
        }

        if (x == "G") {
            return("G")
        }

        if (x == "C") {
            return("C")
        }
    })

    ## create the regex pattern
    regex <- paste0(unlist(converted_motif), collapse = "")

    return(regex)
}


#' create regex patterns in negative strand
#'
#' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
#' @return regex pattern
create_regex_patterns_negative <- function(motif) {
    ## split the motif
    motif_list <- rev(as.list(strsplit(motif, split = "")[[1]]))

    ## check the character in the motif
    legal_character <- c(
        "R",
        "Y",
        "M",
        "K",
        "S",
        "W",
        "H",
        "B",
        "V",
        "D",
        "N",
        "A",
        "T",
        "G",
        "C"
    )

    if (!all(unlist(motif_list) %in% legal_character)) {
        stop("please input legal motif...")
    }

    ## make the connection between degenerate bases and normal bases
    converted_motif <- lapply(motif_list, function(x) {
        ## R->A/G->T/C
        if (x == "R") {
            return("[TC]")
        }

        ## Y->C/T->G/A
        if (x == "Y") {
            return("[GA]")
        }

        ## M->A/C->T/G
        if (x == "M") {
            return("[TG]")
        }

        ## K->G/T->C/A
        if (x == "K") {
            return("[CA]")
        }

        ## S->G/C->C/G
        if (x == "S") {
            return("[CG]")
        }

        ## W->A/T->T/A
        if (x == "W") {
            return("[TA]")
        }

        ## H->A/T/C->T/A/G
        if (x == "H") {
            return("[TAG]")
        }

        ## B->G/T/C->C/A/G
        if (x == "B") {
            return("[CAG]")
        }

        ## V->G/A/C->C/T/G
        if (x == "V") {
            return("[CTG]")
        }

        ## D->G/A/T->C/T/A
        if (x == "D") {
            return("[CTA]")
        }

        ## N->A/T/C/G->T/A/G/C
        if (x == "N") {
            return("[TAGC]")
        }

        ## A->A->T
        if (x == "A") {
            return("T")
        }

        ## T->T->A
        if (x == "T") {
            return("A")
        }

        ## G->G->C
        if (x == "G") {
            return("C")
        }

        ## C->C->G
        if (x == "C") {
            return("G")
        }
    })

    ## create the regex pattern
    regex <- paste0(unlist(converted_motif), collapse = "")

    return(regex)
}


.check_valueNames <- function(valueNames, n0) {
    if (is.null(valueNames)) {
        valueNames <- paste0("value", seq_len(n0))
    }

    ## check the coordination of valueNames and value1/2
    if (length(valueNames) != n0) {
        stop("ValueNames do not match the value...")
    }

    return(valueNames)
}


.check_and_make_sampleNames <- function(data, sampleNames) {
    n0 <- length(data)
    if (!is.null(names(data)) && !any(is.na(names(data)))) {
        return(names(data))
    }

    if (is.null(sampleNames)) {
        sampleNames <- paste("sample", seq_len(n0), sep = "")
        return(sampleNames)
    }

    if (length(data) != length(sampleNames)) {
        stop("sampleNames should have equal length with data...")
    }

    return(sampleNames)
}


.check_variables_names <- function(data) {
    variables_names <- colnames(data[[1]])

    tmp <- unlist(lapply(data, function(x) {
        if (!identical(variables_names, colnames(x))) {
            return(TRUE)
        }

        return(FALSE)
    }))

    if (any(tmp)) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' @importFrom yulab.utils get_cache_item
get_env_genome <- function() {
    current_env <- get_cache_item(item = epiSeekerCache)

    env_txdb <- current_env$TXDB
    env_txdb_meta <- S4Vectors::metadata(env_txdb)
    env_txdb_version <- env_txdb_meta[grep("Genome", env_txdb_meta[, 1]), 2]

    return(env_txdb_version)
}
