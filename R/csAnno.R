#' Class "csAnno"
#' This class represents the output of epiSeeker Annotation
#'
#' @name csAnno-class
#' @aliases csAnno-class show,csAnno-method vennpie,csAnno-method plotDistToTSS,csAnno-method plotAnnoBar,csAnno-method plotAnnoPie,csAnno-method upsetplot,csAnno-method subset,csAnno-method
#'
#' @docType class
#' @slot anno annotation
#' @slot tssRegion TSS region
#' @slot level transcript or gene
#' @slot hasGenomicAnnotation logical
#' @slot detailGenomicAnnotation Genomic Annotation in detail
#' @slot annoStat annotation statistics
#' @slot peakNum number of peaks
#' @exportClass csAnno
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
#' @seealso [annotateSeq()]
#' @return annotation object
#' @keywords classes
setClass("csAnno",
    representation = representation(
        anno = "GRanges",
        tssRegion = "numeric",
        level = "character",
        hasGenomicAnnotation = "logical",
        detailGenomicAnnotation = "data.frame",
        annoStat = "data.frame",
        peakNum = "numeric"
    )
)


#' convert csAnno object to GRanges
#'
#' @title as.GRanges
#' @param x csAnno object
#' @return GRanges object
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
#' @examples
#' data(peakAnno)
#' as.GRanges(peakAnno)
#' @export
as.GRanges <- function(x) {
    if (!is(x, "csAnno")) {
        stop("not supported...")
    }
    return(x@anno)
}

#' getting status of annotation
#'
#' @title getAnnoStat
#' @param x csAnno object
#' @return data frame
#' @examples
#' data(peakAnno)
#' getAnnoStat(peakAnno)
#' @export
getAnnoStat <- function(x) {
    if (!is(x, "csAnno")) {
        stop("not supported...")
    }
    return(x@annoStat)
}


#' Combine csAnno Object
#'
#' https://github.com/YuLab-SMU/ChIPseeker/issues/157
#' @title combine_csAnno
#' @param x csAnno object
#' @param ... csAnno objects
#' @return csAnno object
#' @examples
#' data(peakAnno)
#' combine_csAnno(peakAnno, peakAnno)
#' @importFrom methods new
#' @export
combine_csAnno <- function(x, ...) {
    z <- list(x, ...)

    if (sum(vapply(z, function(x) !is(x, "csAnno"), FUN.VALUE = logical(1))) != 0) {
        stop("not supported...")
    }

    if (length(z) < 2) {
        stop("need two or more csAnno object...")
    }


    if (sum(!duplicated(lapply(z, function(x) x@tssRegion[1]))) != 1 &&
        sum(!duplicated(lapply(z, function(x) x@tssRegion[2]))) != 1) {
        stop("the tss regions of different csAnno objects should be the same...")
    }

    if (sum(!duplicated(lapply(z, function(x) x@level))) != 1) {
        stop("the level of different csAnno object should be the same...")
    }

    if (sum(!duplicated(lapply(z, function(x) x@hasGenomicAnnotation))) != 1) {
        stop("the status of GenomicAnnotation should be the same...")
    }

    combine_tssRegion <- x@tssRegion
    combine_level <- x@level
    combine_hasGenomicAnnotation <- x@hasGenomicAnnotation

    combine_anno <- x@anno
    for (i in 2:length(z)) {
        combine_anno <- c(combine_anno, z[[i]]@anno)
    }

    combine_detailGenomicAnnotation <- lapply(z, function(x) x@detailGenomicAnnotation)
    combine_detailGenomicAnnotation <- do.call("rbind", combine_detailGenomicAnnotation)

    combine_peakNum <- x@peakNum
    for (i in 2:length(z)) {
        combine_peakNum <- combine_peakNum + z[[i]]@peakNum
    }

    feature <- x@annoStat$Feature
    for (i in 2:length(z)) {
        if (length(feature) < length(z[[i]]@annoStat$Feature)) {
            feature_levels <- levels(z[[i]]@annoStat$Feature)
            feature <- c(as.vector(feature), as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature,
                levels = feature_levels
            )
            feature <- sort(feature)
        } else {
            feature_levels <- levels(feature)
            feature <- c(as.vector(feature), as.vector(z[[i]]@annoStat$Feature))
            feature <- feature[!duplicated(feature)]
            feature <- factor(feature,
                levels = feature_levels
            )
            feature <- sort(feature)
        }
    }

    combine_annoStat <- data.frame(Feature = feature)

    for (i in seq_len(length(z))) {
        combine_annoStat <- merge(combine_annoStat, z[[i]]@annoStat,
            by = "Feature", all = TRUE, sort = FALSE
        )
        combine_annoStat[is.na(combine_annoStat)] <- 0
        combine_annoStat <- combine_annoStat[order(combine_annoStat$Feature), ]
    }

    total <- (ncol(combine_annoStat) - 1) * 100
    combine_annoStat$sum <- rowSums(combine_annoStat[, 2:ncol(combine_annoStat)])


    for (i in seq_len(length(combine_annoStat$sum))) {
        combine_annoStat$result[i] <- (combine_annoStat$sum[i] / total) * 100
    }

    annoStat_result <- data.frame(Feature = combine_annoStat[, 1], Frequency = combine_annoStat[, ncol(combine_annoStat)])

    res <- new("csAnno",
        anno = combine_anno,
        tssRegion = combine_tssRegion,
        level = combine_level,
        hasGenomicAnnotation = combine_hasGenomicAnnotation,
        detailGenomicAnnotation = combine_detailGenomicAnnotation,
        annoStat = annoStat_result,
        peakNum = combine_peakNum
    )

    return(res)
}

#' vennpie method generics
#'
#' @name vennpie
#' @docType methods
#' @rdname vennpie-methods
#'
#' @title vennpie method
#' @param x A `csAnno` instance
#' @param r initial radius
#' @param cex value to adjust legend
#' @param ... additional parameter
#' @return plot
#' @usage vennpie(x, r = 0.2, cex=1.2, ...)
#' @exportMethod vennpie
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "vennpie", signature(x = "csAnno"),
    function(x,
             r = 0.2,
             cex = 1.2,
             ...) {
        vennpie.csAnno(x, r, cex, ...)
    }
)


#' upsetplot method generics
#'
#' @name upsetplot
#' @docType methods
#' @rdname upsetplot-methods
#'
#' @title upsetplot method
#' @param x A `csAnno` instance
#' @param n number of categories to be plotted; only used by methods that
#' support it and ignored by the `csAnno` method
#' @param type plot type; only used by methods that support it and ignored by
#' the `csAnno` method
#' @param layer optional layer; only used by methods that support it and
#' ignored by the `csAnno` method
#' @param value score summary to display for overlapping features; only used
#' by methods that support it and ignored by the `csAnno` method
#' @param core_enrichment logical indicating whether only core features should
#' be used; only used by methods that support it and ignored by the `csAnno`
#' method
#' @param ... additional parameter
#' @return plot
#' @usage upsetplot(x, n = 10, type = "boxplot", layer = NULL,
#'     value = c("score", "abs_score"), core_enrichment = FALSE, ...)
#' @importFrom enrichplot upsetplot
#' @exportMethod upsetplot
#' @examples
#' data(peakAnno)
#' upsetplot(peakAnno)
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "upsetplot", signature(x = "csAnno"),
    function(x, ...) {
        upsetplot.csAnno(x, ...)
    }
)

#' convert csAnno object to data.frame
#'
#' @title as.data.frame.csAnno
#' @param x csAnno object
#' @param row.names row names
#' @param optional should be omitted.
#' @param ... additional parameters
#' @return data.frame
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
#' @method as.data.frame csAnno
#' @export
as.data.frame.csAnno <- function(x, row.names = NULL, optional = FALSE, ...) {
    y <- as.GRanges(x)
    if (!(is.null(row.names) || is.character(row.names))) {
        stop("'row.names' must be NULL or a character vector")
    }
    df <- as.data.frame(y)
    rownames(df) <- row.names
    return(df)
}

#' show method for `csAnno` instance
#'
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,csAnno,ANY-method
#' @title show method
#' @param object A `csAnno` instance
#' @return message
#' @importFrom methods show
#' @exportMethod show
#' @usage show(object)
#' @examples
#' peakfile <- system.file("extdata", "sample_peaks.txt", package = "epiSeeker")
#' show(peakfile)
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "show", signature(object = "csAnno"),
    function(object) {
        cat("Annotated peaks generated by epiSeeker\n")
        cat(
            paste(length(object@anno), object@peakNum, sep = "/"),
            " peaks were annotated\n"
        )
        src <- annoSource(object)
        if (!is.na(src$genome) || !is.null(src$TxDb)) {
            cat("Genome:", src$genome,
                if (!is.null(src$TxDb$source)) paste0("(TxDb: ", src$TxDb$source, ", ", src$TxDb$created, ")"),
                "\n")
        }
        chk <- annoChecks(object)
        if (isTRUE(chk$build$suspect)) {
            cat("Check: the peaks may not be on this genome build (see annoChecks())\n")
        } else if (isFALSE(chk$build$reliable)) {
            cat("Check: with", chk$build$n_peaks, "peaks the genome-build check is unreliable (see annoChecks())\n")
        }
        if (isFALSE(chk$assay$consistent)) {
            cat("Check: the peaks do not look like the declared assay,", chk$assay$assay, "(see annoChecks()$assay)\n")
        }
        if (object@hasGenomicAnnotation) {
            cat("Genomic Annotation Summary:\n")
            print(object@annoStat)
        }
    }
)

#' plotAnnoBar method for list of `csAnno` instances
#'
#' @name plotAnnoBar
#' @docType methods
#' @rdname plotAnnoBar-methods
#' @aliases plotAnnoBar,list-method
#' @exportMethod plotAnnoBar
setMethod(
    "plotAnnoBar", signature(x = "list"),
    function(x,
             xlab = "",
             ylab = "Percentage(%)",
             title = "Feature Distribution",
             ...) {
        if (is.null(names(x))) {
            nn <- paste0("Peak", seq_along(x))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
            names(x) <- nn
            ## stop("input object should be a named list...")
        }
        anno <- lapply(x, getAnnoStat)
        ## anno.df <- ldply(anno)
        anno.df <- list_to_dataframe(anno)
        categoryColumn <- ".id"
        plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
    }
)

#' plotAnnoBar method for `csAnno` instance
#'
#' @name plotAnnoBar
#' @docType methods
#' @rdname plotAnnoBar-methods
#' @aliases plotAnnoBar,csAnno,ANY-method
#' @title plotAnnoBar method
#' @param x `csAnno` instance
#' @param xlab xlab
#' @param ylab ylab
#' @param title title
#' @param ... additional paramter
#' @return plot
#' @exportMethod plotAnnoBar
#' @usage plotAnnoBar(x, xlab="", ylab='Percentage(\%)',title="Feature Distribution", ...)
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "plotAnnoBar", signature(x = "csAnno"),
    function(x,
             xlab = "",
             ylab = "Percentage(%)",
             title = "Feature Distribution",
             ...) {
        anno.df <- getAnnoStat(x)
        categoryColumn <- 1
        plotAnnoBar.data.frame(anno.df, xlab, ylab, title, categoryColumn)
    }
)


#' plotAnnoPie method for `csAnno` instance
#'
#' @name plotAnnoPie
#' @docType methods
#' @rdname plotAnnoPie-methods
#' @aliases plotAnnoPie,csAnno,ANY-method
#' @title plotAnnoPie method
#' @param x `csAnno` instance
#' @param ndigit number of digit to round
#' @param cex label cex
#' @param col color
#' @param legend.position topright or other.
#' @param pie3D plot in 3D or not
#' @param radius radius of the pie
#' @param ... extra parameter
#' @return plot
#' @exportMethod plotAnnoPie
#' @usage plotAnnoPie(x, ndigit = 2, cex = 0.9, col = NA,
#'     legend.position = "rightside", pie3D = FALSE,
#'     radius = 0.8, ...)
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "plotAnnoPie", signature(x = "csAnno"),
    function(x,
             ndigit = 2,
             cex = 0.9,
             col = NA,
             legend.position = "rightside",
             pie3D = FALSE,
             radius = 0.8,
             ...) {
        plotAnnoPie.csAnno(x, ndigit, cex, col, legend.position, pie3D, radius, ...)
    }
)


#' plotDistToTSS method for list of `csAnno` instances
#'
#' @name plotDistToTSS
#' @docType methods
#' @rdname plotDistToTSS-methods
#' @aliases plotDistToTSS,list-method
#' @exportMethod plotDistToTSS
setMethod(
    "plotDistToTSS", signature(x = "list"),
    function(x, distanceColumn = "distanceToTSS",
             xlab = "", ylab = "Binding sites (%) (5'->3')",
             title = "Distribution of transcription factor-binding loci relative to TSS",
             distanceBreaks = c(0, 1000, 3000, 5000, 10000, 100000),
             palette = NULL, ...) {
        if (is.null(names(x))) {
            nn <- paste0("Peak", seq_along(x))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
            names(x) <- nn
            ## stop("input object should be a named list...")
        }

        peakAnno <- lapply(x, as.data.frame)
        ## peakDist <- ldply(peakAnno)
        peakDist <- list_to_dataframe(peakAnno)
        categoryColumn <- ".id"
        plotDistToTSS.data.frame(peakDist,
            distanceColumn = distanceColumn,
            distanceBreaks = distanceBreaks, palette = palette,
            xlab = xlab, ylab = ylab, title = title, categoryColumn = categoryColumn
        )
    }
)


#' plotDistToTSS method for `csAnno` instance
#'
#' @name plotDistToTSS
#' @docType methods
#' @rdname plotDistToTSS-methods
#' @aliases plotDistToTSS,csAnno,ANY-method
#' @title plotDistToTSS method
#' @param distanceColumn distance column name
#' @param distanceBreaks breaks of distance, default is `c(0, 1000, 3000, 5000, 10000, 100000)`
#' @param palette palette name for coloring different distances. Run `RColorBrewer::display.brewer.all()` to see all applicable values.
#' @param x `csAnno` instance
#' @param xlab xlab
#' @param ylab ylab
#' @param title title
#' @param ... additional parameter
#' @return plot
#' @exportMethod plotDistToTSS
#' @usage plotDistToTSS(x,distanceColumn="distanceToTSS", xlab="",
#' ylab="Binding sites (\%) (5'->3')",
#' title="Distribution of transcription factor-binding loci relative to TSS",...)
#' @author Guangchuang Yu <https://guangchuangyu.github.io>
setMethod(
    "plotDistToTSS", signature(x = "csAnno"),
    function(x, distanceColumn = "distanceToTSS",
             xlab = "", ylab = "Binding sites (%) (5'->3')",
             title = "Distribution of transcription factor-binding loci relative to TSS",
             distanceBreaks = c(0, 1000, 3000, 5000, 10000, 100000),
             palette = NULL, ...) {
        peakDist <- as.data.frame(x)
        categoryColumn <- 1
        plotDistToTSS.data.frame(peakDist,
            distanceColumn = distanceColumn, distanceBreaks = distanceBreaks, palette = palette,
            xlab = xlab, ylab = ylab, title = title, categoryColumn = categoryColumn
        )
    }
)

#' Where the annotation of a csAnno came from
#'
#' The genome of the regions, and the source, organism and creation time of the `TxDb` and
#' the name and version of the `annoDb` that [annotateSeq()] used. Objects made before this
#' was recorded report only the genome stored with the regions.
#'
#' @param x a `csAnno` object
#' @return a list with `genome`, `TxDb` (`NULL` when not recorded) and `annoDb` (`NULL` when
#'   none was used or not recorded)
#' @examples
#' data(peakAnno)
#' annoSource(peakAnno)
#' @export
annoSource <- function(x) {
    if (!is(x, "csAnno")) stop("`x` must be a csAnno object.")
    src <- S4Vectors::metadata(x@anno)$annoSource
    if (is.null(src)) src <- list()
    if (is.null(src$genome) || is.na(src$genome)) {
        g <- unique(stats::na.omit(GenomeInfoDb::genome(x@anno)))
        src$genome <- if (length(g)) paste(g, collapse = ",") else NA_character_
    }
    src[c("genome", "TxDb", "annoDb")]
}

## Kept in metadata(<anno>) rather than in a slot of csAnno: a new slot would leave every
## csAnno saved by an earlier version without it, and Bioconductor's saveRDS() then fails on
## such an object ("no slot of name ..."). The metadata of a GRanges survives [, subset(),
## c() and sort(), so dropAnno(), subset() and combine_csAnno() keep it (the first object's).
.set_anno_source <- function(gr, TxDb, annoDb) {
    S4Vectors::metadata(gr)$annoSource <- .make_anno_source(TxDb, annoDb)
    gr
}

## what annotateSeq() records; reading metadata never stops an annotation
.make_anno_source <- function(TxDb, annoDb) {
    meta <- function(db) {
        m <- tryCatch(AnnotationDbi::metadata(db), error = function(e) NULL)
        if (is.null(m) || !all(c("name", "value") %in% names(m))) return(list())
        stats::setNames(as.list(m$value), m$name)
    }
    g <- tryCatch(unique(stats::na.omit(GenomeInfoDb::genome(TxDb))), error = function(e) character(0))
    txdb <- if (is(TxDb, "GRanges")) {
        list(source = "user-defined GRanges", organism = NA_character_, created = NA_character_)
    } else {
        m <- meta(TxDb)
        list(source = paste(stats::na.omit(c(m[["Data source"]], m[["UCSC Table"]], m[["UCSC Track"]])), collapse = " "),
             organism = if (is.null(m$Organism)) NA_character_ else m$Organism,
             created = if (is.null(m[["Creation time"]])) NA_character_ else sub(" .*", "", m[["Creation time"]]))
    }
    db <- NULL
    if (!is.null(annoDb)) {
        name <- if (is.character(annoDb)) annoDb else class(annoDb)[1]
        version <- if (is.character(annoDb)) tryCatch(as.character(utils::packageVersion(annoDb)), error = function(e) NA_character_) else NA_character_
        m <- meta(if (is.character(annoDb)) tryCatch(getExportedValue(annoDb, annoDb), error = function(e) NULL) else annoDb)
        db <- list(name = name, version = version,
                   gene_source_date = if (is.null(m$EGSOURCEDATE)) NA_character_ else m$EGSOURCEDATE)
    }
    list(genome = if (length(g)) paste(g, collapse = ",") else NA_character_, TxDb = txdb, annoDb = db)
}

#' Class "csAnnoList"
#'
#' A named list of `csAnno` objects, one per sample or group, for comparing annotations
#' (the counterpart of a `compareClusterResult`). It is a list, so every function that
#' accepts a list of `csAnno` objects (`plotAnnoBar()`, `plotDistToTSS()`, `lapply()`) accepts
#' it as well; the class adds the check that all elements are named `csAnno` objects, and
#' methods such as `interpret()` can recognise it.
#'
#' @name csAnnoList-class
#' @docType class
#' @exportClass csAnnoList
#' @seealso [csAnnoList()]
#' @keywords classes
setClass("csAnnoList", contains = "namedList", validity = function(object) {
    msg <- character(0)
    if (!length(object)) msg <- c(msg, "a csAnnoList needs at least one csAnno object")
    if (!all(vapply(object, is, logical(1), "csAnno"))) msg <- c(msg, "every element must be a csAnno object")
    nm <- names(object)
    if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) msg <- c(msg, "every element needs a unique name")
    if (length(msg)) msg else TRUE
})

#' Make a csAnnoList
#'
#' @param ... named `csAnno` objects, or one named list of them (for example the result of
#'   `lapply(files, annotateSeq, ...)`)
#' @return a `csAnnoList`
#' @examples
#' data(peakAnno)
#' x <- csAnnoList(a = peakAnno, b = peakAnno)
#' plotAnnoBar(x)
#' @export
csAnnoList <- function(...) {
    x <- list(...)
    if (length(x) == 1 && is.list(x[[1]]) && !is(x[[1]], "csAnno")) x <- x[[1]]
    if (is(x, "csAnnoList")) return(x)
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
        stop("every element needs a unique name, e.g. csAnnoList(control = x, treated = y)")
    }
    new("csAnnoList", .Data = unname(x), names = names(x))
}
