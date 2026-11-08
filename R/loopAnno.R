## Chromatin loops: the pair group. A loop is a pair of anchors; each anchor is annotated by
## annotateSeq() (one annotation engine for regions and anchors), and genes are linked to
## a peak through the anchor at the other end of a loop. Loops are represented as
## InteractionSet's GInteractions, the class other Bioconductor packages (HiCExperiment,
## GenomicInteractions, plyinteractions) use, so their loops need no conversion.

#' Class "loopAnno"
#'
#' Annotated chromatin loops, made by [annotateLoops()].
#'
#' @slot loops a `GInteractions` with one row per loop: the input columns, and the feature,
#'   promoter status and genes of both anchors, the loop class (promoter-promoter,
#'   promoter-distal, distal-distal) and the span (distance between anchor centres, `NA`
#'   for loops between chromosomes). Where the loops came from is in
#'   `metadata(x@loops)$loopSource`.
#' @slot anchors a `csAnno`: `annotateSeq()` on the anchors, `regions(loops)`, in the same
#'   order (column `region_id`); anchors that `annotateSeq()` dropped are missing from it.
#' @slot tssRegion the TSS region used to call an anchor a promoter anchor
#' @slot loopStat number and percentage of loops per class
#' @name loopAnno-class
#' @docType class
#' @exportClass loopAnno
#' @importClassesFrom InteractionSet GInteractions
#' @seealso [annotateLoops()], [linkPeaksToGenes()], [plotLoopAnno()]
#' @keywords classes
setClass("loopAnno", representation(
    loops = "GInteractions",
    anchors = "csAnno",
    tssRegion = "numeric",
    loopStat = "data.frame"
))

.loop_classes <- c("promoter-promoter", "promoter-distal", "distal-distal")
.link_classes <- c("promoter", "no loop", "loop, no gene", "loop includes nearest", "loop, other genes")

#' Read chromatin loops
#'
#' Reads loops from a BEDPE file (also gzipped; a header line starting with `#`, as written
#' by Juicer/HiCCUPS, or a plain header line is used for the column names), from a data
#' frame with the BEDPE columns, or takes them from a `GInteractions` or `Pairs` object.
#' BEDPE starts are 0-based and become 1-based, as in [readPeakFile()]. Columns after the
#' six coordinates are kept as `mcols()`.
#'
#' @param x a BEDPE file, a data frame whose first six columns are BEDPE coordinates, a
#'   `GInteractions` or a `Pairs` of `GRanges`
#' @return a `GInteractions`; the source, number of rows read and number of rows dropped
#'   (missing coordinates) are in `metadata(<result>)$loopSource`
#' @examples
#' f <- tempfile(fileext = ".bedpe")
#' writeLines(c("chr1\t1000\t2000\tchr1\t50000\t51000",
#'              "chr1\t1000\t2000\tchr1\t90000\t91000"), f)
#' readLoops(f)
#' @importFrom InteractionSet GInteractions makeGInteractionsFromGRangesPairs
#' @importFrom S4Vectors metadata metadata<- mcols mcols<- DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
readLoops <- function(x) {
    if (is(x, "loopAnno")) return(x@loops)
    if (is(x, "GInteractions")) {
        if (is.null(metadata(x)$loopSource)) {
            metadata(x)$loopSource <- list(source = class(x)[1], rows = length(x),
                                           loops = length(x), dropped = 0L)
        }
        return(x)
    }
    if (is(x, "Pairs")) {
        gi <- makeGInteractionsFromGRangesPairs(x)
        metadata(gi)$loopSource <- list(source = "Pairs", rows = length(x), loops = length(gi),
                                        dropped = 0L)
        return(gi)
    }
    if (is.character(x) && length(x) == 1) {
        if (!file.exists(x)) stop("file '", x, "' does not exist.")
        df <- .read_bedpe(x)
        src <- basename(x)
    } else if (is.data.frame(x)) {
        df <- x
        src <- "data.frame"
    } else {
        stop("`x` must be a BEDPE file, a data frame, a GInteractions or a Pairs object.")
    }
    .bedpe_to_gi(df, src)
}

.read_bedpe <- function(file) {
    head_lines <- readLines(file, n = 50, warn = FALSE)
    hdr <- NULL
    h <- sub("^#\\s*", "", head_lines[grepl("^#", head_lines)])
    h <- h[lengths(strsplit(h, "\t")) >= 6]
    if (length(h)) hdr <- strsplit(h[1], "\t")[[1]]
    df <- utils::read.delim(file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE,
                            check.names = FALSE, colClasses = "character")
    ## a plain header line: its start column is not a number
    if (nrow(df) && !grepl("^\\s*[0-9]+\\s*$", df[1, 2])) {
        hdr <- unlist(df[1, ], use.names = FALSE)
        df <- df[-1, , drop = FALSE]
    }
    df[] <- lapply(df, utils::type.convert, as.is = TRUE)
    if (!is.null(hdr) && length(hdr) == ncol(df)) names(df) <- hdr
    df
}

.bedpe_to_gi <- function(df, src) {
    if (ncol(df) < 6) stop("a BEDPE table needs at least six columns (chrom1, start1, end1, chrom2, start2, end2).")
    n <- nrow(df)
    co <- df[, 1:6]
    ok <- stats::complete.cases(co) & !(as.character(co[[1]]) %in% c("", ".")) &
        !(as.character(co[[4]]) %in% c("", "."))
    df <- df[ok, , drop = FALSE]
    a1 <- GRanges(as.character(df[[1]]), IRanges(as.numeric(df[[2]]) + 1, as.numeric(df[[3]])))
    a2 <- GRanges(as.character(df[[4]]), IRanges(as.numeric(df[[5]]) + 1, as.numeric(df[[6]])))
    gi <- GInteractions(a1, a2)
    if (ncol(df) > 6) mcols(gi) <- DataFrame(df[, -(1:6), drop = FALSE], check.names = FALSE)
    metadata(gi)$loopSource <- list(source = src, rows = n, loops = length(gi), dropped = sum(!ok))
    gi
}

#' Annotate chromatin loops
#'
#' Annotates both anchors of every loop with [annotateSeq()], marks an anchor as a promoter
#' anchor when it overlaps the promoter window (`tssRegion`) of a transcript (or, with
#' `level = "gene"`, a gene) of `TxDb` - the window `annotateSeq()` calls "Promoter" - and
#' records the genes of those promoters, since an anchor of several kilobases can hold more
#' than one. Each loop is then classed as promoter-promoter, promoter-distal or
#' distal-distal.
#'
#' With `gene_ids` only those genes count (for example protein-coding genes, as in the loop
#' figure of the epiSeeker paper): an anchor is a promoter anchor only when it holds the
#' promoter of one of them, and only they are recorded as anchor genes. The ids given are kept
#' in `metadata(<result>@loops)$geneFilter`.
#'
#' @param loops anything [readLoops()] reads
#' @param TxDb `TxDb`, `EnsDb` or a `GRanges` of features with a `gene_id` column
#' @param tssRegion TSS region, as in [annotateSeq()]
#' @param level `"transcript"` or `"gene"`, as in [annotateSeq()]
#' @param annoDb annotation package or `OrgDb`; adds gene symbols to the anchors and the
#'   anchor genes
#' @param gene_ids optional gene ids (as in `TxDb`) to restrict the anchor genes to, e.g.
#'   protein-coding genes; `NULL` (default) keeps every gene of `TxDb`
#' @param verbose logical
#' @param ... further arguments to [annotateSeq()]
#' @return a [loopAnno-class] object
#' @examples
#' \donttest{
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' f <- tempfile(fileext = ".bedpe")
#' writeLines(c("chr17\t7685000\t7690000\tchr17\t7740000\t7745000",
#'              "chr17\t7500000\t7505000\tchr17\t7600000\t7605000"), f)
#' la <- annotateLoops(f, TxDb = txdb, verbose = FALSE)
#' la
#' }
#' @importFrom InteractionSet regions anchorIds pairdist
#' @importFrom IRanges findOverlaps CharacterList
#' @importFrom S4Vectors queryHits subjectHits splitAsList
#' @importFrom IRanges promoters
#' @export
annotateLoops <- function(loops, TxDb = NULL, tssRegion = c(-3000, 3000), level = "transcript",
                          annoDb = NULL, gene_ids = NULL, verbose = TRUE, ...) {
    gi <- readLoops(loops)
    if (!length(gi)) stop("there are no loops to annotate.")
    if (!is(TxDb, "GRanges")) TxDb <- loadTxDb(TxDb)
    level <- match.arg(level, c("transcript", "gene"))

    reg <- regions(gi)
    mcols(reg) <- NULL
    reg$region_id <- seq_along(reg)
    anchors <- annotateSeq(reg, tssRegion = tssRegion, TxDb = TxDb, level = level,
                           annoDb = annoDb, verbose = verbose, ...)
    ag <- anchors@anno
    feature <- rep(NA_character_, length(reg))
    if (!is.null(ag$annotation)) feature[ag$region_id] <- .loop_feature(ag$annotation)

    genes <- .anchor_promoter_genes(reg, TxDb, tssRegion, level)
    if (!is.null(gene_ids)) {
        keep <- unique(as.character(gene_ids))
        ## element-wise: base `%in%` on a List matches whole elements, not the ids inside them
        genes <- genes[relist(unlist(genes, use.names = FALSE) %in% keep, genes)]
        metadata(gi)$geneFilter <- list(n_gene_ids = length(keep))
    }
    promoter <- lengths(genes) > 0
    i1 <- anchorIds(gi, type = "first")
    i2 <- anchorIds(gi, type = "second")
    cls <- ifelse(promoter[i1] & promoter[i2], .loop_classes[1],
                  ifelse(promoter[i1] | promoter[i2], .loop_classes[2], .loop_classes[3]))

    mc <- mcols(gi)
    mc$anchor1_feature <- feature[i1]
    mc$anchor2_feature <- feature[i2]
    mc$anchor1_promoter <- promoter[i1]
    mc$anchor2_promoter <- promoter[i2]
    mc$anchor1_genes <- genes[i1]
    mc$anchor2_genes <- genes[i2]
    if (!is.null(annoDb)) {
        sym <- .gene_symbols(unique(unlist(genes)), TxDb, annoDb)
        to_sym <- function(g) CharacterList(lapply(g, function(z) unique(stats::na.omit(sym[z]))))
        mc$anchor1_symbols <- to_sym(genes[i1])
        mc$anchor2_symbols <- to_sym(genes[i2])
    }
    mc$loop_class <- factor(cls, levels = .loop_classes)
    mc$span <- pairdist(gi, type = "mid")
    mcols(gi) <- mc

    tab <- table(mc$loop_class)
    stat <- data.frame(loop_class = factor(names(tab), levels = .loop_classes),
                       n = as.integer(tab), Frequency = as.numeric(100 * tab / length(gi)))
    new("loopAnno", loops = gi, anchors = anchors, tssRegion = tssRegion, loopStat = stat)
}

## "Intron (uc001clv.1/55182, intron 2 of 14)" -> "Intron"; promoters keep their window
.loop_feature <- function(a) {
    a <- as.character(a)
    keep <- grepl("^Promoter", a)
    a[!keep] <- sub(" \\(.*$", "", a[!keep])
    a
}

## genes whose promoter window (tssRegion around each TSS) overlaps each anchor
.anchor_promoter_genes <- function(reg, TxDb, tssRegion, level) {
    if (is(TxDb, "GRanges")) {
        features <- TxDb
        gid <- as.character(features$gene_id)
    } else if (level == "gene") {
        features <- getGene(TxDb, by = "gene")
        gid <- as.character(features$gene_id)
    } else {
        ## transcriptsBy(TxDb, by = "gene"), unlisted: names are the gene ids
        features <- getGene(TxDb, by = "transcript")
        gid <- names(features)
    }
    prom <- suppressWarnings(promoters(features, upstream = -tssRegion[1], downstream = tssRegion[2]))
    hits <- suppressWarnings(findOverlaps(reg, prom, ignore.strand = TRUE))
    g <- gid[subjectHits(hits)]
    keep <- !is.na(g)
    out <- splitAsList(g[keep], factor(queryHits(hits)[keep], levels = seq_along(reg)))
    names(out) <- NULL
    unique(out)
}

## gene id -> SYMBOL with the same OrgDb lookup annotateSeq() uses
.gene_symbols <- function(ids, TxDb, annoDb) {
    ids <- ids[!is.na(ids)]
    if (!length(ids)) return(stats::setNames(character(0), character(0)))
    type <- if (is(TxDb, "GRanges")) character(0) else IDType(TxDb)
    if (length(type) == 0 || is.na(type)) {
        type <- if (all(grepl("^ENS", ids))) "Ensembl Gene ID" else "Entrez Gene ID"
    }
    ga <- tryCatch(getGeneAnno(annoDb, ids, type, "SYMBOL"), error = function(e) NULL)
    if (!is.data.frame(ga) || is.null(ga$SYMBOL)) return(stats::setNames(rep(NA_character_, length(ids)), ids))
    stats::setNames(as.character(ga$SYMBOL), ids)
}

#' Link peaks to genes through chromatin loops
#'
#' For each peak, finds the loops with an anchor that overlaps the peak and the genes at
#' the promoter anchors on the other end, and compares them with the peak's nearest gene
#' from [annotateSeq()]. Peaks in promoters are reported but not compared.
#'
#' The link of a peak is one of `"promoter"` (the peak is in a promoter), `"no loop"` (no
#' anchor overlaps the peak), `"loop, no gene"` (the other ends are not promoter anchors),
#' `"loop includes nearest"` (the nearest gene is among the loop genes) and
#' `"loop, other genes"` (the loops point to genes other than the nearest one).
#'
#' The defaults use the annotation of [annotateSeq()] and [annotateLoops()]. Stricter
#' definitions, as used to score loops against a CRISPRi enhancer-gene benchmark, are
#' options: restrict the genes (`gene_ids`, e.g. protein-coding genes), call the genes of an
#' anchor those with a transcript start inside the anchor extended by `anchor_flank` bp
#' (instead of the promoter window of `annotateLoops()`), take the nearest gene from the
#' transcript starts of those genes and call a peak a promoter peak within `distal_distance`
#' bp of one (instead of the annotation), and leave out the genes of the peak's own anchors
#' (`exclude_own_anchor`). `TxDb` is needed for `anchor_flank` and `distal_distance`. With
#' these options, symbols are reported only when `annoDb` is given.
#'
#' @param peak a `csAnno` from [annotateSeq()], made with the same `TxDb` as `loops`
#' @param loops a [loopAnno-class] object from [annotateLoops()]
#' @param TxDb `TxDb`, `EnsDb` or `GRanges` of features with `gene_id`; needed for
#'   `anchor_flank` and `distal_distance`
#' @param gene_ids optional gene ids to consider (as in `TxDb`), e.g. protein-coding genes:
#'   `AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db), "GENETYPE")`, those with
#'   `GENETYPE == "protein-coding"`
#' @param anchor_flank optional; genes of an anchor = genes with a transcript start within the
#'   anchor extended by this many bp on each side
#' @param distal_distance optional; nearest gene = gene of the nearest transcript start to the
#'   peak centre, and a peak within this many bp of it is a promoter peak
#' @param exclude_own_anchor logical; leave out of a peak's loop genes the genes of the
#'   anchors the peak is in. Default `FALSE`.
#' @param annoDb optional annotation package or `OrgDb`; with `anchor_flank` or
#'   `distal_distance`, adds `loop_symbols` and `nearest_symbol` for the gene ids
#' @return a data frame with one row per peak, in the order of `peak`: coordinates,
#'   `feature`, `nearest_gene` (and `nearest_symbol` when available), `n_loops`,
#'   `loop_genes` (and `loop_symbols`), separated by "/", and `link`
#' @examples
#' \donttest{
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' f <- tempfile(fileext = ".bedpe")
#' writeLines("chr17\t7685000\t7690000\tchr17\t7740000\t7745000", f)
#' la <- annotateLoops(f, TxDb = txdb, verbose = FALSE)
#' pk <- GenomicRanges::GRanges("chr17", IRanges::IRanges(7741000, 7741500))
#' linkPeaksToGenes(annotateSeq(pk, TxDb = txdb, verbose = FALSE), la)
#' }
#' @export
linkPeaksToGenes <- function(peak, loops, TxDb = NULL, gene_ids = NULL, anchor_flank = NULL,
                             distal_distance = NULL, exclude_own_anchor = FALSE, annoDb = NULL) {
    if (!is(peak, "csAnno")) {
        stop("`peak` must be a csAnno: annotate it with annotateSeq() and the TxDb used for the loops.")
    }
    if (!is(loops, "loopAnno")) stop("`loops` must be a loopAnno from annotateLoops().")
    pg <- peak@anno
    gi <- loops@loops
    mc <- mcols(gi)
    i1 <- anchorIds(gi, type = "first")
    i2 <- anchorIds(gi, type = "second")
    reg <- regions(gi)
    strict <- !is.null(anchor_flank) || !is.null(distal_distance)
    if (strict && is.null(TxDb)) stop("`anchor_flank` and `distal_distance` need `TxDb`.")
    if (strict) tss <- .tss_of(TxDb, gene_ids)

    ## genes of every anchor: stored by annotateLoops(), or from transcript starts
    anchor_genes <- if (!is.null(anchor_flank)) {
        .genes_near(reg, tss, anchor_flank)
    } else {
        ag <- vector("list", length(reg))
        ag[i1] <- as.list(mc$anchor1_genes); ag[i2] <- as.list(mc$anchor2_genes)
        lapply(ag, function(z) unique(as.character(unlist(z))))
    }
    if (!is.null(gene_ids)) anchor_genes <- lapply(anchor_genes, intersect, as.character(gene_ids))

    ## for every anchor: the genes at the other end of its loops, and how many loops
    a <- c(i1, i2)
    other <- c(anchor_genes[i2], anchor_genes[i1])
    by_anchor <- split(other, factor(a, levels = seq_along(reg)))
    partner <- lapply(by_anchor, function(z) unique(unlist(z)))
    n_by_anchor <- tabulate(a, nbins = length(reg))
    have_sym <- !is.null(mc$anchor1_symbols) && is.null(anchor_flank) && is.null(gene_ids)
    if (have_sym) {
        other_sym <- c(as.list(mc$anchor2_symbols), as.list(mc$anchor1_symbols))
        partner_sym <- lapply(split(other_sym, factor(a, levels = seq_along(reg))),
                              function(z) unique(unlist(z)))
    }

    hits <- suppressWarnings(findOverlaps(pg, reg, ignore.strand = TRUE))
    per_peak <- split(subjectHits(hits), factor(queryHits(hits), levels = seq_along(pg)))
    loop_genes <- lapply(per_peak, function(ids) unique(unlist(partner[ids])))
    if (exclude_own_anchor) {
        ## only anchors that are on a loop: a subset GInteractions keeps unused regions
        own <- lapply(per_peak, function(ids) unique(unlist(anchor_genes[ids[n_by_anchor[ids] > 0]])))
        loop_genes <- mapply(setdiff, loop_genes, own, SIMPLIFY = FALSE)
    }
    n_loops <- vapply(per_peak, function(ids) sum(n_by_anchor[ids]), numeric(1))

    feature <- if (!is.null(pg$annotation)) .loop_feature(pg$annotation) else rep(NA_character_, length(pg))
    in_promoter <- if (!is.null(pg$annotation)) {
        grepl("^Promoter", pg$annotation)
    } else {
        d <- pg$distanceToTSS
        !is.na(d) & d >= peak@tssRegion[1] & d <= peak@tssRegion[2]
    }
    ## a self-defined GRanges TxDb keeps its own gene_id column instead of geneId
    nearest <- if (!is.null(pg$geneId)) as.character(pg$geneId) else
        if (!is.null(pg$gene_id)) as.character(pg$gene_id) else rep(NA_character_, length(pg))
    if (!is.null(distal_distance)) {
        nn <- .nearest_tss(pg, tss)
        nearest <- nn$gene
        in_promoter <- !is.na(nn$distance) & nn$distance <= distal_distance
    }
    link <- ifelse(in_promoter, "promoter",
            ifelse(n_loops == 0, "no loop",
            ifelse(lengths(loop_genes) == 0, "loop, no gene",
            ifelse(mapply(function(n, g) !is.na(n) && n %in% g, nearest, loop_genes),
                   "loop includes nearest", "loop, other genes"))))

    out <- data.frame(seqnames = as.character(GenomicRanges::seqnames(pg)),
                      start = GenomicRanges::start(pg), end = GenomicRanges::end(pg),
                      feature = feature, nearest_gene = nearest, stringsAsFactors = FALSE)
    if (!is.null(pg$SYMBOL) && is.null(distal_distance)) out$nearest_symbol <- as.character(pg$SYMBOL)
    out$n_loops <- as.integer(n_loops)
    out$loop_genes <- vapply(loop_genes, paste, character(1), collapse = "/")
    if (have_sym) {
        out$loop_symbols <- vapply(per_peak, function(ids) paste(unique(unlist(partner_sym[ids])), collapse = "/"),
                                   character(1))
    }
    if ((!is.null(anchor_flank) || !is.null(gene_ids) || !is.null(distal_distance)) && !is.null(annoDb)) {
        sym <- .gene_symbols(unique(c(unlist(loop_genes), nearest)), TxDb, annoDb)
        to_sym <- function(g) paste(unique(stats::na.omit(sym[g])), collapse = "/")
        out$loop_symbols <- vapply(loop_genes, to_sym, character(1))
        if (!is.null(distal_distance)) out$nearest_symbol <- unname(sym[nearest])
    }
    out$link <- factor(link, levels = .link_classes)
    rownames(out) <- NULL
    out
}

## transcript starts (one per gene and position) of TxDb, EnsDb or a GRanges of features,
## optionally only of `gene_ids`
.tss_of <- function(TxDb, gene_ids = NULL) {
    if (is(TxDb, "GRanges")) {
        if (is.null(TxDb$gene_id)) stop("a GRanges `TxDb` needs a gene_id column.")
        tss <- GenomicRanges::resize(GenomicRanges::granges(TxDb), 1, fix = "start")
        tss$gene_id <- as.character(TxDb$gene_id)
    } else {
        tx <- GenomicFeatures::transcripts(TxDb, columns = "gene_id")
        gid <- tx$gene_id
        tss <- GenomicRanges::resize(GenomicRanges::granges(rep(tx, lengths(gid))), 1, fix = "start")
        tss$gene_id <- as.character(unlist(gid))
    }
    tss <- tss[!is.na(tss$gene_id)]
    if (!is.null(gene_ids)) tss <- tss[tss$gene_id %in% as.character(gene_ids)]
    unique(tss)
}

## genes with a transcript start within each region extended by `flank` on both sides
.genes_near <- function(reg, tss, flank) {
    hits <- suppressWarnings(findOverlaps(reg + flank, tss, ignore.strand = TRUE))
    g <- splitAsList(tss$gene_id[subjectHits(hits)], factor(queryHits(hits), levels = seq_along(reg)))
    lapply(as.list(g), unique)
}

## gene of the nearest transcript start to each region centre, and the distance
.nearest_tss <- function(x, tss) {
    ctr <- GenomicRanges::resize(GenomicRanges::granges(x), 1, fix = "center")
    d <- suppressWarnings(GenomicRanges::distanceToNearest(ctr, tss, ignore.strand = TRUE))
    gene <- rep(NA_character_, length(x)); dist <- rep(NA_real_, length(x))
    gene[queryHits(d)] <- tss$gene_id[subjectHits(d)]
    dist[queryHits(d)] <- S4Vectors::mcols(d)$distance
    list(gene = gene, distance = dist)
}

#' Plot annotated loops
#'
#' @param x a [loopAnno-class] object
#' @param type `"class"`: number of loops per class; `"span"`: distribution of the loop span
#'   per class (loops within one chromosome)
#' @return a ggplot object
#' @examples
#' \donttest{
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' f <- tempfile(fileext = ".bedpe")
#' writeLines(c("chr17\t7685000\t7690000\tchr17\t7740000\t7745000",
#'              "chr17\t7500000\t7505000\tchr17\t7600000\t7605000"), f)
#' plotLoopAnno(annotateLoops(f, TxDb = txdb, verbose = FALSE))
#' }
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_boxplot labs theme_bw scale_y_log10
#' @export
plotLoopAnno <- function(x, type = c("class", "span")) {
    if (!is(x, "loopAnno")) stop("`x` must be a loopAnno from annotateLoops().")
    type <- match.arg(type)
    if (type == "class") {
        d <- x@loopStat
        d$label <- paste0(d$n, " (", round(d$Frequency, 1), "%)")
        return(ggplot(d, aes(x = .data$loop_class, y = .data$n)) +
                   geom_col(fill = "steelblue") +
                   geom_text(aes(label = .data$label), vjust = -0.3, size = 3.5) +
                   labs(x = NULL, y = "Number of loops") + theme_bw())
    }
    mc <- mcols(x@loops)
    d <- data.frame(loop_class = mc$loop_class, span = mc$span)
    d <- d[!is.na(d$span) & d$span > 0, , drop = FALSE]
    ggplot(d, aes(x = .data$loop_class, y = .data$span)) +
        geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5) +
        scale_y_log10() + labs(x = NULL, y = "Loop span (bp)") + theme_bw()
}

#' Plot how peaks are linked to genes through loops
#'
#' @param links the data frame returned by [linkPeaksToGenes()]
#' @param distal_only logical; plot only the peaks outside promoters (default `TRUE`)
#' @return a ggplot object
#' @examples
#' links <- data.frame(link = factor(c("no loop", "loop, other genes", "promoter"),
#'                                   levels = c("promoter", "no loop", "loop, no gene",
#'                                              "loop includes nearest", "loop, other genes")))
#' plotLoopLinks(links)
#' @export
plotLoopLinks <- function(links, distal_only = TRUE) {
    if (!is.data.frame(links) || is.null(links$link)) stop("`links` must come from linkPeaksToGenes().")
    l <- factor(links$link, levels = .link_classes)
    if (distal_only) l <- factor(l[l != "promoter"], levels = setdiff(.link_classes, "promoter"))
    tab <- table(l)
    d <- data.frame(link = factor(names(tab), levels = names(tab)), n = as.integer(tab))
    d$label <- paste0(d$n, " (", round(100 * d$n / max(1, sum(d$n)), 1), "%)")
    ggplot(d, aes(x = .data$link, y = .data$n)) +
        geom_col(fill = "darkorange") +
        geom_text(aes(label = .data$label), vjust = -0.3, size = 3.5) +
        labs(x = NULL, y = if (distal_only) "Peaks outside promoters" else "Peaks") + theme_bw()
}

#' @rdname loopAnno-class
#' @param object a `loopAnno` object
#' @exportMethod show
setMethod("show", signature(object = "loopAnno"), function(object) {
    src <- metadata(object@loops)$loopSource
    cat("Annotated loops generated by epiSeeker\n")
    cat(length(object@loops), "loops between", length(regions(object@loops)), "anchors",
        if (!is.null(src)) paste0("(from ", src$source, ")"), "\n")
    g <- annoSource(object@anchors)$genome
    if (!is.na(g)) cat("Genome:", g, "\n")
    cat("Loop classes:\n")
    print(object@loopStat, row.names = FALSE)
})

#' @method as.data.frame loopAnno
#' @export
as.data.frame.loopAnno <- function(x, row.names = NULL, optional = FALSE, ...) {
    gi <- x@loops
    mc <- mcols(gi)
    for (cn in names(mc)) {
        if (is(mc[[cn]], "CharacterList")) mc[[cn]] <- vapply(as.list(mc[[cn]]), paste, character(1), collapse = "/")
    }
    a1 <- InteractionSet::anchors(gi, type = "first")
    a2 <- InteractionSet::anchors(gi, type = "second")
    data.frame(seqnames1 = as.character(GenomicRanges::seqnames(a1)), start1 = GenomicRanges::start(a1),
               end1 = GenomicRanges::end(a1), seqnames2 = as.character(GenomicRanges::seqnames(a2)),
               start2 = GenomicRanges::start(a2), end2 = GenomicRanges::end(a2),
               as.data.frame(mc, optional = TRUE), stringsAsFactors = FALSE, check.names = FALSE)
}
