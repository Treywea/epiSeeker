## Adapters from single-cell and spatial chromatin data to the range result classes.
##
## The upstream tools (Signac, ArchR, SnapATAC2, a spatial pipeline) do the clustering; these
## functions take their output -- a fragments file with a barcode -> group table, or the peaks
## of each group -- and turn each group into a pseudo-bulk region set, annotated as a csAnno
## and collected in a named csAnnoList. What the groups are (resolution, number of cells or
## spots, spot coordinates) is kept in metadata(<anno>)$group, like annoSource, so the classes
## are unchanged and annotation, comparison, overlap and interpret() work as for bulk data.

#' Read a fragments file
#'
#' Reads the fragments file of a single-cell or spatial chromatin assay (10x Genomics, ArchR,
#' Signac and spatial ATAC / CUT&Tag pipelines): tab-separated `chr start end barcode count`,
#' 0-based start, optionally gzipped, `#` lines ignored. The whole file is read into memory.
#'
#' @param file path of the fragments file
#' @param which optional `GRanges`; only fragments overlapping it are kept
#' @param barcodes optional character vector; only fragments of these barcodes are kept
#' @return a `GRanges` with metadata columns `barcode` and `count`
#' @seealso [fragmentsToGroups()]
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c("chr1\t100\t300\tAAAC-1\t1", "chr1\t500\t700\tAAAG-1\t2"), f)
#' readFragments(f)
#' @export
readFragments <- function(file, which = NULL, barcodes = NULL) {
    df <- utils::read.delim(file, header = FALSE, comment.char = "#", quote = "",
                            stringsAsFactors = FALSE)
    if (ncol(df) < 4) stop("a fragments file needs at least 4 columns: chr, start, end, barcode")
    count <- if (ncol(df) >= 5) as.integer(df[[5]]) else rep(1L, nrow(df))
    gr <- GenomicRanges::GRanges(df[[1]], IRanges::IRanges(df[[2]] + 1L, df[[3]]),
                                 barcode = as.character(df[[4]]), count = count)
    if (!is.null(barcodes)) gr <- gr[gr$barcode %in% barcodes]
    if (!is.null(which)) gr <- IRanges::subsetByOverlaps(gr, which)
    gr
}

## barcode -> group as a named character vector, from a named vector or a data.frame
.as_group_map <- function(groups) {
    if (is.data.frame(groups)) {
        bc <- intersect(c("barcode", "cell", "spot"), names(groups))[1]
        gc <- intersect(c("group", "cluster", "domain", "cell_type", "celltype"), names(groups))[1]
        if (is.na(bc) || is.na(gc)) {
            stop("`groups` as a data.frame needs a barcode column (barcode, cell or spot) and a ",
                 "group column (group, cluster, domain or cell_type)")
        }
        return(stats::setNames(as.character(groups[[gc]]), as.character(groups[[bc]])))
    }
    if (is.null(names(groups))) stop("`groups` must be named by barcode, or a data.frame")
    stats::setNames(as.character(groups), names(groups))
}

## match fragment barcodes to the group table; when fewer than half match exactly, compare
## them without a "-1"-like suffix, which fragments files carry and many tables do not
.match_barcodes <- function(frag_bc, map) {
    strip <- function(b) sub("-[0-9]+$", "", b)
    idx <- match(frag_bc, names(map))
    used <- "exact"
    ub <- unique(frag_bc)
    if (mean(names(map) %in% ub) < 0.5) {
        idx2 <- match(strip(frag_bc), strip(names(map)))
        if (sum(!is.na(idx2)) > sum(!is.na(idx))) { idx <- idx2; used <- "without the -N suffix" }
    }
    found <- unique(names(map)[idx[!is.na(idx)]])
    list(group = unname(map[idx]), matched = length(found), total = length(map), how = used)
}

#' Pseudo-bulk region sets from fragments and cell or spot groups
#'
#' Sums the fragments of the barcodes of each group in genomic tiles and keeps, for each group,
#' the tiles with enough fragments -- and, with `specific = TRUE`, more fragments per million
#' than the other groups. This turns clusters of cells, or spatial domains of spots, into one
#' region set per group, ready for [annotateGroups()]. It does not call peaks; for peaks called
#' per cluster use [peaksToGroups()].
#'
#' Fragment barcodes are matched to `groups` exactly; when fewer than half of the barcodes of
#' `groups` match, they are matched without a `-1`-like suffix, which fragments files often
#' carry and spot or cell tables often do not. How many matched is reported.
#'
#' @param fragments a `GRanges` from [readFragments()] or the path of a fragments file
#' @param groups barcode -> group: a character vector named by barcode, or a data.frame with a
#'   barcode column (`barcode`, `cell` or `spot`) and a group column (`group`, `cluster`,
#'   `domain` or `cell_type`)
#' @param tile tile width in bp; e.g. 500 for ATAC, 5000 for a broad histone mark
#' @param min_fragments minimum fragments of a group in a tile
#' @param specific keep only tiles enriched in the group over the other groups
#' @param fold minimum enrichment, `(CPM in group + 1) / (CPM in the other groups + 1)`, when
#'   `specific = TRUE`
#' @param max_regions keep at most this many tiles per group, highest CPM first
#' @return a named `GRangesList`, one element per group, with metadata columns `count`, `cpm`
#'   and `fold`; `metadata()` holds `members` (barcodes per group with fragments), `tile` and
#'   `match` (barcode matching)
#' @seealso [annotateGroups()], [peaksToGroups()]
#' @examples
#' fr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(101, 151, 5101, 5201), width = 50),
#'                              barcode = c("a-1", "b-1", "c-1", "c-1"), count = 1L)
#' fragmentsToGroups(fr, c(a = "A", b = "A", c = "B"), tile = 1000, min_fragments = 1)
#' @export
fragmentsToGroups <- function(fragments, groups, tile = 500, min_fragments = 5, specific = TRUE,
                              fold = 2, max_regions = 5000) {
    if (is.character(fragments)) fragments <- readFragments(fragments)
    map <- .as_group_map(groups)
    m <- .match_barcodes(fragments$barcode, map)
    message("Matched ", m$matched, " of ", m$total, " barcodes of `groups` to fragments (",
            m$how, ").")
    keep <- !is.na(m$group)
    if (!any(keep)) stop("no fragment barcode matches `groups`")
    fr <- fragments[keep]
    grp <- factor(m$group[keep])
    bin <- (BiocGenerics::start(fr) - 1L) %/% tile
    key <- paste(as.character(GenomicRanges::seqnames(fr)), bin, sep = ":")
    tk <- unique(key)
    ti <- match(key, tk)
    ng <- nlevels(grp)
    counts <- matrix(tabulate(ti + length(tk) * (as.integer(grp) - 1L), nbins = length(tk) * ng),
                     ncol = ng, dimnames = list(NULL, levels(grp)))
    tot <- colSums(counts)
    chr <- sub(":[0-9]+$", "", tk)
    st <- as.integer(sub("^.*:", "", tk)) * tile + 1L
    members <- vapply(levels(grp), function(g) length(unique(fr$barcode[grp == g])), 0L)
    out <- lapply(levels(grp), function(g) {
        cg <- counts[, g]
        cpm <- cg / tot[g] * 1e6
        rest <- rowSums(counts[, setdiff(levels(grp), g), drop = FALSE])
        rest_cpm <- if (ng > 1) rest / sum(tot[setdiff(levels(grp), g)]) * 1e6 else 0
        f <- (cpm + 1) / (rest_cpm + 1)
        ok <- cg >= min_fragments & (!specific | ng == 1 | f >= fold)
        o <- which(ok)[order(-cpm[ok])]
        o <- utils::head(o, max_regions)
        GenomicRanges::GRanges(chr[o], IRanges::IRanges(st[o], width = tile),
                               count = cg[o], cpm = cpm[o], fold = f[o])
    })
    res <- GenomicRanges::GRangesList(stats::setNames(out, levels(grp)))
    S4Vectors::metadata(res) <- list(members = members, tile = tile, match = m[c("matched", "total", "how")])
    res
}

#' Region sets of groups from per-group peaks or marker-peak tables
#'
#' Turns the peaks of each cluster or domain into one region set per group for
#' [annotateGroups()]. Accepted: a named list of `GRanges`; a `GRanges` with a group column; or
#' a data.frame such as the marker table of ArchR `getMarkers()` (`seqnames`, `start`, `end`),
#' a table with `chr`, `start`, `end`, or Signac `FindMarkers()` / 10x differential tables whose
#' peaks are written `"chr1-100-200"` or `"chr1:100-200"` in a `feature` or `peak` column or in
#' the row names. The wide differential table of Cell Ranger ATAC
#' (`analysis/enrichment/graphclust/differential_expression.csv`: one row per peak,
#' `Cluster N Log2 fold change` and `Cluster N Adjusted p value` columns) is read as well; a peak
#' then belongs to every cluster where it passes `min_log2fc` and `max_padj`.
#'
#' @param x the peaks, see Description
#' @param group_col name of the group column; `NULL` looks for `cluster`, `group`,
#'   `group_name` or `Cluster`
#' @param min_log2fc,max_padj thresholds for the wide Cell Ranger ATAC table
#' @return a named `GRangesList`, one element per group; other columns are kept (for the Cell
#'   Ranger ATAC table: `log2fc` and `padj` of that cluster)
#' @seealso [annotateGroups()], [fragmentsToGroups()]
#' @examples
#' df <- data.frame(feature = c("chr1-100-200", "chr1-500-600", "chr2-50-150"),
#'                  cluster = c("1", "1", "2"))
#' peaksToGroups(df)
#' @export
peaksToGroups <- function(x, group_col = NULL, min_log2fc = 1, max_padj = 0.05) {
    if (is.list(x) && !is.data.frame(x) && !is(x, "GRanges")) {
        if (is.null(names(x)) || any(!nzchar(names(x)))) stop("a list of peaks must be named by group")
        return(GenomicRanges::GRangesList(lapply(x, function(p) if (is(p, "GRanges")) p else loadPeak(p))))
    }
    df <- if (is(x, "GRanges")) as.data.frame(x) else as.data.frame(x, check.names = FALSE)
    fc_cols <- grep("^Cluster .+ Log2 fold change$", names(df), value = TRUE)
    if (is.null(group_col) && length(fc_cols)) df <- .cellranger_long(df, fc_cols, min_log2fc, max_padj)
    if (is.null(group_col)) group_col <- intersect(c("cluster", "group", "group_name", "Cluster"), names(df))[1]
    if (is.na(group_col) || !group_col %in% names(df)) stop("no group column found; give `group_col`")
    if (all(c("seqnames", "start", "end") %in% names(df))) {
        chr <- df$seqnames; st <- df$start; en <- df$end
    } else if (all(c("chr", "start", "end") %in% names(df))) {
        chr <- df$chr; st <- df$start; en <- df$end
    } else {
        id <- if ("feature" %in% names(df)) df$feature else if ("peak" %in% names(df)) df$peak else rownames(df)
        p <- regmatches(id, regexec("^(.+)[-:]([0-9]+)-([0-9]+)$", id))
        bad <- lengths(p) != 4
        if (all(bad)) stop("could not find peak coordinates: give seqnames/start/end, chr/start/end, ",
                           "or peaks written chr-start-end in a feature or peak column or the row names")
        if (any(bad)) warning(sum(bad), " rows without parseable coordinates were dropped")
        df <- df[!bad, , drop = FALSE]; p <- p[!bad]
        chr <- vapply(p, `[`, "", 2); st <- as.integer(vapply(p, `[`, "", 3)); en <- as.integer(vapply(p, `[`, "", 4))
    }
    drop_cols <- intersect(names(df), c("seqnames", "chr", "start", "end", "width", "strand"))
    extra <- df[, setdiff(names(df), drop_cols), drop = FALSE]
    gr <- GenomicRanges::GRanges(as.character(chr), IRanges::IRanges(as.integer(st), as.integer(en)))
    S4Vectors::mcols(gr) <- extra
    grp <- as.character(df[[group_col]])
    split(gr, factor(grp, levels = unique(grp)))
}

## wide Cell Ranger ATAC differential table -> long: one row per peak and passing cluster
.cellranger_long <- function(df, fc_cols, min_log2fc, max_padj) {
    id <- df[[intersect(c("Feature ID", "Feature Name"), names(df))[1]]]
    rows <- lapply(fc_cols, function(fc) {
        cl <- sub("^Cluster (.+) Log2 fold change$", "\\1", fc)
        pc <- paste0("Cluster ", cl, " Adjusted p value")
        padj <- if (pc %in% names(df)) df[[pc]] else rep(NA_real_, nrow(df))
        ok <- !is.na(df[[fc]]) & df[[fc]] >= min_log2fc & (is.na(padj) | padj <= max_padj)
        if (!any(ok)) return(NULL)
        data.frame(feature = id[ok], cluster = cl, log2fc = df[[fc]][ok], padj = padj[ok],
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out)) stop("no peak passes min_log2fc and max_padj in any cluster")
    out$cluster <- factor(out$cluster, levels = unique(sub("^Cluster (.+) Log2 fold change$", "\\1", fc_cols)))
    out <- out[order(out$cluster), ]
    out$cluster <- as.character(out$cluster)
    out
}

#' Annotate the region sets of groups of cells or spots
#'
#' Calls [annotateSeq()] on the region set of every group and returns a [csAnnoList], so that
#' clusters of single cells or spatial domains are compared, tested and interpreted like bulk
#' samples. What each group is -- `resolution`, number of cells or spots (`members`) and, for
#' spatial data, the spot coordinates -- is kept in `metadata(<anno>)$group` of every element
#' and read with [groupInfo()]; `interpret()` states it in its evidence digest.
#'
#' @param x a named `GRangesList` or list of `GRanges` (from [fragmentsToGroups()] or
#'   [peaksToGroups()])
#' @param ... passed to [annotateSeq()], e.g. `TxDb`, `annoDb`, `tssRegion`
#' @param resolution `"single-cell"`, `"spatial"` or `"bulk"`
#' @param members optional named numeric vector, cells or spots per group; taken from
#'   `metadata(x)$members` when `x` comes from [fragmentsToGroups()]
#' @param coords optional data.frame of spot coordinates with columns `barcode`, `x`, `y` and
#'   `group`; the rows of each group are kept with it
#' @return a `csAnnoList`
#' @seealso [groupInfo()], [plotSpatialGroups()]
#' @examples
#' \dontrun{
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' groups <- fragmentsToGroups("fragments.tsv.gz", clusters)   # clusters: barcode -> cluster
#' x <- annotateGroups(groups, TxDb = txdb, annoDb = "org.Hs.eg.db")
#' groupInfo(x); plotAnnoBar(x)
#' }
#' @export
annotateGroups <- function(x, ..., resolution = c("single-cell", "spatial", "bulk"), members = NULL,
                           coords = NULL) {
    resolution <- match.arg(resolution)
    if (is.null(members) && is(x, "GRangesList")) members <- S4Vectors::metadata(x)$members
    x <- as.list(x)
    if (is.null(names(x)) || any(!nzchar(names(x)))) stop("the groups of `x` must be named")
    empty <- vapply(x, length, 0L) == 0
    if (any(empty)) {
        message("Skipping groups without regions: ", paste(names(x)[empty], collapse = ", "))
        x <- x[!empty]
    }
    if (!length(x)) stop("no group has regions")
    res <- lapply(names(x), function(g) {
        a <- annotateSeq(x[[g]], ..., verbose = FALSE)
        S4Vectors::metadata(a@anno)$group <- list(
            name = g, resolution = resolution,
            n_members = if (!is.null(members) && g %in% names(members)) unname(members[[g]]) else NA_real_,
            coords = if (!is.null(coords) && "group" %in% names(coords)) coords[coords$group == g, , drop = FALSE] else NULL
        )
        a
    })
    csAnnoList(stats::setNames(res, names(x)))
}

#' What the groups of a csAnnoList are
#'
#' @param x a `csAnnoList` (or a `csAnno`) from [annotateGroups()]
#' @return a data.frame with `group`, `resolution`, `n_members` and `n_regions`; `NA` where not
#'   recorded (for example objects not made by [annotateGroups()])
#' @examples
#' data(peakAnnoList)
#' groupInfo(csAnnoList(peakAnnoList))
#' @export
groupInfo <- function(x) {
    if (is(x, "csAnno")) x <- list(x)
    nm <- names(x); if (is.null(nm)) nm <- rep(NA_character_, length(x))
    rows <- lapply(seq_along(x), function(i) {
        g <- S4Vectors::metadata(x[[i]]@anno)$group
        data.frame(group = if (is.null(g$name)) nm[i] else g$name,
                   resolution = if (is.null(g$resolution)) NA_character_ else g$resolution,
                   n_members = if (is.null(g$n_members)) NA_real_ else g$n_members,
                   n_regions = length(x[[i]]@anno), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' Read spot positions of a spatial assay
#'
#' Reads a Visium-style `tissue_positions_list.csv` (`barcode, in_tissue, row, col, pxl_row,
#' pxl_col`, with or without a header line), as written by Space Ranger and by spatial ATAC /
#' CUT&Tag pipelines.
#'
#' @param file path of the positions file
#' @return a data.frame with `barcode`, `in_tissue`, `row`, `col`, `x` and `y`; `x`, `y` are the
#'   pixel coordinates when given, otherwise `col` and `row`
#' @examples
#' f <- tempfile(fileext = ".csv")
#' writeLines(c("AAAC,1,0,0,100,200", "AAAG,1,0,1,100,210"), f)
#' readSpotPositions(f)
#' @export
readSpotPositions <- function(file) {
    first <- readLines(file, n = 1)
    header <- !grepl("[0-9]", strsplit(first, ",")[[1]][2])
    df <- utils::read.csv(file, header = header, stringsAsFactors = FALSE)
    if (ncol(df) < 4) stop("expected at least barcode, in_tissue, row, col")
    out <- data.frame(barcode = as.character(df[[1]]), in_tissue = as.integer(df[[2]]),
                      row = as.numeric(df[[3]]), col = as.numeric(df[[4]]), stringsAsFactors = FALSE)
    if (ncol(df) >= 6) { out$x <- as.numeric(df[[6]]); out$y <- as.numeric(df[[5]]) } else {
        out$x <- out$col; out$y <- out$row
    }
    out
}

#' Plot spots of a spatial assay by group or value
#'
#' @param coords data.frame with `barcode`, `x`, `y` (e.g. from [readSpotPositions()]), and
#'   optionally `group`
#' @param groups optional barcode -> group (named vector or data.frame, as in
#'   [fragmentsToGroups()]); overrides `coords$group`
#' @param value optional named numeric vector by barcode, e.g. fragments per spot; colours the
#'   spots instead of the groups
#' @param size point size
#' @return a `ggplot`
#' @importFrom rlang .data
#' @examples
#' pos <- data.frame(barcode = c("a", "b", "c"), x = 1:3, y = 1:3)
#' plotSpatialGroups(pos, groups = c(a = "A", b = "A", c = "B"))
#' @export
plotSpatialGroups <- function(coords, groups = NULL, value = NULL, size = 1.5) {
    df <- coords
    if (!is.null(groups)) {
        map <- .as_group_map(groups)
        m <- .match_barcodes(df$barcode, map)
        df$group <- m$group
    }
    if (!is.null(value)) {
        df$value <- unname(value[df$barcode])
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$value)) +
            ggplot2::scale_colour_viridis_c(name = NULL, na.value = "grey85")
    } else {
        if (is.null(df$group)) stop("give `groups`, `value`, or a `group` column in `coords`")
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
            ggplot2::labs(colour = NULL)
    }
    p + ggplot2::geom_point(size = size) + ggplot2::scale_y_reverse() + ggplot2::coord_fixed() +
        ggplot2::theme_void()
}
