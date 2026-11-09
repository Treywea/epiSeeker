#' Get the information of motif in a range
#'
#' Every region is scanned separately and the hits are placed at the genomic coordinates of
#' the region they were found in (column `region`, the index of the region in `region`).
#' Regions on the `-` strand are scanned on that strand (the sequence `getSeq()` returns), and
#' their hits are mapped back to genomic coordinates, with the strand of the hit flipped.
#' `attr(, "range")` is the start and end of the region for one region, and the span of all
#' regions when they lie on one chromosome (`NA` otherwise); `attr(, "regions")` keeps the
#' regions. Earlier versions recycled the region starts over the rows of the hits, which put
#' the hits of several regions at wrong coordinates.
#'
#' @param region region object in granges (one or more regions).
#' @param pwm PFMatrixList.
#' @param ref_obj seq reference object. e.g. BSgenome object.
#' @param by show the motif by name or ID.
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rlang check_installed
#' @return score matrix
#' @examples
#' require(BSgenome.Hsapiens.UCSC.hg38)
#' data(pwm_obj)
#'
#' region_oi <- GRanges(
#'     seqnames = "chr22",
#'     ranges = IRanges(start = 10525891, end = 10525991)
#' )
#' motifMatrix <- getMotifMatrix(
#'     region = region_oi,
#'     pwm = pwm_obj[c(45, 120, 170)],
#'     ref_obj = BSgenome.Hsapiens.UCSC.hg38
#' )
#'
#' @export
getMotifMatrix <- function(region,
                           pwm,
                           ref_obj,
                           by = "name") {
    rlang::check_installed(c("BSgenome", "universalmotif"), reason = "Matching motif...")

    by <- match.arg(by, c("name", "ID"))

    region <- GenomicRanges::granges(region)
    n_reg <- length(region)
    if (n_reg == 0) stop("`region` contains no region.")
    start_pos <- start(region)
    end_pos <- end(region)
    chr_name <- as.character(seqnames(region))
    minus <- as.character(GenomicRanges::strand(region)) == "-"

    # get seq base
    if (requireNamespace("BSgenome", quietly = TRUE)) {
        regionSeqs <- BSgenome::getSeq(ref_obj, region)
        names(regionSeqs) <- paste0("region", seq_len(n_reg))
    }

    # get motif position
    if (requireNamespace("universalmotif", quietly = TRUE)) {
        motif_df <- as.data.frame(universalmotif::scan_sequences(pwm, regionSeqs, respect.strand = TRUE))
    }

    ## the region every hit was found in; offsets are relative to that region's sequence
    ri <- if (!is.null(motif_df$sequence.i)) as.integer(motif_df$sequence.i) else
        match(as.character(motif_df$sequence), names(regionSeqs))
    motif_df$len <- motif_df$stop - motif_df$start
    to_genome <- function(p) ifelse(minus[ri], end_pos[ri] - p + 1, start_pos[ri] + p - 1)
    motif_df$start <- to_genome(motif_df$start)
    motif_df$stop <- to_genome(motif_df$stop)
    motif_df$strand <- as.character(motif_df$strand)
    flip <- minus[ri] & motif_df$strand %in% c("+", "-")
    motif_df$strand[flip] <- ifelse(motif_df$strand[flip] == "+", "-", "+")
    motif_df$chr <- chr_name[ri]
    motif_df$region <- ri

    if (by == "name") {
        motif_df$index <- motif_df$motif
    } else {
        motif_df$index <- names(pwm)[motif_df$motif.i]
    }

    motif_list <- lapply(seq_len(nrow(motif_df)), function(i) {
        s <- motif_df$start[i]
        e <- motif_df$stop[i]
        sc <- motif_df$score[i]
        st <- motif_df$strand[i]
        mt <- motif_df$index[i]
        seg_id <- i

        coords <- seq.int(s, e)

        df <- data.frame(
            chr        = motif_df$chr[i],
            coordinate = coords,
            score      = sc,
            strand     = st,
            motif      = mt,
            seg_id     = seg_id,
            region     = motif_df$region[i]
        )
        if (st == "-") {
            df$score <- -df$score
        }

        df
    })

    # Efficient rbind using data.table if available
    if (requireNamespace("data.table", quietly = TRUE)) {
        result <- data.table::rbindlist(motif_list)
    } else {
        result <- yulab.utils::rbindlist(motif_list)
    }

    attr(result, "range") <- if (n_reg == 1) c(start_pos, end_pos) else
        if (length(unique(chr_name)) == 1) c(min(start_pos), max(end_pos)) else c(NA_integer_, NA_integer_)
    attr(result, "regions") <- region

    return(result)
}
