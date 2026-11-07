## Single-molecule modification data: an adapter from per-read tool output to epiSeeker.
##
## Reads come in as the tables that modkit and fibertools write, not as BAM/CRAM: decoding the
## MM/ML tags (skip encoding, probabilities, reverse strand, soft clips, CIGAR, and the old and
## new fibertools tags for nucleosomes and MSPs) is what those tools already do well. Per-read
## calls become a GRanges, one row per call, with the read in the mcols and in metadata(); they
## summarise to the existing bmData, so getBmMatrix() and plotBmProf() work unchanged.

## call types: base-level modifications and fibre-level intervals
.sm_base_types <- c("5mC", "5hmC", "m6A")
.sm_interval_types <- c("nuc", "msp")

## modkit mod codes to call types
.modkit_codes <- c(m = "5mC", h = "5hmC", a = "m6A", "21839" = "4mC")

#' Read per-read modification calls
#'
#' Reads the per-read table of `ft extract --all` (fibertools, for Fiber-seq) or of
#' `modkit extract full` / `modkit extract calls` (any modBAM, ONT or PacBio) into a `GRanges`
#' with one row per call. Base-level calls (`5mC`, `5hmC`, `m6A`) are one base wide; for
#' Fiber-seq, nucleosomes (`nuc`) and methyltransferase-sensitive patches (`msp`) are intervals,
#' and an MSP whose FIRE score is above 230 (of 255) is flagged as `modified`.
#'
#' Coordinates of both tools are 0-based and are returned 1-based. Calls without a reference
#' position (`-1`, unmapped or in an insertion) are dropped and counted in
#' `metadata(x)$modSource$n_unplaced`.
#'
#' @param file path to the table (plain or gzipped).
#' @param format `"fibertools"`, `"modkit"` or `"auto"` (from the header).
#' @param region optional `GRanges`; only reads overlapping it, and only their calls overlapping
#'   it, are kept.
#' @param min_prob calls with a probability of at least this are `modified`. For `modkit extract
#'   calls` the tool's own call (`call_code`) is used instead, and calls flagged `fail` are
#'   dropped.
#' @return a `GRanges` with mcols `read_id`, `type`, `prob`, `modified` and `HP` (haplotype, `NA`
#'   when absent); `metadata(x)$reads` has one row per read (`read_id`, `chr`, `start`, `end`,
#'   `strand`, `HP`) and `metadata(x)$modSource` the file, format and threshold.
#' @seealso [summarizeModCalls()], [plotModReads()], [modCooccurrence()]
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c(
#'     "#ct\tst\ten\tfiber\tstrand\tHP\tref_m6a\tm6a_qual\tref_5mC\t5mC_qual\tref_nuc_starts\tref_nuc_lengths\tref_msp_starts\tref_msp_lengths\tfire",
#'     "chr1\t100\t400\tr1\t+\t1\t120,130,150\t250,250,250\t200\t240\t300\t80\t110\t60\t240"
#' ), f)
#' x <- readModCalls(f)
#' table(x$type)
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom S4Vectors metadata metadata<- mcols
#' @export
readModCalls <- function(file, format = c("auto", "fibertools", "modkit"), region = NULL,
                         min_prob = 0.5) {
    format <- match.arg(format)
    df <- .sm_read_table(file)
    if (format == "auto") {
        format <- if ("read_id" %in% names(df) && any(c("ref_position", "mod_qual", "call_prob") %in% names(df))) {
            "modkit"
        } else if (all(c("ct", "st", "en", "fiber") %in% names(df))) {
            "fibertools"
        } else {
            stop("cannot tell the format from the header; give `format`.")
        }
    }
    res <- if (format == "fibertools") .sm_from_fibertools(df, min_prob) else .sm_from_modkit(df, min_prob)
    gr <- res$calls
    reads <- res$reads
    if (!is.null(region)) {
        rg <- GRanges(reads$chr, IRanges(reads$start, reads$end))
        reads <- reads[overlapsAny(rg, region, ignore.strand = TRUE), , drop = FALSE]
        gr <- gr[gr$read_id %in% reads$read_id & overlapsAny(gr, region, ignore.strand = TRUE)]
    }
    metadata(gr)$reads <- reads
    metadata(gr)$modSource <- list(file = normalizePath(file), format = format, min_prob = min_prob,
                                   n_unplaced = res$n_unplaced,
                                   region = if (is.null(region)) NULL else as.character(region))
    gr
}

## read a tab-separated table whose header may start with "#"
.sm_read_table <- function(file) {
    con <- if (grepl("\\.gz$", file)) gzfile(file) else file(file)
    header <- readLines(con, n = 1)
    close(con)
    cols <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1]]
    df <- if (requireNamespace("data.table", quietly = TRUE)) {
        as.data.frame(data.table::fread(file, sep = "\t", header = FALSE, skip = 1, colClasses = "character",
                                        col.names = cols, na.strings = NULL, showProgress = FALSE))
    } else {
        utils::read.delim(file, header = FALSE, skip = 1, col.names = cols, colClasses = "character",
                          check.names = FALSE, quote = "", comment.char = "", na.strings = NULL)
    }
    names(df) <- cols
    df
}

## "1,2,3" -> list of integer vectors; "." or "" -> integer(0)
.sm_split <- function(v) {
    lapply(strsplit(ifelse(is.na(v) | v %in% c(".", ""), "", v), ",", fixed = TRUE), as.numeric)
}

.sm_from_fibertools <- function(df, min_prob) {
    hp <- if ("HP" %in% names(df)) suppressWarnings(as.integer(df$HP)) else rep(NA_integer_, nrow(df))
    strand <- if ("strand" %in% names(df)) df$strand else rep("*", nrow(df))
    reads <- data.frame(read_id = df$fiber, chr = df$ct, start = as.numeric(df$st) + 1,
                        end = as.numeric(df$en), strand = strand, HP = hp, stringsAsFactors = FALSE)
    out <- list()
    n_unplaced <- 0
    ## base-level calls: ref positions with a 0-255 quality
    for (b in list(c("ref_m6a", "m6a_qual", "m6A"), c("ref_5mC", "5mC_qual", "5mC"))) {
        if (!all(b[1:2] %in% names(df))) next
        pos <- .sm_split(df[[b[1]]])
        q <- .sm_split(df[[b[2]]])
        n <- lengths(pos)
        if (!sum(n)) next
        p <- unlist(pos)
        keep <- p >= 0
        n_unplaced <- n_unplaced + sum(!keep)
        prob <- unlist(q) / 255
        out[[b[3]]] <- data.frame(chr = rep(df$ct, n)[keep], start = p[keep] + 1, end = p[keep] + 1,
                                  strand = rep(strand, n)[keep], read_id = rep(df$fiber, n)[keep],
                                  type = b[3], prob = prob[keep], modified = prob[keep] >= min_prob,
                                  HP = rep(hp, n)[keep], stringsAsFactors = FALSE)
    }
    ## intervals: nucleosomes, MSPs (FIRE score 0-255 on MSPs)
    for (iv in list(c("ref_nuc_starts", "ref_nuc_lengths", "nuc"), c("ref_msp_starts", "ref_msp_lengths", "msp"))) {
        if (!all(iv[1:2] %in% names(df))) next
        st <- .sm_split(df[[iv[1]]])
        ln <- .sm_split(df[[iv[2]]])
        n <- lengths(st)
        if (!sum(n)) next
        s <- unlist(st)
        l <- unlist(ln)
        keep <- s >= 0 & l > 0
        n_unplaced <- n_unplaced + sum(!keep)
        fire <- if (iv[3] == "msp" && "fire" %in% names(df)) unlist(.sm_split(df$fire)) else rep(NA_real_, length(s))
        if (length(fire) != length(s)) fire <- rep(NA_real_, length(s))
        prob <- fire / 255
        out[[iv[3]]] <- data.frame(chr = rep(df$ct, n)[keep], start = s[keep] + 1, end = s[keep] + l[keep],
                                   strand = rep(strand, n)[keep], read_id = rep(df$fiber, n)[keep],
                                   type = iv[3], prob = prob[keep],
                                   modified = if (iv[3] == "msp") !is.na(fire[keep]) & fire[keep] > 230 else NA,
                                   HP = rep(hp, n)[keep], stringsAsFactors = FALSE)
    }
    list(calls = .sm_to_granges(do.call(rbind, out)), reads = reads, n_unplaced = n_unplaced)
}

.sm_from_modkit <- function(df, min_prob) {
    pos <- as.numeric(df$ref_position)
    keep <- !is.na(pos) & pos >= 0
    n_unplaced <- sum(!keep)
    df <- df[keep, , drop = FALSE]
    pos <- pos[keep]
    if ("call_code" %in% names(df)) {
        ## extract calls: the tool's own call; "-" is canonical; failed calls are dropped. The
        ## table has no mod_code: a modified call is named by its call code, a canonical one by
        ## the modification usually tested on its base (C -> 5mC, A -> m6A)
        ok <- if ("fail" %in% names(df)) !(tolower(df$fail) %in% c("true", "1")) else rep(TRUE, nrow(df))
        df <- df[ok, , drop = FALSE]
        pos <- pos[ok]
        prob <- as.numeric(df$call_prob)
        modified <- df$call_code != "-"
        by_code <- unname(.modkit_codes[df$call_code])
        by_code[is.na(by_code)] <- df$call_code[is.na(by_code)]
        base <- if ("canonical_base" %in% names(df)) df$canonical_base else rep(NA_character_, nrow(df))
        by_base <- c(C = "5mC", A = "m6A")[base]
        type <- ifelse(modified, by_code, ifelse(is.na(by_base), "canonical", by_base))
    } else {
        ## extract full: one row per tested base and modification, with its probability
        code <- df$mod_code
        type <- unname(.modkit_codes[code])
        type[is.na(type)] <- code[is.na(type)]
        prob <- as.numeric(df$mod_qual)
        modified <- prob >= min_prob
    }
    strand <- if ("ref_strand" %in% names(df)) df$ref_strand else "*"
    calls <- data.frame(chr = df$chrom, start = pos + 1, end = pos + 1, strand = strand,
                        read_id = df$read_id, type = unname(type), prob = prob, modified = modified,
                        HP = NA_integer_, stringsAsFactors = FALSE)
    ## a read spans its alignment when the table says so, otherwise its first to last call
    if (all(c("alignment_start", "alignment_end") %in% names(df))) {
        r <- df[!duplicated(df$read_id), , drop = FALSE]
        reads <- data.frame(read_id = r$read_id, chr = r$chrom, start = as.numeric(r$alignment_start) + 1,
                            end = as.numeric(r$alignment_end), strand = strand[!duplicated(df$read_id)],
                            HP = NA_integer_, stringsAsFactors = FALSE)
    } else {
        s <- split(calls$start, calls$read_id)
        first <- calls[!duplicated(calls$read_id), , drop = FALSE]
        reads <- data.frame(read_id = names(s), chr = first$chr[match(names(s), first$read_id)],
                            start = vapply(s, min, 0), end = vapply(s, max, 0),
                            strand = first$strand[match(names(s), first$read_id)], HP = NA_integer_,
                            stringsAsFactors = FALSE)
    }
    list(calls = .sm_to_granges(calls), reads = reads, n_unplaced = n_unplaced)
}

.sm_to_granges <- function(d) {
    if (is.null(d) || !nrow(d)) {
        return(GRanges(read_id = character(0), type = character(0), prob = numeric(0),
                       modified = logical(0), HP = integer(0)))
    }
    d$strand[!d$strand %in% c("+", "-")] <- "*"
    GRanges(d$chr, IRanges(d$start, d$end), strand = d$strand, read_id = d$read_id, type = d$type,
            prob = d$prob, modified = d$modified, HP = d$HP)
}

#' Read a modkit pileup bedMethyl file as bmData
#'
#' The summarised level: one row per position, `Cov` = valid coverage (`Nvalid_cov`) and
#' `Methylation` = fraction modified (percent / 100), the value columns of the bundled
#' `demo_bmdata`, so [getBmMatrix()] and [plotBmProf()] work on the result.
#'
#' @param file a bedMethyl file from `modkit pileup` (plain or gzipped).
#' @param mod modification code to keep (`"m"` for 5mC, `"h"` for 5hmC, `"a"` for 6mA); `NULL`
#'   keeps all and requires that there is only one.
#' @param sampleName sample name; the file name by default.
#' @return a `bmData`.
#' @examples
#' f <- tempfile(fileext = ".bed")
#' writeLines(c(
#'     "chr1 99 100 m 10 + 99 100 255,0,0 10 80.00 8 2 0 0 0 0 0",
#'     "chr1 149 150 m 5 + 149 150 255,0,0 5 20.00 1 4 0 0 0 0 0"
#' ), f)
#' readBedMethyl(f, mod = "m", sampleName = "s1")
#' @export
readBedMethyl <- function(file, mod = NULL, sampleName = NULL) {
    df <- utils::read.table(file, header = FALSE, sep = "", colClasses = "character", comment.char = "",
                            quote = "", fill = TRUE)
    if (ncol(df) < 11) stop("expected at least 11 columns (modkit pileup bedMethyl).")
    code <- sub(",.*", "", df[[4]])
    if (is.null(mod)) {
        if (length(unique(code)) > 1) stop("several modification codes (", paste(unique(code), collapse = ", "),
                                          "); choose one with `mod`.")
    } else {
        df <- df[code %in% mod, , drop = FALSE]
        if (!nrow(df)) stop("no rows with modification code ", paste(mod, collapse = ", "), ".")
    }
    if (is.null(sampleName)) sampleName <- sub("\\.(bed|bedmethyl)(\\.gz)?$", "", basename(file))
    d <- data.frame(chr = df[[1]], pos = as.numeric(df[[2]]) + 1, Cov = as.numeric(df[[10]]),
                    Methylation = as.numeric(df[[11]]) / 100, stringsAsFactors = FALSE)
    d <- d[!duplicated(d[, c("chr", "pos")]), , drop = FALSE]
    makeBmDataFromData(data = list(d), sampleNames = sampleName)
}

## does the table list every candidate base (calls), or only modified ones (positives)?
.sm_lists_negatives <- function(x) any(!x$modified, na.rm = TRUE)

#' Summarise per-read calls to bmData
#'
#' Counts, for each position with a call of `type`, how many reads cover it (`Cov`) and the
#' fraction of them modified (`Methylation`). When the tool lists every tested base (5mC from
#' PacBio or modkit), coverage is the number of calls at the position; when it lists only
#' modified bases (Fiber-seq m6A), coverage is the number of reads spanning the position.
#'
#' @param x per-read calls from [readModCalls()].
#' @param type a base-level type, e.g. `"5mC"` or `"m6A"`.
#' @param coverage `"auto"` (calls when unmodified calls are present, otherwise reads), `"calls"`
#'   or `"reads"`.
#' @param sampleName sample name.
#' @return a `bmData` with `Cov` and `Methylation`.
#' @export
summarizeModCalls <- function(x, type = "5mC", coverage = c("auto", "calls", "reads"), sampleName = "sample") {
    coverage <- match.arg(coverage)
    if (!type %in% .sm_base_types && !type %in% x$type) stop("`type` must be a base-level call type in `x`.")
    y <- x[x$type == type]
    if (!length(y)) stop("no calls of type ", type, ".")
    if (coverage == "auto") coverage <- if (.sm_lists_negatives(y)) "calls" else "reads"
    key <- paste(as.character(GenomicRanges::seqnames(y)), GenomicRanges::start(y), sep = ":")
    pos <- y[!duplicated(key)]
    nmod <- as.vector(tapply(y$modified, factor(key, levels = unique(key)), sum))
    if (coverage == "calls") {
        cov <- as.vector(table(factor(key, levels = unique(key))))
    } else {
        reads <- metadata(x)$reads
        rg <- GRanges(reads$chr, IRanges(reads$start, reads$end))
        cov <- GenomicRanges::countOverlaps(pos, rg, ignore.strand = TRUE)
    }
    d <- data.frame(chr = as.character(GenomicRanges::seqnames(pos)), pos = GenomicRanges::start(pos),
                    Cov = cov, Methylation = nmod / cov, stringsAsFactors = FALSE)
    makeBmDataFromData(data = list(d), sampleNames = sampleName)
}

#' Per-read values over regions
#'
#' One value per read and region, for reads that span the whole region (`NA` otherwise): for an
#' interval type (`"msp"`, `"nuc"`) the fraction of the region the read covers with it (with
#' `fire_only = TRUE`, only MSPs that are FIREs); for a base-level type the fraction of calls
#' modified when the table lists unmodified calls, otherwise modified calls per kb.
#'
#' @param x per-read calls from [readModCalls()].
#' @param regions a `GRanges`; names, if present, become column names.
#' @param type call type.
#' @param fire_only for `type = "msp"`, count only FIREs.
#' @return a numeric matrix, reads by regions.
#' @export
modReadMatrix <- function(x, regions, type = "msp", fire_only = FALSE) {
    reads <- metadata(x)$reads
    rg <- GRanges(reads$chr, IRanges(reads$start, reads$end))
    y <- x[x$type == type]
    if (type == "msp" && fire_only) y <- y[!is.na(y$modified) & y$modified]
    interval <- type %in% .sm_interval_types
    negatives <- !interval && .sm_lists_negatives(y)
    m <- matrix(NA_real_, nrow = length(rg), ncol = length(regions),
                dimnames = list(reads$read_id, if (is.null(names(regions))) as.character(regions) else names(regions)))
    for (j in seq_along(regions)) {
        r <- regions[j]
        spans <- which(GenomicRanges::start(rg) <= GenomicRanges::start(r) & GenomicRanges::end(rg) >= GenomicRanges::end(r) &
                           as.character(GenomicRanges::seqnames(rg)) == as.character(GenomicRanges::seqnames(r)))
        if (!length(spans)) next
        yr <- y[y$read_id %in% reads$read_id[spans] & overlapsAny(y, r, ignore.strand = TRUE)]
        for (i in spans) {
            yi <- yr[yr$read_id == reads$read_id[i]]
            m[i, j] <- if (interval) {
                sum(GenomicRanges::width(GenomicRanges::intersect(GenomicRanges::reduce(yi, ignore.strand = TRUE), r,
                                                                  ignore.strand = TRUE))) / GenomicRanges::width(r)
            } else if (negatives) {
                if (length(yi)) mean(yi$modified) else NA_real_
            } else {
                1000 * sum(yi$modified) / GenomicRanges::width(r)
            }
        }
    }
    m
}

#' Do two regions open on the same molecules?
#'
#' For the reads that span both regions, whether each region is open (or modified) on that read
#' -- value of [modReadMatrix()] at least `cutoff` -- crossed in a 2-by-2 table and tested with
#' Fisher's exact test. A positive association means the two regions tend to open together on
#' the same molecule, which a bulk profile cannot show.
#'
#' @param x per-read calls from [readModCalls()].
#' @param region1,region2 `GRanges` of length one.
#' @param type,fire_only passed to [modReadMatrix()].
#' @param cutoff value at or above which a region counts as open on a read.
#' @return a list with `table`, `n_reads`, `odds_ratio` and `p_value` (`NA` when fewer than 2
#'   reads span both regions or a margin is empty).
#' @export
modCooccurrence <- function(x, region1, region2, type = "msp", fire_only = FALSE, cutoff = 0.5) {
    m <- modReadMatrix(x, c(region1, region2), type = type, fire_only = fire_only)
    ok <- stats::complete.cases(m)
    a <- factor(m[ok, 1] >= cutoff, levels = c(TRUE, FALSE), labels = c("open", "closed"))
    b <- factor(m[ok, 2] >= cutoff, levels = c(TRUE, FALSE), labels = c("open", "closed"))
    tab <- table(region1 = a, region2 = b)
    res <- list(table = tab, n_reads = sum(ok), odds_ratio = NA_real_, p_value = NA_real_)
    if (sum(ok) >= 2 && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
        ft <- stats::fisher.test(tab)
        res$odds_ratio <- unname(ft$estimate)
        res$p_value <- ft$p.value
    }
    res
}

#' Plot reads over a region
#'
#' One row per read over `region`: the read extent in grey, nucleosomes, MSPs (FIREs darker),
#' and modified bases as points; optional `peaks` (for example FIRE peaks or ChIP peaks) are
#' drawn as a band above the reads.
#'
#' @param x per-read calls from [readModCalls()].
#' @param region a `GRanges` of length one.
#' @param types call types to draw.
#' @param peaks optional `GRanges` drawn above the reads.
#' @param max_reads at most this many reads are drawn (the first by start).
#' @param order_by `"start"`, or `"open"` to sort reads by the fraction of the region covered by
#'   MSPs.
#' @return a ggplot object.
#' @importFrom ggplot2 ggplot aes geom_segment geom_rect geom_point scale_fill_manual scale_colour_manual labs theme_classic theme element_blank
#' @export
plotModReads <- function(x, region, types = c("nuc", "msp", "m6A", "5mC"), peaks = NULL,
                         max_reads = 200, order_by = c("start", "open")) {
    order_by <- match.arg(order_by)
    reads <- metadata(x)$reads
    rg <- GRanges(reads$chr, IRanges(reads$start, reads$end))
    keep <- overlapsAny(rg, region, ignore.strand = TRUE)
    reads <- reads[keep, , drop = FALSE]
    if (!nrow(reads)) stop("no reads overlap `region`.")
    if (order_by == "open") {
        open <- modReadMatrix(x[x$read_id %in% reads$read_id], region, type = "msp")
        reads$score <- open[match(reads$read_id, rownames(open)), 1]
        reads <- reads[order(-reads$score, reads$start, na.last = TRUE), , drop = FALSE]
    } else {
        reads <- reads[order(reads$start), , drop = FALSE]
    }
    reads <- utils::head(reads, max_reads)
    reads$y <- seq_len(nrow(reads))
    lo <- GenomicRanges::start(region)
    hi <- GenomicRanges::end(region)
    y <- x[x$read_id %in% reads$read_id & x$type %in% types & overlapsAny(x, region, ignore.strand = TRUE)]
    d <- as.data.frame(y)
    d$y <- reads$y[match(d$read_id, reads$read_id)]
    d$start <- pmax(d$start, lo)
    d$end <- pmin(d$end, hi)
    ivl <- d[d$type %in% .sm_interval_types, , drop = FALSE]
    ivl$what <- ifelse(ivl$type == "nuc", "nucleosome", ifelse(!is.na(ivl$modified) & ivl$modified, "FIRE", "MSP"))
    pts <- d[d$type %in% .sm_base_types & d$modified %in% TRUE, , drop = FALSE]
    reads$s <- pmax(reads$start, lo)
    reads$e <- pmin(reads$end, hi)
    start <- end <- s <- e <- what <- type <- NULL
    p <- ggplot() +
        geom_segment(data = reads, aes(x = s, xend = e, y = y, yend = y), colour = "grey80", linewidth = 0.3)
    if (nrow(ivl)) {
        p <- p + geom_rect(data = ivl, aes(xmin = start, xmax = end, ymin = y - 0.4, ymax = y + 0.4, fill = what)) +
            scale_fill_manual(values = c(nucleosome = "grey60", MSP = "#9ecae1", FIRE = "#08519c"), name = NULL)
    }
    if (nrow(pts)) {
        p <- p + geom_point(data = pts, aes(x = start, y = y, colour = type), size = 0.3, shape = 15) +
            scale_colour_manual(values = c(m6A = "#e6550d", `5mC` = "#756bb1", `5hmC` = "#31a354"), name = NULL)
    }
    if (!is.null(peaks)) {
        pk <- as.data.frame(peaks[overlapsAny(peaks, region, ignore.strand = TRUE)])
        if (nrow(pk)) {
            top <- nrow(reads) + 1.5
            p <- p + geom_rect(data = pk, aes(xmin = pmax(start, lo), xmax = pmin(end, hi)),
                               ymin = top, ymax = top + 1.5, fill = "black", inherit.aes = FALSE)
        }
    }
    p + labs(x = paste0(as.character(GenomicRanges::seqnames(region)), " position"), y = "reads") +
        theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
