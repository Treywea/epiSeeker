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
## does the table list unmodified calls? `ft extract` keeps only calls at or above its ML filter
## (125 of 255 by default), so a fibertools table never does, even though calls of 125-127 fall
## under min_prob = 0.5
.sm_lists_negatives <- function(x) {
    if (identical(metadata(x)$modSource$format, "fibertools")) return(FALSE)
    any(!x$modified, na.rm = TRUE)
}

#' Summarise per-read calls to bmData
#'
#' Counts, for each position with a call of `type`, how many reads cover it (`Cov`) and the
#' fraction of them modified (`Methylation`). When the tool lists every tested base (modkit),
#' coverage is the number of calls at the position; when it lists only calls above a threshold
#' (fibertools: m6A and 5mC), coverage is the number of reads spanning the position -- for 5mC
#' and 5hmC only reads of the strand that gave calls there (PacBio calls CpGs on the read's own
#' strand, at the C on `+` reads and the G on `-` reads), for m6A reads of both strands (Fiber-seq
#' reads carry m6A of both strands). Positions without any listed call are absent.
#'
#' @param x per-read calls from [readModCalls()].
#' @param type a base-level type, e.g. `"5mC"` or `"m6A"`.
#' @param coverage `"auto"` (reads for fibertools tables; otherwise calls when unmodified calls
#'   are present, else reads), `"calls"` or `"reads"`.
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
        if (type == "m6A" || !all(reads$strand %in% c("+", "-"))) {
            cov <- GenomicRanges::countOverlaps(pos, rg, ignore.strand = TRUE)
        } else {
            ks <- unique(paste(key, as.character(GenomicRanges::strand(y))))
            ku <- key[!duplicated(key)]
            cov <- (paste(ku, "+") %in% ks) * GenomicRanges::countOverlaps(pos, rg[reads$strand == "+"], ignore.strand = TRUE) +
                (paste(ku, "-") %in% ks) * GenomicRanges::countOverlaps(pos, rg[reads$strand == "-"], ignore.strand = TRUE)
        }
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
#' modified when the table lists unmodified calls (modkit), otherwise (fibertools) modified calls
#' per kb.
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
    ## (read i, region j) pairs where the read spans the whole region; values come from overlap
    ## joins keyed by read id and region (a per-read loop took hours on a real 120-kb window)
    sp <- findOverlaps(regions, rg, type = "within", ignore.strand = TRUE)
    if (!length(sp)) return(m)
    si <- subjectHits(sp)
    sj <- queryHits(sp)
    key <- paste(reads$read_id[si], sj, sep = "\r")
    rw <- GenomicRanges::width(regions)
    agg <- function(v, k, f) {
        if (!length(v)) return(stats::setNames(numeric(0), character(0)))
        a <- tapply(v, k, f)
        stats::setNames(as.vector(a), names(a))
    }
    if (interval) {
        yr <- unlist(GenomicRanges::reduce(split(GenomicRanges::granges(y), y$read_id), ignore.strand = TRUE))
        ov <- findOverlaps(yr, regions, ignore.strand = TRUE)
        q <- queryHits(ov)
        s <- subjectHits(ov)
        w <- pmin(GenomicRanges::end(yr)[q], GenomicRanges::end(regions)[s]) -
            pmax(GenomicRanges::start(yr)[q], GenomicRanges::start(regions)[s]) + 1
        v <- agg(w, paste(names(yr)[q], s, sep = "\r"), sum)[key]
        v[is.na(v)] <- 0
        m[cbind(si, sj)] <- v / rw[sj]
    } else {
        ov <- findOverlaps(y, regions, ignore.strand = TRUE)
        k <- paste(y$read_id[queryHits(ov)], subjectHits(ov), sep = "\r")
        mod <- y$modified[queryHits(ov)]
        if (negatives) {
            m[cbind(si, sj)] <- agg(mod, k, mean)[key]
        } else {
            v <- agg(mod, k, sum)[key]
            v[is.na(v)] <- 0
            m[cbind(si, sj)] <- 1000 * v / rw[sj]
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
#' Molecules differ in how much of them is open or labelled at all (e.g. uneven m6A labelling
#' in Fiber-seq): a read with many MSPs everywhere is more likely to be open at both regions,
#' which makes unrelated regions look co-open. On real COLO829BL Fiber-seq, 7 of 44 pairs of
#' FIRE elements were positive by the unadjusted test and 23% of distance-matched non-FIRE
#' pairs as well. `adjust = "read"` tests the association given each read's overall openness
#' (its fraction covered by `type`, or for a base-level type its modified fraction or
#' modifications per kb): by logistic regression of region 2 on region 1 and the covariate
#' (`method = "logistic"`), or by a Cochran-Mantel-Haenszel test stratified by quantiles of the
#' covariate (`method = "cmh"`). By default the covariate leaves out the two regions
#' (`covariate = "outside"`), so it does not contain the outcome; `"all"` uses the whole read;
#' a numeric vector named by read id gives any other per-read covariate. The unadjusted table,
#' odds ratio and p-value are always returned; the default stays `adjust = "none"`, so results
#' made before the option are unchanged -- use `adjust = "read"` before calling regions
#' co-open.
#'
#' @param x per-read calls from [readModCalls()].
#' @param region1,region2 `GRanges` of length one.
#' @param type,fire_only passed to [modReadMatrix()].
#' @param cutoff value at or above which a region counts as open on a read.
#' @param adjust `"none"` or `"read"` (adjust for each read's overall openness).
#' @param method with `adjust = "read"`: `"logistic"` or `"cmh"`.
#' @param covariate with `adjust = "read"`: `"outside"` (the read outside the two regions),
#'   `"all"` (the whole read), or a numeric vector named by read id.
#' @param strata number of quantile strata of the covariate for `method = "cmh"`.
#' @return a list with `table`, `n_reads`, `odds_ratio` and `p_value` (`NA` when fewer than 2
#'   reads span both regions or a margin is empty), and `adjust`; with `adjust = "read"` also
#'   `adj_odds_ratio` and `adj_p_value` (conditional on the covariate), `adjust_method` and
#'   `covariate` (how it was made, and its median and IQR over the reads tested).
#' @export
modCooccurrence <- function(x, region1, region2, type = "msp", fire_only = FALSE, cutoff = 0.5,
                            adjust = c("none", "read"), method = c("logistic", "cmh"),
                            covariate = "outside", strata = 4) {
    adjust <- match.arg(adjust)
    method <- match.arg(method)
    x <- x[x$type == type]  # once, for both the matrix and the covariate (metadata is kept)
    m <- modReadMatrix(x, c(region1, region2), type = type, fire_only = fire_only)
    ok <- stats::complete.cases(m)
    a <- factor(m[ok, 1] >= cutoff, levels = c(TRUE, FALSE), labels = c("open", "closed"))
    b <- factor(m[ok, 2] >= cutoff, levels = c(TRUE, FALSE), labels = c("open", "closed"))
    tab <- table(region1 = a, region2 = b)
    res <- list(table = tab, n_reads = sum(ok), odds_ratio = NA_real_, p_value = NA_real_, adjust = adjust)
    testable <- sum(ok) >= 2 && all(rowSums(tab) > 0) && all(colSums(tab) > 0)
    if (testable) {
        ft <- stats::fisher.test(tab)
        res$odds_ratio <- unname(ft$estimate)
        res$p_value <- ft$p.value
    }
    if (adjust == "none") return(res)

    if (is.numeric(covariate)) {
        if (is.null(names(covariate))) stop("a numeric `covariate` must be named by read id.")
        cv <- unname(covariate[rownames(m)[ok]])
        cv_desc <- "given by the user"
    } else {
        covariate <- match.arg(covariate, c("outside", "all"))
        exclude <- if (covariate == "outside") c(region1, region2) else NULL
        cv <- .sm_read_openness(x, type, fire_only, exclude)[rownames(m)[ok]]
        cv_desc <- paste0(if (.sm_is_interval(type)) "fraction of the read covered by " else
            "per-read level of ", type, if (fire_only) " (FIREs only)", ", ",
            if (covariate == "outside") "outside the two regions" else "whole read")
    }
    res$adjust_method <- method
    res$covariate <- list(how = cv_desc, median = stats::median(cv, na.rm = TRUE),
                          IQR = stats::quantile(cv, c(0.25, 0.75), na.rm = TRUE, names = FALSE))
    res$adj_odds_ratio <- NA_real_
    res$adj_p_value <- NA_real_
    o1 <- a == "open"
    o2 <- b == "open"
    use <- !is.na(cv)
    if (!testable || sum(use) < 3) return(res)
    if (method == "logistic") {
        f <- suppressWarnings(stats::glm(o2[use] ~ o1[use] + cv[use], family = stats::binomial()))
        cf <- summary(f)$coefficients
        if ("o1[use]TRUE" %in% rownames(cf)) {
            res$adj_odds_ratio <- exp(cf["o1[use]TRUE", 1])
            res$adj_p_value <- cf["o1[use]TRUE", 4]
        }
    } else {
        br <- unique(stats::quantile(cv[use], seq(0, 1, length.out = strata + 1), names = FALSE))
        st <- if (length(br) > 1) cut(cv[use], br, include.lowest = TRUE) else factor(rep(1, sum(use)))
        t3 <- table(factor(o1[use], c(TRUE, FALSE)), factor(o2[use], c(TRUE, FALSE)), st)
        t3 <- t3[, , apply(t3, 3, sum) > 1, drop = FALSE]
        if (dim(t3)[3] >= 1) {
            mh <- tryCatch(stats::mantelhaen.test(t3), error = function(e) NULL)
            if (!is.null(mh)) {
                res$adj_odds_ratio <- unname(mh$estimate)
                res$adj_p_value <- mh$p.value
            }
        }
    }
    res
}

.sm_is_interval <- function(type) type %in% .sm_interval_types

## Overall openness of every read: for an interval type the fraction of the read covered by it,
## for a base-level type the modified fraction (tables listing unmodified calls) or modifications
## per kb. `exclude` (GRanges) is left out of both numerator and denominator. When the calls
## were read for a region only, the read is clipped to that region.
.sm_read_openness <- function(x, type, fire_only = FALSE, exclude = NULL) {
    reads <- metadata(x)$reads
    ids <- reads$read_id
    ## each read as its own sequence, so overlaps are only ever within a read
    per_read <- function(g) {
        g <- GenomicRanges::granges(g)
        k <- which(as.character(GenomicRanges::seqnames(g))[rep(seq_along(g), each = length(ids))] ==
                       rep(reads$chr, length(g)))
        gi <- (k - 1) %/% length(ids) + 1
        ri <- (k - 1) %% length(ids) + 1
        GRanges(ids[ri], IRanges(GenomicRanges::start(g)[gi], GenomicRanges::end(g)[gi]))
    }
    rr <- GRanges(ids, IRanges(reads$start, reads$end))
    reg <- metadata(x)$modSource$region
    if (!is.null(reg)) rr <- GenomicRanges::intersect(rr, per_read(GRanges(reg)), ignore.strand = TRUE)
    if (!is.null(exclude)) rr <- GenomicRanges::setdiff(rr, per_read(exclude), ignore.strand = TRUE)
    sp <- tapply(GenomicRanges::width(rr), factor(as.character(GenomicRanges::seqnames(rr)), levels = ids), sum)
    out <- stats::setNames(as.numeric(sp), ids)
    has_span <- !is.na(out) & out > 0
    y <- x[x$type == type]
    if (type == "msp" && fire_only) y <- y[!is.na(y$modified) & y$modified]
    yr <- GRanges(y$read_id, IRanges(GenomicRanges::start(y), GenomicRanges::end(y)))
    rf <- function(v) factor(v, levels = ids)
    if (.sm_is_interval(type)) {
        iv <- GenomicRanges::intersect(GenomicRanges::reduce(yr), rr, ignore.strand = TRUE)
        num <- tapply(GenomicRanges::width(iv), rf(as.character(GenomicRanges::seqnames(iv))), sum)
        num[is.na(num)] <- 0
        val <- as.numeric(num) / out
    } else {
        inside <- overlapsAny(yr, rr, ignore.strand = TRUE)
        if (.sm_lists_negatives(y)) {
            val <- as.numeric(tapply(y$modified[inside], rf(y$read_id[inside]), mean))
        } else {
            num <- tapply(y$modified[inside], rf(y$read_id[inside]), sum)
            num[is.na(num)] <- 0
            val <- 1000 * as.numeric(num) / out
        }
    }
    val[!has_span] <- NA_real_
    stats::setNames(val, ids)
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
#' @param max_reads at most this many reads are drawn: those covering most of `region`.
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
    if (nrow(reads) > max_reads) {
        ## keep the reads covering most of the region, not the first by start (all clipped at its left edge)
        cover <- pmin(reads$end, GenomicRanges::end(region)) - pmax(reads$start, GenomicRanges::start(region))
        reads <- reads[order(-cover, reads$start)[seq_len(max_reads)], , drop = FALSE]
    }
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
