## Evidence digest for interpret(<csAnno>): everything the model sees, built without a model.
##
## The order follows what the clusterProfiler interpret layer learned: specific evidence first
## (top regions, genes, overlap with reference samples), summaries after (genomic distribution,
## distance to TSS, width). Nothing is cut silently: every section states how many items it
## sends and how many there are.

## columns tried, in this order, to rank regions when `signal` is NULL
.signal_columns <- c("signalValue", "signal", "score", "fold_enrichment", "V7", "V5")

## columns tried, in this order, to rank reference samples of an overlap table. The p-value of
## enrichPeakOverlap() cannot rank: its smallest value is 1/(nShuffle+1), shared by most samples.
.overlap_rank_columns <- c("jaccard", "frac_query", "fold_enrichment", "N_OL")

.tss_breaks <- c(0, 1e3, 3e3, 5e3, 1e4, 1e5, Inf)

#' @noRd
.seq_digest <- function(x, context = NULL, assay = NULL, question = NULL, overlap = NULL,
                        n_regions = 200, signal = NULL, max_genes = 300, n_distal_genes = 50,
                        n_overlap = 30) {
    if (!is(x, "csAnno")) stop("`x` must be a csAnno object from annotateSeq().")
    gr <- x@anno
    df <- as.data.frame(gr)
    n <- nrow(df)
    if (n == 0) stop("`x` contains no regions.")

    ord <- .rank_regions(df, signal)
    df <- df[ord$index, , drop = FALSE]
    gene <- .region_gene(df)
    feature <- .short_annotation(df$annotation)

    sections <- list()
    sections$overview <- .digest_overview(x, n, context, assay, question)
    sections$top_regions <- .digest_regions(df, gene, feature, ord$column, n_regions)
    genes <- .digest_genes(gene, feature, max_genes, n_distal_genes)
    sections$genes <- genes$text
    ov <- .digest_overlap(overlap, n_overlap)
    if (!is.null(ov)) sections$overlap <- ov$text
    sections$distribution <- .digest_distribution(x)
    sections$tss_distance <- .digest_tss(df$distanceToTSS)
    sections$width <- .digest_width(df$width)

    structure(list(
        text = paste(unlist(sections), collapse = "\n\n"),
        sections = sections,
        sha256 = vapply(sections, .sha256, ""),
        evidence = list(
            regions = utils::head(data.frame(rank = seq_len(n), seqnames = df$seqnames,
                                             start = df$start, end = df$end, gene = gene,
                                             feature = feature, stringsAsFactors = FALSE),
                                  n_regions),
            promoter_genes = genes$promoter,
            distal_genes = genes$distal,
            all_genes = unique(gene[!is.na(gene)]),
            overlap = if (is.null(ov)) NULL else ov$table,
            annoStat = x@annoStat
        ),
        totals = list(n_regions = n, n_regions_sent = min(n, n_regions), rank_by = ord$column)
    ), class = "seq_digest")
}

## rank regions by a signal column (descending) or, without one, by distance to TSS
.rank_regions <- function(df, signal) {
    if (!is.null(signal)) {
        if (!signal %in% names(df)) stop("`signal` column '", signal, "' not found in the regions.")
    } else {
        signal <- intersect(.signal_columns, names(df))[1]
        if (!is.na(signal) && !is.numeric(df[[signal]])) signal <- NA
    }
    if (is.null(signal) || is.na(signal)) {
        return(list(index = order(abs(df$distanceToTSS)), column = "|distanceToTSS| (no signal column)"))
    }
    list(index = order(-df[[signal]], na.last = TRUE), column = signal)
}

## the gene a region is assigned to: symbol when annotated, otherwise the gene id
.region_gene <- function(df) {
    g <- if ("SYMBOL" %in% names(df)) as.character(df$SYMBOL) else rep(NA_character_, nrow(df))
    if ("geneId" %in% names(df)) g[is.na(g)] <- as.character(df$geneId[is.na(g)])
    g
}

## "Intron (uc001clv.1/55182, intron 2 of 14)" -> "Intron"; "Promoter (<=1kb)" is kept whole
.short_annotation <- function(a) {
    a <- as.character(a)
    keep <- grepl("^Promoter", a)
    a[!keep] <- sub(" \\(.*$", "", a[!keep])
    a
}

.digest_overview <- function(x, n, context, assay, question) {
    paste0(
        "## Overview\n",
        "Regions: ", n, " (annotated at ", x@level, " level, TSS region ",
        x@tssRegion[1], " to ", x@tssRegion[2], " bp)\n",
        "Genome: ", .source_line(x), "\n",
        .checks_line(x),
        "Declared assay: ", if (is.null(assay)) "not given" else assay, "\n",
        "Context: ", if (is.null(context)) "not given" else context, "\n",
        "Question: ", if (is.null(question)) "none" else question
    )
}

## What the checks of annotateSeq() found (annoChecks()), as one line of the overview, so that
## a model or an agent reading the digest sees a suspected genome-build mismatch or dropped
## peaks before any other evidence. Empty for objects made before the checks existed.
.checks_line <- function(x, prefix = "Checks: ") {
    chk <- annoChecks(x)
    if (is.null(chk)) return("")
    b <- chk$build; d <- chk$dropped
    parts <- c(
        if (isTRUE(b$suspect)) paste0("GENOME BUILD SUSPECT - ", paste(b$reasons, collapse = "; "),
                                      "; annotation-based evidence below may be wrong") else
        if (isFALSE(b$reliable)) paste0("genome build not contradicted, but the check is unreliable with ",
                                        b$n_peaks, " peaks") else
        if (!is.null(b)) "genome build consistent with TxDb",
        if (!is.null(d) && isTRUE(d$n > 0)) paste0(d$n, " of ", d$n_total, " peaks dropped (",
            paste0(d$by_seqlevel$reason, ": ", d$by_seqlevel$seqlevel, " ", d$by_seqlevel$n,
                   collapse = ", "), ")")
    )
    if (!length(parts)) return("")
    paste0(prefix, paste(parts, collapse = "; "), "\n")
}

## "hg38 (TxDb: UCSC knownGene GENCODE V48, Homo sapiens; annoDb: org.Hs.eg.db 3.23.1)"
.source_line <- function(x) {
    src <- annoSource(x)
    bits <- c(
        if (!is.null(src$TxDb)) paste0("TxDb: ", src$TxDb$source,
                                       if (!is.na(src$TxDb$organism)) paste0(", ", src$TxDb$organism)),
        if (!is.null(src$annoDb)) paste0("annoDb: ", src$annoDb$name,
                                         if (!is.na(src$annoDb$version)) paste0(" ", src$annoDb$version))
    )
    paste0(if (is.na(src$genome)) "not recorded" else src$genome,
           if (length(bits)) paste0(" (", paste(bits, collapse = "; "), ")"))
}

.digest_regions <- function(df, gene, feature, rank_by, n_regions) {
    k <- min(nrow(df), n_regions)
    d <- df[seq_len(k), , drop = FALSE]
    sig <- if (rank_by %in% names(d)) paste0("\t", signif(d[[rank_by]], 4)) else ""
    lines <- paste0(seq_len(k), "\t", d$seqnames, ":", d$start, "-", d$end, "\t", d$width, "\t",
                    feature[seq_len(k)], "\t", gene[seq_len(k)], "\t", d$distanceToTSS, sig)
    header <- paste0("rank\tregion\twidth\tfeature\tnearest_gene\tdistance_to_TSS",
                     if (nzchar(sig[1])) paste0("\t", rank_by) else "")
    paste0("## Top regions (", k, " of ", nrow(df), ", ranked by ", rank_by, ")\n",
           header, "\n", paste(lines, collapse = "\n"))
}

.digest_genes <- function(gene, feature, max_genes, n_distal_genes) {
    prom <- unique(gene[grepl("^Promoter", feature) & !is.na(gene)])
    distal <- gene[!grepl("^Promoter", feature) & !is.na(gene)]
    tab <- sort(table(distal), decreasing = TRUE)
    tab <- tab[tab >= 2]
    k <- min(length(prom), max_genes)
    kd <- min(length(tab), n_distal_genes)
    text <- paste0(
        "## Genes\n",
        "Promoter genes (", k, " of ", length(prom), ", in region rank order): ",
        if (k) paste(prom[seq_len(k)], collapse = ", ") else "none", "\n",
        "Nearest genes of non-promoter regions with >= 2 regions (", kd, " of ", length(tab), "): ",
        if (kd) paste0(names(tab)[seq_len(kd)], " (", as.integer(tab)[seq_len(kd)], ")", collapse = ", ") else "none"
    )
    list(text = text, promoter = prom, distal = names(tab))
}

#' @noRd
.digest_overlap <- function(overlap, n_overlap) {
    if (is.null(overlap)) return(NULL)
    if (!is.data.frame(overlap) || !"tSample" %in% names(overlap)) {
        stop("`overlap` must be the data frame returned by enrichPeakOverlap() (with a tSample column).")
    }
    ov <- overlap
    total <- nrow(ov)
    if ("p.adjust" %in% names(ov)) ov <- ov[!is.na(ov$p.adjust) & ov$p.adjust < 0.05, , drop = FALSE]
    rank_by <- intersect(.overlap_rank_columns, names(ov))[1]
    if (!is.na(rank_by)) ov <- ov[order(-ov[[rank_by]]), , drop = FALSE]
    ov$label <- if ("label" %in% names(ov)) as.character(ov$label) else .gsm_label(ov$tSample)
    ov$rank <- seq_len(nrow(ov))
    k <- min(nrow(ov), n_overlap)
    shown <- ov[seq_len(k), , drop = FALSE]
    cols <- intersect(c("N_OL", setdiff(.overlap_rank_columns, "N_OL"), "p.adjust"), names(shown))
    ## counts are printed whole, ratios and p-values to 3 significant digits
    fmt <- lapply(cols, function(cl) {
        v <- if (cl == "N_OL") format(shown[[cl]], scientific = FALSE, trim = TRUE) else signif(shown[[cl]], 3)
        paste0(cl, "=", v)
    })
    vals <- if (length(cols)) do.call(paste, c(fmt, sep = ", ")) else rep("", k)
    lines <- paste0(shown$rank, "\t", shown$tSample, "\t", shown$label, "\t", vals)
    prm <- attr(overlap, "params")
    null <- if (is.null(prm$null) || is.na(prm$null)) {
        if (isTRUE(prm$pool)) "pooled over reference sets, shuffled over whole chromosomes" else NA_character_
    } else prm$null
    how <- if (is.null(prm)) "" else paste0(
        "p-values from ", prm$nShuffle, " shuffles; null: ", if (is.na(null)) "not recorded" else null,
        "; the smallest possible p is ", signif(1 / (prm$nShuffle + 1), 3), "\n",
        if (!is.na(null) && grepl("whole chromosomes", null))
            paste("A whole-genome null is anti-conservative for peaks in open chromatin: most reference",
                  "sets pass p.adjust < 0.05, so rely on the ranking, not on the p-values.\n") else "")
    ranking <- if (!is.na(rank_by) && rank_by == "jaccard")
        "jaccard (fold_enrichment is reported as an effect size, not used for ranking)" else
        if (is.na(rank_by)) "input order" else rank_by
    text <- paste0(
        "## Overlap with reference samples (", k, " of ", nrow(ov), " with p.adjust < 0.05, ",
        total, " tested; ranked by ", ranking, ")\n",
        how,
        if (k) paste(lines, collapse = "\n") else "none"
    )
    list(text = text, table = ov)
}

## tSample of downloaded GEO files ("GSM1295076_xxx.bed.gz") -> "<title> (<organism>)"
.gsm_label <- function(tsample) {
    gsm <- regmatches(tsample, regexpr("GSM[0-9]+", tsample))
    lab <- as.character(tsample)
    has <- grepl("GSM[0-9]+", tsample)
    if (!any(has)) return(lab)
    info <- tryCatch(get_gsminfo(), error = function(e) NULL)
    if (is.null(info)) return(lab)
    m <- match(gsm, info$gsm)
    title <- ifelse(is.na(m), NA_character_, paste0(info$title[m], " (", info$organism[m], ")"))
    lab[has] <- ifelse(is.na(title), lab[has], title)
    lab
}

.digest_distribution <- function(x) {
    st <- x@annoStat
    paste0("## Genomic distribution (% of regions)\n",
           paste0(st$Feature, ": ", round(st$Frequency, 1), collapse = "\n"))
}

.digest_tss <- function(d) {
    side <- ifelse(d < 0, "upstream", "downstream")
    bin <- cut(abs(d), .tss_breaks, right = FALSE,
               labels = c("0-1kb", "1-3kb", "3-5kb", "5-10kb", "10-100kb", ">100kb"))
    tab <- table(side, bin)
    pct <- round(100 * tab / length(d), 1)
    rows <- vapply(colnames(pct), function(b) {
        paste0(b, ": upstream ", if ("upstream" %in% rownames(pct)) pct["upstream", b] else 0,
               ", downstream ", if ("downstream" %in% rownames(pct)) pct["downstream", b] else 0)
    }, "")
    paste0("## Distance to the nearest TSS (% of regions)\n", paste(rows, collapse = "\n"))
}

.digest_width <- function(w) {
    q <- stats::quantile(w, c(0.25, 0.5, 0.75), names = FALSE)
    paste0("## Region width (bp)\n",
           "median ", q[2], ", IQR ", q[1], "-", q[3], ", min ", min(w), ", max ", max(w), "\n",
           "wider than 1 kb: ", round(100 * mean(w > 1e3), 1), "%, wider than 5 kb: ",
           round(100 * mean(w > 5e3), 1), "%")
}

.sha256 <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
    digest::digest(paste(as.character(x), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

## ---- csAnnoList: groups side by side ----------------------------------------------------

#' @noRd
.seq_digest_list <- function(x, context = NULL, assay = NULL, question = NULL, n_regions = 30,
                             signal = NULL, max_genes = 50) {
    if (!is(x, "csAnnoList")) stop("`x` must be a csAnnoList.")
    grp <- names(x)
    per <- lapply(x, function(a) {
        df <- as.data.frame(a@anno)
        if (!nrow(df)) stop("a group of `x` contains no regions.")
        ord <- .rank_regions(df, signal)
        df <- df[ord$index, , drop = FALSE]
        gene <- .region_gene(df)
        feature <- .short_annotation(df$annotation)
        list(a = a, df = df, gene = gene, feature = feature, rank_by = ord$column,
             prom = unique(gene[grepl("^Promoter", feature) & !is.na(gene)]))
    })
    sections <- list()
    sections$overview <- paste0(
        "## Overview\n",
        "Groups: ", paste0(grp, " (", vapply(per, function(p) nrow(p$df), 0), " regions)", collapse = ", "), "\n",
        "Genome: ", paste(unique(vapply(x, .source_line, "")), collapse = " | "), "\n",
        paste(vapply(grp, function(g) .checks_line(x[[g]], prefix = paste0("Checks (", g, "): ")), ""),
              collapse = ""),
        "Declared assay: ", if (is.null(assay)) "not given" else assay, "\n",
        "Context: ", if (is.null(context)) "not given" else context, "\n",
        "Question: ", if (is.null(question)) "none" else question
    )
    sections$distribution <- .digest_distribution_list(x)
    sections$tss_distance <- paste0(
        "## Distance to the nearest TSS (% of regions)\ngroup\twithin 1kb\t1-3kb\t3-10kb\t10-100kb\t>100kb\n",
        paste(vapply(grp, function(g) {
            d <- abs(per[[g]]$df$distanceToTSS)
            b <- table(cut(d, c(0, 1e3, 3e3, 1e4, 1e5, Inf), right = FALSE))
            paste(c(g, round(100 * as.numeric(b) / length(d), 1)), collapse = "\t")
        }, ""), collapse = "\n"))
    sections$width <- paste0(
        "## Region width (bp)\ngroup\tmedian\tIQR\twider than 1 kb (%)\n",
        paste(vapply(grp, function(g) {
            w <- per[[g]]$df$width
            q <- stats::quantile(w, c(0.25, 0.5, 0.75), names = FALSE)
            paste0(g, "\t", q[2], "\t", q[1], "-", q[3], "\t", round(100 * mean(w > 1e3), 1))
        }, ""), collapse = "\n"))
    shared <- Reduce(intersect, lapply(per, `[[`, "prom"))
    specific <- lapply(grp, function(g) setdiff(per[[g]]$prom, unlist(lapply(per[setdiff(grp, g)], `[[`, "prom"))))
    names(specific) <- grp
    sections$genes <- paste0(
        "## Promoter genes\n",
        "Shared by all groups: ", length(shared), if (length(shared)) paste0(" (first ", min(max_genes, length(shared)), ": ",
            paste(utils::head(shared, max_genes), collapse = ", "), ")"), "\n",
        paste(vapply(grp, function(g) {
            sp <- specific[[g]]
            paste0("Only in ", g, ": ", length(sp), " of ", length(per[[g]]$prom),
                   if (length(sp)) paste0(" (first ", min(max_genes, length(sp)), " in region rank order: ",
                                          paste(utils::head(sp, max_genes), collapse = ", "), ")"))
        }, ""), collapse = "\n"))
    sections$top_regions <- paste(vapply(grp, function(g) {
        p <- per[[g]]
        sub("^## Top regions", paste0("## Top regions of ", g),
            .digest_regions(p$df, p$gene, p$feature, p$rank_by, n_regions))
    }, ""), collapse = "\n\n")
    structure(list(
        text = paste(unlist(sections), collapse = "\n\n"),
        sections = sections,
        sha256 = vapply(sections, .sha256, ""),
        evidence = list(shared_promoter_genes = shared, specific_promoter_genes = specific,
                        promoter_genes = lapply(per, `[[`, "prom"),
                        all_genes = lapply(per, function(p) unique(p$gene[!is.na(p$gene)]))),
        totals = list(groups = grp, n_regions = vapply(per, function(p) nrow(p$df), 0),
                      n_regions_sent = vapply(per, function(p) min(nrow(p$df), n_regions), 0))
    ), class = "seq_digest")
}

.digest_distribution_list <- function(x) {
    feats <- unique(unlist(lapply(x, function(a) as.character(a@annoStat$Feature))))
    rows <- vapply(feats, function(f) {
        v <- vapply(x, function(a) {
            st <- a@annoStat
            i <- match(f, as.character(st$Feature))
            if (is.na(i)) 0 else round(st$Frequency[i], 1)
        }, 0)
        paste(c(f, v), collapse = "\t")
    }, "")
    paste0("## Genomic distribution (% of regions)\nfeature\t", paste(names(x), collapse = "\t"), "\n",
           paste(rows, collapse = "\n"))
}

#' @method print seq_digest
#' @export
print.seq_digest <- function(x, ...) {
    cat(x$text, "\n")
    invisible(x)
}
