## interpret() for the base group (bmData) and the pair group (loopAnno). Both use the shared
## skeleton of R/interpret.R and only add an evidence digest and a schema. The file name puts
## it after bmData-class.R and loopAnno.R, whose classes the methods are for (no Collate).
## Design: checks/20261107_interpret-base-pair/DESIGN.md in the epiSeeker-LLM records.

#' Interpret single-base data and chromatin loops with a language model (experimental)
#'
#' Methods of [clusterProfiler::interpret()] for the base group ([bmData-class]) and the pair
#' group ([loopAnno-class]). Like the methods for annotated regions
#' ([interpret-epiSeeker]), each builds an evidence digest without a model and makes one
#' structured call; every digest section says how many items it sends and how many there are,
#' and is hashed into the execution record ([clusterProfiler::provenance()]). Both methods are
#' **experimental**: neither has been compared with a baseline that uses no model.
#'
#' * `bmData`: sites are merged into clusters (sites within `gap` bp); the digest lists the
#'   clusters with the mean level per sample (and, with `TxDb`, their feature and nearest gene
#'   from [annotateSeq()]), the differences between samples, the distribution of levels, the
#'   coverage, and the genomic distribution of the clusters.
#' * `loopAnno`: the genes at promoter anchors ranked by the number of loops on them, the
#'   peak links of [linkPeaksToGenes()] when `links` is given, the loop classes, the anchor
#'   distribution and the span of each class, and a fixed caveat: loop partners are candidate
#'   genes (in the package's CRISPRi benchmark evaluation real loops were no more precise than
#'   distance-matched rewired loops, 0.220 vs 0.219).
#'
#' @param x a `bmData` or a `loopAnno`
#' @param context,question,model,max_tokens,temperature,reasoning as in [interpret-epiSeeker]
#' @param assay declared assay, e.g. "WGBS", "oxBS-seq" (for a `bmData`)
#' @param value,coverage names of the assays holding the level (in \[0, 1\]) and the coverage.
#'   `NULL` takes the first assay whose values all lie in \[0, 1\] as the level and the first
#'   other assay of non-negative whole numbers as the coverage.
#' @param TxDb optional `TxDb` to annotate the clusters of a `bmData`
#' @param annoDb optional annotation package, adds gene symbols
#' @param gap sites closer than this (bp) form one cluster. Default 500.
#' @param min_cov sites with less coverage are left out of levels and differences. Default 5.
#' @param n_regions number of clusters sent. Default 30.
#' @param n_sites number of sites with the largest difference sent per pair of samples. Default 20.
#' @param links optional data frame from [linkPeaksToGenes()] for the same loops
#' @param n_genes number of genes sent in each gene list of a `loopAnno` digest. Default 50.
#' @param ... not used
#' @return a `seq_interpretation` list. For a `bmData`: `answer`, `candidates`,
#'   `methylation_profile`, `sample_differences`, `genomic_context`, `data_quality`,
#'   `narrative`. For a `loopAnno`: `answer`, `candidates`, `loop_profile`, `gene_hubs`,
#'   `peak_links`, `narrative`. The digest is in `attr(, "evidence")`, the execution record
#'   in `attr(, "provenance")`.
#' @name interpret-base-pair
#' @aliases interpret,bmData-method interpret,loopAnno-method
#' @examples
#' data(demo_bmdata)
#' ## the evidence the model would see, without calling one
#' cat(epiSeeker:::.bm_digest(demo_bmdata)$text)
#' if (requireNamespace("aisdk", quietly = TRUE) && nzchar(Sys.getenv("OPENAI_API_KEY"))) {
#'     res <- interpret(demo_bmdata, assay = "WGBS", model = "openai:gpt-4o-mini")
#' }
NULL

#' @rdname interpret-base-pair
setMethod("interpret", "bmData", function(x, context = NULL, assay = NULL, question = NULL,
                                          value = NULL, coverage = NULL, TxDb = NULL,
                                          annoDb = NULL, gap = 500, min_cov = 5,
                                          n_regions = 30, n_sites = 20, model = NULL,
                                          max_tokens = 4096, temperature = NULL,
                                          reasoning = NULL, ...) {
    ev <- .bm_digest(x, context = context, assay = assay, question = question, value = value,
                     coverage = coverage, TxDb = TxDb, annoDb = annoDb, gap = gap,
                     min_cov = min_cov, n_regions = n_regions, n_sites = n_sites)
    .seq_interpret_run(ev, schema = .bm_schema(), entry_point = "interpret,bmData",
                       model = model, max_tokens = max_tokens, temperature = temperature,
                       reasoning = reasoning,
                       system = .seq_system_prompt("methylation levels at single bases"))
})

#' @rdname interpret-base-pair
setMethod("interpret", "loopAnno", function(x, context = NULL, question = NULL, links = NULL,
                                            n_genes = 50, model = NULL, max_tokens = 4096,
                                            temperature = NULL, reasoning = NULL, ...) {
    ev <- .loop_digest(x, context = context, question = question, links = links,
                       n_genes = n_genes)
    .seq_interpret_run(ev, schema = .loop_schema(), entry_point = "interpret,loopAnno",
                       model = model, max_tokens = max_tokens, temperature = temperature,
                       reasoning = reasoning,
                       system = .seq_system_prompt("chromatin loops between pairs of anchors"))
})

## ---- bmData digest ------------------------------------------------------------------------

## which assay is the level and which the coverage
.bm_assays <- function(x, value, coverage) {
    an <- SummarizedExperiment::assayNames(x)
    if (is.null(an) || !length(an)) stop("the bmData has no named assays.")
    vals <- function(nm) as.matrix(SummarizedExperiment::assay(x, nm))
    if (is.null(value)) {
        is_level <- vapply(an, function(nm) {
            v <- vals(nm); v <- v[is.finite(v)]
            length(v) > 0 && all(v >= 0 & v <= 1)
        }, logical(1))
        if (!any(is_level)) {
            stop("no assay of the bmData holds levels in [0, 1]; name the methylation-level ",
                 "assay with `value =` (assays: ", paste(an, collapse = ", "), ").")
        }
        value <- an[is_level][1]
    } else if (!value %in% an) {
        stop("`value` '", value, "' is not an assay of the bmData (", paste(an, collapse = ", "), ").")
    }
    if (is.null(coverage)) {
        is_cov <- vapply(setdiff(an, value), function(nm) {
            v <- vals(nm); v <- v[is.finite(v)]
            length(v) > 0 && all(v >= 0 & v == round(v))
        }, logical(1))
        coverage <- if (any(is_cov)) names(is_cov)[is_cov][1] else NA_character_
    } else if (!coverage %in% an) {
        stop("`coverage` '", coverage, "' is not an assay of the bmData (", paste(an, collapse = ", "), ").")
    }
    list(value = value, coverage = coverage)
}

#' @noRd
.bm_digest <- function(x, context = NULL, assay = NULL, question = NULL, value = NULL,
                       coverage = NULL, TxDb = NULL, annoDb = NULL, gap = 500, min_cov = 5,
                       n_regions = 30, n_sites = 20) {
    if (!is(x, "bmData")) stop("`x` must be a bmData.")
    gr <- SummarizedExperiment::rowRanges(x)
    if (!length(gr)) stop("`x` contains no sites.")
    an <- .bm_assays(x, value, coverage)
    lev <- as.matrix(SummarizedExperiment::assay(x, an$value))
    cov <- if (is.na(an$coverage)) NULL else as.matrix(SummarizedExperiment::assay(x, an$coverage))
    smp <- colnames(x)
    if (is.null(smp)) smp <- paste0("sample", seq_len(ncol(x)))
    ## sites below min_cov do not count for levels and differences
    ok <- if (is.null(cov)) matrix(TRUE, nrow(lev), ncol(lev)) else !is.na(cov) & cov >= min_cov
    lv <- lev
    lv[!ok | !is.finite(lv)] <- NA

    ## clusters of nearby sites
    cl <- GenomicRanges::reduce(gr, min.gapwidth = gap + 1L, ignore.strand = TRUE, with.revmap = TRUE)
    rev <- cl$revmap
    n_sites_cl <- lengths(rev)
    cl_mean <- t(vapply(as.list(rev), function(i) colMeans(lv[i, , drop = FALSE], na.rm = TRUE),
                        numeric(ncol(lv))))
    if (ncol(lv) == 1) cl_mean <- matrix(cl_mean, ncol = 1)
    colnames(cl_mean) <- paste0("mean_level_", smp)
    feature <- gene <- rep(NA_character_, length(cl))
    anno <- NULL
    if (!is.null(TxDb)) {
        cl_plain <- GenomicRanges::granges(cl)
        anno <- tryCatch(annotateSeq(cl_plain, TxDb = TxDb, annoDb = annoDb, verbose = FALSE),
                         error = function(e) NULL)
        if (!is.null(anno)) {
            a <- anno@anno
            hit <- match(paste(GenomicRanges::seqnames(a), GenomicRanges::start(a)),
                         paste(GenomicRanges::seqnames(cl), GenomicRanges::start(cl)))
            ## a self-defined GRanges TxDb gives no annotation column and a gene_id column
            if (!is.null(a$annotation)) feature[hit] <- .short_annotation(a$annotation)
            g <- if (!is.null(a$SYMBOL)) a$SYMBOL else if (!is.null(a$geneId)) a$geneId else a$gene_id
            if (!is.null(g)) gene[hit] <- as.character(g)
        }
    }

    sections <- list()
    chr_tab <- table(as.character(GenomicRanges::seqnames(gr)))
    rng <- vapply(names(chr_tab), function(ch) {
        s <- GenomicRanges::start(gr)[as.character(GenomicRanges::seqnames(gr)) == ch]
        paste0(ch, ":", min(s), "-", max(s), " (", chr_tab[[ch]], " sites)")
    }, character(1))
    g <- unique(stats::na.omit(GenomeInfoDb::genome(gr)))
    sections$overview <- paste0(
        "## Overview\n",
        "Sites: ", length(gr), " single bases in ", length(cl), " clusters (sites within ", gap,
        " bp); samples: ", paste(smp, collapse = ", "), "\n",
        "Covered: ", paste(rng, collapse = "; "), "\n",
        "Genome: ", if (length(g)) paste(g, collapse = ",") else "not recorded", "\n",
        "Level assay: ", an$value, "; coverage assay: ",
        if (is.na(an$coverage)) "none" else paste0(an$coverage, " (sites with coverage < ", min_cov,
                                                   " are left out of levels and differences)"), "\n",
        "Declared assay: ", if (is.null(assay)) "not given" else assay, "\n",
        "Context: ", if (is.null(context)) "not given" else context, "\n",
        "Question: ", if (is.null(question)) "none" else question
    )

    ord <- order(-n_sites_cl, as.character(GenomicRanges::seqnames(cl)), GenomicRanges::start(cl))
    k <- min(length(cl), n_regions)
    o <- ord[seq_len(k)]
    fmt <- function(v) ifelse(is.na(v), "NA", formatC(v, digits = 2, format = "f"))
    lines <- vapply(o, function(i) {
        paste0(GenomicRanges::seqnames(cl)[i], ":", GenomicRanges::start(cl)[i], "-",
               GenomicRanges::end(cl)[i], "\t", n_sites_cl[i], "\t",
               paste(fmt(cl_mean[i, ]), collapse = "\t"),
               if (!is.null(anno)) paste0("\t", feature[i], "\t", gene[i]) else "")
    }, character(1))
    sections$regions <- paste0(
        "## Clusters (", k, " of ", length(cl), ", most sites first)\n",
        "cluster\tsites\t", paste(paste0("mean_level_", smp), collapse = "\t"),
        if (!is.null(anno)) "\tfeature\tnearest_gene" else "", "\n",
        paste(lines, collapse = "\n"))

    if (ncol(lv) >= 2) {
        pairs <- utils::combn(seq_len(ncol(lv)), 2, simplify = FALSE)
        txt <- vapply(pairs, function(p) {
            a <- lv[, p[1]]; b <- lv[, p[2]]
            both <- !is.na(a) & !is.na(b)
            d <- (b - a)[both]
            big <- sum(abs(d) >= 0.2)
            idx <- which(both)[order(-abs(d))][seq_len(min(n_sites, sum(both)))]
            top <- if (length(idx)) paste0(
                GenomicRanges::seqnames(gr)[idx], ":", GenomicRanges::start(gr)[idx], " ",
                fmt(a[idx]), " -> ", fmt(b[idx]), collapse = "\n") else "none"
            paste0("### ", smp[p[1]], " vs ", smp[p[2]], "\n",
                   "sites covered in both: ", sum(both), "; mean level ", fmt(mean(a[both])), " vs ",
                   fmt(mean(b[both])), "; mean difference ", fmt(mean(d)), "\n",
                   "sites with |difference| >= 0.2: ", big, " of ", sum(both), "\n",
                   "largest differences (", length(idx), " of ", sum(both), " sites, ",
                   smp[p[1]], " -> ", smp[p[2]], "):\n", top)
        }, character(1))
        sections$differences <- paste0("## Differences between samples\n", paste(txt, collapse = "\n"))
    } else {
        sections$differences <- "## Differences between samples\none sample; no comparison"
    }

    lv_txt <- vapply(seq_len(ncol(lv)), function(j) {
        v <- lv[, j]; v <- v[!is.na(v)]
        if (!length(v)) return(paste0(smp[j], ": no site with a value"))
        paste0(smp[j], ": median ", fmt(stats::median(v)), "; < 0.2: ",
               round(100 * mean(v < 0.2), 1), "%; 0.2-0.8: ",
               round(100 * mean(v >= 0.2 & v <= 0.8), 1), "%; > 0.8: ",
               round(100 * mean(v > 0.8), 1), "% (", length(v), " sites; ",
               nrow(lv) - length(v), " without a value)")
    }, character(1))
    sections$levels <- paste0("## Methylation levels per sample\n", paste(lv_txt, collapse = "\n"))

    sections$coverage <- if (is.null(cov)) "## Coverage\nno coverage assay" else {
        paste0("## Coverage per sample\n", paste(vapply(seq_len(ncol(cov)), function(j) {
            v <- cov[, j]; v <- v[!is.na(v)]
            q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
            paste0(smp[j], ": median ", q[2], ", IQR ", q[1], "-", q[3], "; below ", min_cov, ": ",
                   round(100 * mean(v < min_cov), 1), "%")
        }, character(1)), collapse = "\n"))
    }

    sections$genomic_context <- if (is.null(anno)) {
        "## Genomic context of the clusters\nnot annotated (no TxDb given)"
    } else {
        st <- anno@annoStat
        paste0("## Genomic context of the clusters (% of ", length(anno@anno), " clusters)\n",
               "Genome: ", .source_line(anno), "\n",
               ## a self-defined GRanges TxDb gives nearest genes but no genomic features
               if (nrow(st)) paste0(st$Feature, ": ", round(st$Frequency, 1), collapse = "\n")
               else "no genomic features (self-defined TxDb); nearest genes are in the clusters section")
    }

    structure(list(
        text = paste(unlist(sections), collapse = "\n\n"),
        sections = sections,
        sha256 = vapply(sections, .sha256, ""),
        evidence = list(clusters = data.frame(
            seqnames = as.character(GenomicRanges::seqnames(cl)), start = GenomicRanges::start(cl),
            end = GenomicRanges::end(cl), sites = n_sites_cl, cl_mean, feature = feature, gene = gene,
            stringsAsFactors = FALSE, check.names = FALSE),
            value = an$value, coverage = an$coverage, samples = smp),
        totals = list(n_sites = length(gr), n_clusters = length(cl), n_clusters_sent = k,
                      n_samples = ncol(lv))
    ), class = "seq_digest")
}

.bm_schema <- function() {
    aisdk::z_object(
        answer = .answer_field(),
        candidates = .candidates_field(),
        methylation_profile = aisdk::z_string(paste(
            "What the methylation looks like (high, low, intermediate; where), citing the levels",
            "and percentages given.")),
        sample_differences = aisdk::z_string(paste(
            "How the samples differ, citing the sites and differences given; 'None' with one sample.")),
        genomic_context = aisdk::z_string(paste(
            "Which features and genes the clusters are in, as given; 'Not annotated' when no",
            "annotation was given.")),
        data_quality = aisdk::z_string(
            "Coverage and missing values, citing the numbers given, and what they limit."),
        narrative = aisdk::z_string("A short paragraph that ties the evidence together.")
    )
}

## ---- loopAnno digest ----------------------------------------------------------------------

#' @noRd
.loop_digest <- function(x, context = NULL, question = NULL, links = NULL, n_genes = 50) {
    if (!is(x, "loopAnno")) stop("`x` must be a loopAnno from annotateLoops().")
    gi <- x@loops
    mc <- S4Vectors::mcols(gi)
    n <- length(gi)
    if (!n) stop("`x` contains no loops.")
    src <- S4Vectors::metadata(gi)$loopSource
    span <- mc$span
    cis <- !is.na(span)
    q <- if (any(cis)) stats::quantile(span[cis], c(0.25, 0.5, 0.75), names = FALSE) else rep(NA, 3)

    sections <- list()
    sections$overview <- paste0(
        "## Overview\n",
        "Loops: ", n, " between ", length(InteractionSet::regions(gi)), " anchors",
        if (!is.null(src)) paste0(" (from ", src$source, "; ", src$dropped, " rows without coordinates dropped)"),
        "\n",
        "Genome: ", .source_line(x@anchors), "\n",
        "Promoter anchor: an anchor overlapping the window ", x@tssRegion[1], " to ", x@tssRegion[2],
        " bp around a TSS",
        if (!is.null(S4Vectors::metadata(gi)$geneFilter))
            paste0(" of one of ", S4Vectors::metadata(gi)$geneFilter$n_gene_ids,
                   " genes given (gene_ids); other genes are not counted"),
        "\n",
        "Loops between chromosomes: ", sum(!cis), "\n",
        "Span (bp, loops within one chromosome): median ", q[2], ", IQR ", q[1], "-", q[3], "\n",
        "Context: ", if (is.null(context)) "not given" else context, "\n",
        "Question: ", if (is.null(question)) "none" else question
    )

    ## genes at promoter anchors: number of loops, and what is at the other end
    has_sym <- !is.null(mc$anchor1_symbols)
    ## symbols when the loops carry them; a gene without a symbol keeps its id
    pick <- function(sym, ids) {
        ids <- lapply(as.list(ids), function(z) unique(as.character(z)))
        if (!has_sym) return(ids)
        sym <- lapply(as.list(sym), function(z) unique(as.character(z)))
        Map(function(s, i) if (length(s)) s else i, sym, ids)
    }
    g1 <- pick(mc$anchor1_symbols, mc$anchor1_genes)
    g2 <- pick(mc$anchor2_symbols, mc$anchor2_genes)
    gene <- c(unlist(g1), unlist(g2))
    other_prom <- c(rep(mc$anchor2_promoter, lengths(g1)), rep(mc$anchor1_promoter, lengths(g2)))
    if (length(gene)) {
        n_loop <- table(gene)
        to_prom <- tapply(other_prom, gene, sum)[names(n_loop)]
        tab <- data.frame(gene = names(n_loop), loops = as.integer(n_loop),
                          to_promoter = as.integer(to_prom), stringsAsFactors = FALSE)
        tab$to_distal <- tab$loops - tab$to_promoter
        tab <- tab[order(-tab$loops, tab$gene), , drop = FALSE]
        rownames(tab) <- NULL
    } else {
        tab <- data.frame(gene = character(0), loops = integer(0), to_promoter = integer(0),
                          to_distal = integer(0))
    }
    k <- min(nrow(tab), n_genes)
    sections$gene_hubs <- paste0(
        "## Genes at promoter anchors (", k, " of ", nrow(tab), ", most loops first)\n",
        "gene\tloops\tother end promoter\tother end distal\n",
        if (k) paste(paste(tab$gene[seq_len(k)], tab$loops[seq_len(k)], tab$to_promoter[seq_len(k)],
                           tab$to_distal[seq_len(k)], sep = "\t"), collapse = "\n") else "none")

    sections$peak_links <- if (is.null(links)) "## Peak links\nno peaks given" else .loop_links_section(links, n_genes)

    st <- x@loopStat
    sections$loop_classes <- paste0("## Loop classes (% of ", n, " loops)\n",
                                    paste0(st$loop_class, ": ", st$n, " (", round(st$Frequency, 1), "%)",
                                           collapse = "\n"))
    ast <- x@anchors@annoStat
    sections$anchor_distribution <- paste0(
        "## Genomic distribution of the anchors (% of ", length(x@anchors@anno), " anchors)\n",
        if (nrow(ast)) paste0(ast$Feature, ": ", round(ast$Frequency, 1), collapse = "\n") else "not available")
    sp <- vapply(levels(mc$loop_class), function(cl) {
        s <- span[cis & mc$loop_class == cl]
        if (!length(s)) return(paste0(cl, ": no loop within one chromosome"))
        qq <- stats::quantile(s, c(0.25, 0.5, 0.75), names = FALSE)
        paste0(cl, ": ", length(s), " loops, median ", qq[2], ", IQR ", qq[1], "-", qq[3])
    }, character(1))
    sections$span <- paste0("## Span per class (bp)\n", paste(sp, collapse = "\n"))
    sections$caveat <- .loop_caveat


    structure(list(
        text = paste(unlist(sections), collapse = "\n\n"),
        sections = sections,
        sha256 = vapply(sections, .sha256, ""),
        evidence = list(gene_hubs = tab),
        totals = list(n_loops = n, n_anchors = length(InteractionSet::regions(gi)),
                      n_genes = nrow(tab), n_genes_sent = k)
    ), class = "seq_digest")
}

## A fact the digest always carries (checks/20261109_fig4-distmatch in the epiSeeker-LLM
## records): what a loop partner is worth as evidence of regulation.
.loop_caveat <- paste0(
    "## Caveat on loop partners\n",
    "Genes at the other end of a loop are candidate genes, not established targets. In the ",
    "package's evaluation on a CRISPRi enhancer-gene benchmark (K562, 9,830 tested pairs), genes ",
    "linked to elements by real loops were regulated no more often than genes linked by ",
    "distance-matched rewired loops (precision 0.220 vs 0.219; enrichment 1.00, p = 0.47): loops ",
    "added recall over the nearest gene, not specificity.")

.loop_links_section <- function(links, n_genes) {
    if (!is.data.frame(links) || is.null(links$link)) {
        stop("`links` must be the data frame returned by linkPeaksToGenes().")
    }
    l <- as.character(links$link)
    d <- l != "promoter"
    tab <- table(factor(l[d], levels = setdiff(.link_classes, "promoter")))
    cnt <- paste0(names(tab), ": ", as.integer(tab), " (", round(100 * as.integer(tab) / max(1, sum(d)), 1), "%)",
                  collapse = "\n")
    other <- which(l == "loop, other genes")
    gcol <- if (!is.null(links$loop_symbols)) "loop_symbols" else "loop_genes"
    ncol <- if (!is.null(links$nearest_symbol)) "nearest_symbol" else "nearest_gene"
    via <- character(0)
    if (length(other)) {
        gl <- strsplit(as.character(links[[gcol]][other]), "/", fixed = TRUE)
        near <- as.character(links[[ncol]][other])
        df <- data.frame(gene = unlist(gl), nearest = rep(near, lengths(gl)), stringsAsFactors = FALSE)
        df <- df[nzchar(df$gene), , drop = FALSE]
        if (nrow(df)) {
            agg <- split(df$nearest, df$gene)
            npk <- vapply(agg, length, integer(1))
            o <- order(-npk, names(npk))
            k <- min(length(o), n_genes)
            via <- vapply(o[seq_len(k)], function(i) {
                nn <- sort(table(agg[[i]]), decreasing = TRUE)
                paste0(names(npk)[i], "\t", npk[i], "\t",
                       paste0(names(nn)[seq_len(min(3, length(nn)))], collapse = "/"))
            }, character(1))
            via <- c(paste0("Genes reached only through loops (", k, " of ", length(o),
                            ", most peaks first):"),
                     "gene\tpeaks\tnearest gene of those peaks (up to 3)", via)
        }
    }
    paste0("## Peak links (", sum(d), " of ", length(l), " peaks are outside promoters)\n", cnt,
           if (length(via)) paste0("\n", paste(via, collapse = "\n")) else "")
}

.loop_schema <- function() {
    aisdk::z_object(
        answer = .answer_field(),
        candidates = .candidates_field(),
        loop_profile = aisdk::z_string(paste(
            "What the loops look like (classes, spans, where the anchors are), citing the numbers given.")),
        gene_hubs = aisdk::z_string(paste(
            "Which genes sit on the most loops and what they connect to, naming genes as given.")),
        peak_links = aisdk::z_string(paste(
            "What the loops add to nearest-gene assignment of peaks, citing the counts and genes",
            "given; 'No peaks given' when there are none.")),
        narrative = aisdk::z_string("A short paragraph that ties the evidence together.")
    )
}
