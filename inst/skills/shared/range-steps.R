## Steps shared by the skills for bulk region sets (range-tf, range-histone, range-atac). Like
## skill-helpers.R (which must be sourced first) they hold no analysis of their own: they order
## the epiSeeker calls the same way in every skill and record each check. A skill script
## sources both files with system.file(), so a copied skill folder still runs.

## file names -> sample names (without the usual peak-file extensions)
range_clean_names <- function(files) sub("\\.(bed|narrowPeak|broadPeak|txt)(\\.gz)?$", "", basename(files))

## ---- 1. inputs ------------------------------------------------------------------------------

## --peaks (+ --compare), --names, --txdb, --genome: read every region set and run the input
## checks (TxDb named, genome = TxDb, genome recorded on the regions, seqlevels style).
range_inputs <- function(a) {
    if (is.null(a$peaks)) skill_check("input", "peaks given", "STOP", "give --peaks <file>")
    if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP",
                                     "give --txdb <TxDb package>; annotateSeq() would otherwise fall back to hg38")
    txdb <- skill_load(a$txdb)
    skill_check_genome(a$genome, txdb)
    files <- c(a$peaks, skill_split(a$compare))
    nm <- skill_split(a$names)
    if (is.null(nm)) nm <- range_clean_names(files)
    if (length(nm) != length(files) || anyDuplicated(nm)) {
        skill_check("input", "one unique name per peak file", "STOP", "--names must name --peaks and --compare in order")
    }
    skill_input(files, "peaks")
    peaks <- stats::setNames(lapply(files, readPeakFile), nm)
    for (s in nm) {
        skill_check_recorded_genome(peaks[[s]], a$genome, what = s)
        peaks[[s]] <- skill_check_seqlevels(peaks[[s]], txdb, fix = a$fix_seqlevels, what = s)
    }
    w <- vapply(peaks, function(p) stats::median(GenomicRanges::width(p)), 0)
    list(txdb = txdb, peaks = peaks, names = nm, width = w)
}

## ---- 2. annotation ----------------------------------------------------------------------------

## annotateSeq() of every set; annoChecks() of each before anything is drawn from it; then the
## query's annotation table, the feature table of all sets and the two standard figures.
range_annotate <- function(regions, txdb, annodb, tss, allow_suspect = FALSE, what = "") {
    nm <- names(regions)
    anno <- lapply(nm, function(s) skill_capture(
        annotateSeq(regions[[s]], TxDb = txdb, annoDb = annodb, tssRegion = c(-tss, tss), verbose = FALSE),
        paste("annotateSeq", s)))
    names(anno) <- nm
    for (s in nm) skill_check_anno(anno[[s]], paste0(what, s), allow_suspect = allow_suspect)
    skill_table(as.data.frame(anno[[1]]), paste0("annotation_", nm[1], ".tsv"))
    stat <- do.call(rbind, lapply(nm, function(s) cbind(set = s, getAnnoStat(anno[[s]]))))
    skill_table(stat, "annotation_stat.tsv")
    skill_figure(plotAnnoBar(csAnnoList(anno)), "annotation_bar.pdf")
    skill_figure(plotDistToTSS(csAnnoList(anno)), "distance_to_tss.pdf")
    anno
}

## ---- assay label check (seqAssayQC) ------------------------------------------------------------

## seqAssayQC() of every peak set against ENCODE peak sets of the declared --assay: a flag (WARN)
## when the promoter fraction, median width or gene-body fraction is outside the range of that
## assay. The check is partial: on held-out ENCODE files it flagged 3.6% of correctly labelled
## sets and 59% of narrow <-> broad swaps, and no broad mark declared as "TF". So no flag is a
## NOTE ("not contradicted"), never a PASS: it does not verify the assay. The ranges come from
## peak files annotated with annotateSeq()'s defaults, so the peaks are annotated that way for the
## check, whatever --tss the analysis uses (the build and dropped-peak checks are those of the
## analysed annotation, recorded there). range_keep_assay_qc() then keeps each result in
## annoChecks(<annotation analysed>)$assay, so that interpret(assay =) reports the same result
## (for range-histone: the peaks', not the merged domains').
range_assay_qc <- function(peaks, txdb, assay) {
    nm <- names(peaks)
    partial <- paste("partial check: on held-out ENCODE files 59% of narrow/broad swaps flagged,",
                     "no broad mark declared as TF; no flag does not verify the assay")
    if (is.null(assay) || isTRUE(assay) || !nzchar(trimws(assay))) {
        skill_check("input", "assay label (seqAssayQC)", "NOTE", "no --assay given; not checked")
        return(NULL)
    }
    qc <- lapply(nm, function(s) {
        an <- suppressWarnings(annotateSeq(peaks[[s]], TxDb = txdb, verbose = FALSE))
        seqAssayQC(an, assay)
    })
    names(qc) <- nm
    for (s in nm) {
        q <- qc[[s]]
        met <- sprintf("promoter fraction (<= 1 kb) %.1f%%, median width %s bp, gene-body fraction %.1f%%",
                       100 * q$metrics[["promoter_1kb"]], format(round(q$metrics[["width_median"]]), big.mark = ","),
                       100 * q$metrics[["gene_body"]])
        fits <- if (length(q$fits)) paste0("fits ENCODE ", paste(q$fits, collapse = ", ")) else "fits no ENCODE class"
        if (isFALSE(q$consistent)) {
            skill_check("input", paste("assay label of", s, "(seqAssayQC)"), "WARN",
                        paste0("does not look like ", q$assay, ": ", paste(q$reasons, collapse = "; "), "; ", fits,
                               ". Check the assay label of this file"))
        } else if (isTRUE(q$consistent)) {
            skill_check("input", paste("assay label of", s, "(seqAssayQC)"), "NOTE",
                        paste0("not contradicted for ", q$class, " (", met, "); ", partial,
                               if (length(q$note)) paste0("; ", paste(q$note, collapse = "; ")) else ""))
        } else {
            skill_check("input", paste("assay label of", s, "(seqAssayQC)"), "NOTE",
                        paste0("not assessable: ", paste(q$reasons, collapse = "; ")))
        }
    }
    tab <- do.call(rbind, lapply(nm, function(s) {
        q <- qc[[s]]
        data.frame(set = s, assay = q$assay, class = if (is.na(q$class)) NA else q$class,
                   assessable = q$assessable, flagged = isFALSE(q$consistent),
                   n_peaks = q$metrics[["n_peaks"]], promoter_1kb = round(q$metrics[["promoter_1kb"]], 4),
                   width_median = q$metrics[["width_median"]], gene_body = round(q$metrics[["gene_body"]], 4),
                   fits = paste(q$fits, collapse = ","), reasons = paste(q$reasons, collapse = "; "))
    }))
    skill_table(tab, "assay_qc.tsv")
    qc
}

## keep the assay check with the annotation interpret() reads
range_keep_assay_qc <- function(anno, qc) {
    for (s in intersect(names(anno), names(qc))) S4Vectors::metadata(anno[[s]]@anno)$annoChecks$assay <- qc[[s]]
    anno
}

## the report lines of the assay check
range_assay_section <- function(qc) {
    if (is.null(qc)) return("not checked (no --assay)")
    c(vapply(names(qc), function(s) {
        q <- qc[[s]]
        paste0("- ", s, ": ", if (isFALSE(q$consistent)) paste0("**flagged** - does not look like ", q$assay, ": ",
                                                               paste(q$reasons, collapse = "; ")) else
            if (isTRUE(q$consistent)) paste0("no flag for ", q$class) else
                paste0("not assessable (", paste(q$reasons, collapse = "; "), ")"))
    }, ""),
    paste("- seqAssayQC() compares promoter fraction, median width and gene-body fraction with ENCODE peak sets;",
          "it is partial (held-out: 59% of narrow/broad swaps flagged, none of the broad marks declared as TF),",
          "so no flag does not verify the assay. tables/assay_qc.tsv"))
}

## % of each set within the promoter window, from getAnnoStat()
range_promoter_pct <- function(anno) {
    vapply(names(anno), function(s) {
        st <- getAnnoStat(anno[[s]])
        sum(st$Frequency[grepl("^Promoter", st$Feature)])
    }, 0)
}

## profile of the region sets around TSSs (or another window type) and its figure
range_profile <- function(regions, txdb, upstream, downstream, file, type = "start_site", nbin = NULL,
                          step = "getTagMatrix") {
    tm <- skill_capture(getTagMatrix(regions, upstream = upstream, downstream = downstream, type = type,
                                     by = "gene", TxDb = txdb, nbin = nbin, verbose = FALSE), step)
    skill_figure(plotPeakProf(tm, facet = "none", conf = NULL), file)
    invisible(tm)
}

## ---- overlap with reference sets in a background ----------------------------------------------

## enrichPeakOverlap(method = "fisher") of `query` against --targets in --background. The
## background is required; the test is limited to the chromosomes the background covers (units
## outside it would count as query-only); a background made of little more than the tested sets
## is flagged; the result is ranked by its effect size, never by p-value.
range_overlap <- function(query, a, txdb) {
    targets <- skill_split(a$targets)
    if (length(targets) == 1 && dir.exists(targets)) targets <- list.files(targets, full.names = TRUE)
    if (!length(targets)) return(NULL)
    skill_input(targets, "target")
    if (is.null(a$background) || isTRUE(a$background)) {
        skill_check("overlap", "background given", "STOP", paste(
            "give --background: open chromatin (DNase/ATAC) of the same cell type, or the union of a large",
            "reference collection. A whole-genome null calls almost every pair significant",
            "(FPR 0.996 in the package's evaluation)"))
    }
    skill_input(a$background, "background")
    bg <- readPeakFile(a$background)
    skill_note("background", a$background)
    bg <- skill_check_seqlevels(bg, txdb, fix = a$fix_seqlevels, what = "background")
    tg <- stats::setNames(lapply(targets, readPeakFile), range_clean_names(targets))
    for (s in names(tg)) tg[[s]] <- skill_check_seqlevels(tg[[s]], txdb, fix = a$fix_seqlevels, what = paste("target", s))

    ## limit the test to the chromosomes of the background
    bchr <- unique(as.character(GenomicRanges::seqnames(bg)))
    inb <- function(g) g[as.character(GenomicRanges::seqnames(g)) %in% bchr]
    q2 <- inb(GenomicRanges::granges(query))
    tg2 <- lapply(tg, function(g) inb(GenomicRanges::granges(g)))
    n_out <- (length(query) - length(q2)) + (sum(lengths(tg)) - sum(lengths(tg2)))
    if (n_out > 0) {
        skill_check("overlap", "test limited to the chromosomes of the background", "NOTE",
                    sprintf("background covers %s; the test uses %d of %d query regions and %d of %d target regions on them",
                            paste(utils::head(bchr, 5), collapse = ", "), length(q2), length(query),
                            sum(lengths(tg2)), sum(lengths(tg))))
    }
    if (!length(q2)) skill_check("overlap", "query regions on the chromosomes of the background", "STOP", "none")
    empty <- lengths(tg2) == 0
    if (any(empty)) {
        skill_check("overlap", "targets with regions on the chromosomes of the background", "NOTE",
                    paste("left out:", paste(names(tg2)[empty], collapse = ", ")))
        tg2 <- tg2[!empty]
    }
    if (!length(tg2)) skill_check("overlap", "targets with regions on the chromosomes of the background", "STOP", "none")
    ## a universe made of little more than the tested sets: every unit is hit by something, so
    ## p-values are conservative (p near 1) -- rank, do not test, and say so
    own <- mean(IRanges::overlapsAny(GenomicRanges::reduce(bg), c(q2, Reduce(c, unname(tg2)))))
    skill_check("overlap", "background is broader than the tested sets", if (own > 0.9) "WARN" else "PASS",
                sprintf("%.0f%% of background regions overlap the query or a target%s", 100 * own,
                        if (own > 0.9) "; p-values will be conservative: rank by effect size only" else ""))
    ov <- skill_capture(enrichPeakOverlap(q2, tg2, TxDb = txdb, method = "fisher", background = bg,
                                          verbose = FALSE), "enrichPeakOverlap")
    rk <- attr(ov, "params")$rank_by
    ov <- ov[order(-ov[[rk]]), , drop = FALSE]
    skill_table(as.data.frame(ov), "overlap.tsv")
    skill_check("overlap", "ranked by effect size, not p-value", "PASS",
                paste0("rank_by = ", rk, "; ", attr(ov, "params")$null))
    ov
}

range_overlap_section <- function(ov) {
    if (is.null(ov)) return("not run (no --targets)")
    c(paste0("Background: ", .skill$notes$background, "; ranked by ", attr(ov, "params")$rank_by, "."),
      "", "| target | N_OL | jaccard | odds_ratio | p.adjust |", "|---|---|---|---|---|",
      paste0("| ", ov$tSample, " | ", ov$N_OL, " | ", signif(ov$jaccard, 3), " | ",
             signif(ov$odds_ratio, 3), " | ", signif(ov$p.adjust, 3), " |"))
}

## the Input section of a region-set report
range_input_section <- function(inp, anno, a) {
    src <- annoSource(anno[[1]])
    c(paste0("- ", inp$names, ": ", vapply(inp$peaks, length, 0L), " regions, median width ", inp$width, " bp"),
      paste0("- genome ", a$genome, "; TxDb ", src$TxDb$source, " (", src$TxDb$created, "); annoDb ",
             if (is.null(src$annoDb)) "none" else paste(src$annoDb$name, src$annoDb$version)))
}
