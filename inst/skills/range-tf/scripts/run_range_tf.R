## range-tf: narrow peaks (TF ChIP-seq, CUT&RUN, ATAC summits) -> input checks -> annotation,
## profile, comparison, overlap with reference sets in a background -> correctness checks ->
## interpret() last -> provenance. See ../SKILL.md.
##
## Rscript run_range_tf.R --example --out <dir> [--model provider:model]
## Rscript run_range_tf.R --peaks q.bed --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     --annodb org.Hs.eg.db --assay TF --context "..." [--question "..."] \
##     [--compare a.bed,b.bed --names q,a,b] [--targets dir_or_files --background open.bed] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

a <- skill_args(list(assay = "TF", tss = "3000"))
if (isTRUE(a$example)) {
    ## the package's own example: AR ChIP-seq peaks at 0, 1 and 100 nM androgen and CBX6/CBX7
    ## peaks (GEO GSM1174480-2, GSM1295076-7; lifted from hg19 to hg38 for the package)
    f <- getSampleFiles()
    a$peaks <- f$ARmo_1nM
    a$compare <- paste(f$ARmo_0M, f$ARmo_100nM, sep = ",")
    a$names <- "ARmo_1nM,ARmo_0M,ARmo_100nM"
    a$targets <- paste(unlist(f[c("ARmo_0M", "ARmo_100nM", "CBX6_BF", "CBX7_BF")]), collapse = ",")
    a$background <- "union_of_examples"
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$context <- "ChIP-seq in a prostate cancer cell line treated with 1 nM androgen"
    a$question <- "Which transcription factor was profiled?"
}
skill_start("range-tf", a, a$out)
tss <- as.numeric(a$tss)

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$peaks)) skill_check("input", "peaks given", "STOP", "give --peaks <file>")
if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP",
                                 "give --txdb <TxDb package>; annotateSeq() would otherwise fall back to hg38")
txdb <- skill_load(a$txdb)
skill_check_genome(a$genome, txdb)

files <- c(a$peaks, skill_split(a$compare))
nm <- skill_split(a$names)
if (is.null(nm)) nm <- sub("\\.(bed|narrowPeak|txt)(\\.gz)?$", "", basename(files))
if (length(nm) != length(files) || anyDuplicated(nm)) {
    skill_check("input", "one unique name per peak file", "STOP", "--names must name --peaks and --compare in order")
}
skill_input(files, "peaks")
peaks <- stats::setNames(lapply(files, readPeakFile), nm)
for (s in nm) {
    skill_check_recorded_genome(peaks[[s]], a$genome, what = s)
    peaks[[s]] <- skill_check_seqlevels(peaks[[s]], txdb, fix = a$fix_seqlevels, what = s)
}
## assay: TF peaks are narrow. Reference only, not a calibrated threshold: the median peak width
## of 20 ENCODE TF ChIP-seq sets (P1 pilot) was 184-556 bp.
w <- vapply(peaks, function(p) stats::median(GenomicRanges::width(p)), 0)
for (s in nm) {
    skill_check("input", paste("peak width of", s, "fits a narrow-peak assay"),
                if (w[[s]] > 556) "WARN" else "PASS",
                paste0("median ", w[[s]], " bp (20 ENCODE TF sets: 184-556 bp)",
                       if (w[[s]] > 556) "; wider than any of them: broad mark? range-tf may not fit" else ""))
}

## ---- 2. analysis with epiSeeker functions ----------------------------------------------
anno <- lapply(nm, function(s) skill_capture(
    annotateSeq(peaks[[s]], TxDb = txdb, annoDb = a$annodb, tssRegion = c(-tss, tss), verbose = FALSE),
    paste("annotateSeq", s)))
names(anno) <- nm
## correctness checks of the annotation come before anything is drawn from it
for (s in nm) skill_check_anno(anno[[s]], s, allow_suspect = a$allow_suspect)

q <- anno[[1]]
skill_table(as.data.frame(q), paste0("annotation_", nm[1], ".tsv"))
stat <- do.call(rbind, lapply(nm, function(s) cbind(set = s, getAnnoStat(anno[[s]]))))
skill_table(stat, "annotation_stat.tsv")
skill_figure(plotAnnoBar(csAnnoList(anno)), "annotation_bar.pdf")
skill_figure(plotDistToTSS(csAnnoList(anno)), "distance_to_tss.pdf")
tm <- skill_capture(getTagMatrix(peaks, upstream = tss, downstream = tss, type = "start_site",
                                 by = "gene", TxDb = txdb, verbose = FALSE), "getTagMatrix")
skill_figure(plotPeakProf(tm, facet = "none"), "tss_profile.pdf")

ov <- NULL
targets <- skill_split(a$targets)
if (length(targets) == 1 && dir.exists(targets)) targets <- list.files(targets, full.names = TRUE)
if (length(targets)) {
    skill_input(targets, "target")
    if (is.null(a$background)) {
        skill_check("overlap", "background given", "STOP", paste(
            "give --background: open chromatin (DNase/ATAC) of the same cell type, or the union of a large",
            "reference collection. A whole-genome null calls almost every pair significant",
            "(FPR 0.996 in the package's evaluation)"))
    }
    if (identical(a$background, "union_of_examples")) {
        bg <- Reduce(c, lapply(unlist(getSampleFiles()), readPeakFile))
        skill_note("background", "union of the five bundled example sets (example only)")
    } else {
        skill_input(a$background, "background")
        bg <- readPeakFile(a$background)
        skill_note("background", a$background)
    }
    bg <- skill_check_seqlevels(bg, txdb, fix = a$fix_seqlevels, what = "background")
    ## a universe made of little more than the tested sets: every unit is hit by something,
    ## so p-values are conservative (p near 1) -- rank, do not test, and say so
    tgt <- Reduce(c, lapply(targets, function(f) GenomicRanges::granges(readPeakFile(f))))
    own <- mean(IRanges::overlapsAny(GenomicRanges::reduce(bg), c(GenomicRanges::granges(q@anno), tgt)))
    skill_check("overlap", "background is broader than the tested sets", if (own > 0.9) "WARN" else "PASS",
                sprintf("%.0f%% of background regions overlap the query or a target%s", 100 * own,
                        if (own > 0.9) "; p-values will be conservative: rank by effect size only" else ""))
    ov <- skill_capture(enrichPeakOverlap(peaks[[1]], targets, TxDb = txdb, method = "fisher",
                                          background = bg, verbose = FALSE), "enrichPeakOverlap")
    rk <- attr(ov, "params")$rank_by
    ov <- ov[order(-ov[[rk]]), , drop = FALSE]
    skill_table(as.data.frame(ov), "overlap.tsv")
    skill_check("overlap", "ranked by effect size, not p-value", "PASS",
                paste0("rank_by = ", rk, "; ", attr(ov, "params")$null))
}

## ---- 3. report (numbers from the functions above) --------------------------------------
prom <- vapply(nm, function(s) { st <- getAnnoStat(anno[[s]]); sum(st$Frequency[grepl("^Promoter", st$Feature)]) }, 0)
src <- annoSource(q)
sec <- list(
    "Input" = c(paste0("- ", nm, ": ", vapply(peaks, length, 0L), " peaks, median width ", w, " bp"),
                paste0("- genome ", a$genome, "; TxDb ", src$TxDb$source, " (", src$TxDb$created, "); annoDb ",
                       if (is.null(src$annoDb)) "none" else paste(src$annoDb$name, src$annoDb$version))),
    "Annotation" = c(paste0("- ", nm, ": ", sprintf("%.1f", prom), "% in promoters (+/-", tss, " bp)"),
                     "- tables/annotation_stat.tsv, figures/annotation_bar.pdf, distance_to_tss.pdf, tss_profile.pdf"),
    "Overlap" = if (is.null(ov)) "not run (no --targets)" else
        c(paste0("Background: ", .skill$notes$background, "; ranked by ", attr(ov, "params")$rank_by, "."),
          "", "| target | N_OL | jaccard | odds_ratio | p.adjust |", "|---|---|---|---|---|",
          paste0("| ", ov$tSample, " | ", ov$N_OL, " | ", signif(ov$jaccard, 3), " | ",
                 signif(ov$odds_ratio, 3), " | ", signif(ov$p.adjust, 3), " |"))
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(q, model = a$model, context = a$context, assay = a$assay,
                       question = a$question, overlap = ov)
skill_flag_phrases(res, c(
    "causal claim from binding" = "\\b(regulates|activates|represses|drives|controls)\\b",
    "nearest gene called the target" = "\\btarget genes?\\b|\\bis the target\\b",
    "overlap as identity proof" = "\\b(proves?|confirms?|demonstrates?)\\b"))
skill_report("range-tf report", sec)
saveRDS(list(anno = anno, overlap = ov), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-tf done: ", skill_path("report", "report.md"))
