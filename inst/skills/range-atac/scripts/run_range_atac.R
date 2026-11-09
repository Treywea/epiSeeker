## range-atac: ATAC-seq / DNase-seq peaks (and fragments, if given) -> input checks ->
## accessibility at promoters, distance to TSS, TSS profile, optional overlap in a background;
## with fragments: TSS enrichment, FRiP, fragment lengths -> correctness checks -> interpret()
## last -> provenance. See ../SKILL.md. Shared steps: ../../shared/range-steps.R.
##
## Rscript run_range_atac.R --example --out <dir> [--model provider:model]             (DNase peaks)
## Rscript run_range_atac.R --example_fragments --out <dir> [--model provider:model]   (+ fragments)
## Rscript run_range_atac.R --peaks peaks.narrowPeak --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     --annodb org.Hs.eg.db --assay ATAC [--fragments fragments.tsv.gz] [--compare b.bed --names a,b] \
##     [--targets refs --background universe.bed] [--context "..."] [--question "..."] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
source(system.file("skills", "shared", "range-steps.R", package = "epiSeeker"))

a <- skill_args(list(assay = "ATAC", tss = "3000"))
if (isTRUE(a$example) || isTRUE(a$example_fragments)) {
    ## DNase-seq peaks of LNCaP after 1 nM R1881 (ENCODE ENCSR000EKU, ENCFF570MWX), chr1-chr2 only
    ## (inst/extdata/encode/SOURCE.txt)
    a$peaks <- system.file("extdata", "encode", "LNCaP_R1881_DNase_ENCFF570MWX_chr1-2.bed.gz", package = "epiSeeker")
    a$names <- "LNCaP_R1881_DNase"
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$assay <- "DNase-seq"
    a$context <- "DNase-seq of LNCaP prostate cancer cells treated with 1 nM R1881, chromosomes 1 and 2 only"
    a$question <- "Where is chromatin accessible, and near which genes?"
}
if (isTRUE(a$example_fragments)) {
    ## 60,000 fragments simulated from those peaks (seed 1): 60% start inside a peak, the rest
    ## anywhere on chr1-chr2; lengths from a nucleosome-free (40-150 bp) and a mono-nucleosome
    ## (180-250 bp) mode. It shows the fragment path; its TSS enrichment and FRiP are not real.
    dir.create(file.path(a$out, "provenance"), recursive = TRUE, showWarnings = FALSE)
    pk <- readPeakFile(a$peaks)
    sl <- GenomeInfoDb::seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)[c("chr1", "chr2")]
    set.seed(1)
    n <- 60000; n_in <- 0.6 * n
    j <- sample(length(pk), n_in, replace = TRUE)
    s_in <- GenomicRanges::start(pk)[j] + floor(stats::runif(n_in) * GenomicRanges::width(pk)[j])
    ch <- sample(names(sl), n - n_in, replace = TRUE, prob = sl)
    s_out <- floor(stats::runif(n - n_in) * (sl[ch] - 1000))
    len <- ifelse(stats::runif(n) < 0.6, sample(40:150, n, replace = TRUE), sample(180:250, n, replace = TRUE))
    st <- c(s_in, s_out) - floor(len / 2)
    fr <- data.frame(chr = c(as.character(GenomicRanges::seqnames(pk))[j], ch), start = st, end = st + len,
                     barcode = "bulk", count = 1L)
    fr <- fr[order(fr$chr, fr$start), ]
    a$fragments <- file.path(a$out, "provenance", "example_fragments.tsv.gz")
    con <- gzfile(a$fragments, "w"); utils::write.table(fr, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(con)
    a$context <- paste(a$context, "(fragments simulated from the peaks for demonstration)")
}
skill_start("range-atac", a, a$out)
tss <- as.numeric(a$tss)

## ---- 1. input checks -------------------------------------------------------------------
inp <- range_inputs(a)
nm <- inp$names
for (s in nm) {
    w <- inp$width[[s]]
    skill_check("input", paste("peak width of", s, "for an accessibility assay"), if (w > 1000) "NOTE" else "PASS",
                paste0("median ", w, " bp", if (w > 1000) "; wider than usual ATAC/DNase peak calls: broad-mode caller or merged regions?" else ""))
}
## the declared assay against ENCODE peak sets (seqAssayQC; partial)
qc <- range_assay_qc(inp$peaks, inp$txdb, a$assay)
fr <- NULL
if (!is.null(a$fragments)) {
    skill_input(a$fragments, "fragments")
    fr <- skill_capture(readFragments(a$fragments), "readFragments")
    skill_check_recorded_genome(fr, a$genome, what = "fragments")
    fr <- skill_check_seqlevels(fr, inp$txdb, fix = a$fix_seqlevels, what = "fragments")
}

## ---- 2. analysis with epiSeeker functions (correctness checks of each annotation inside) ----
anno <- range_annotate(inp$peaks, inp$txdb, a$annodb, tss, allow_suspect = a$allow_suspect)
anno <- range_keep_assay_qc(anno, qc)
range_profile(inp$peaks, inp$txdb, tss, tss, "tss_profile.pdf")
fq <- NULL
if (!is.null(fr)) {
    ## Tn5 insertion sites = the two ends of each fragment; their profile around TSSs
    ins <- c(GenomicRanges::resize(GenomicRanges::granges(fr), 1, fix = "start"),
             GenomicRanges::resize(GenomicRanges::granges(fr), 1, fix = "end"))
    tm <- range_profile(ins, inp$txdb, 2000, 2000, "tss_insertions.pdf", step = "getTagMatrix insertions")
    prof <- colMeans(tm)
    k <- length(prof); mid <- (k + 1) / 2
    tsse <- mean(prof[(mid - 50):(mid + 50)]) / mean(c(prof[1:100], prof[(k - 99):k]))
    frip <- mean(IRanges::overlapsAny(fr, inp$peaks[[1]], ignore.strand = TRUE))
    fw <- GenomicRanges::width(fr)
    fq <- data.frame(fragments = length(fr), insertions = length(ins), tss_enrichment = round(tsse, 3),
                     frip = round(frip, 4), fragment_length_median = stats::median(fw),
                     fraction_below_147bp = round(mean(fw < 147), 4))
    skill_table(fq, "fragments_qc.tsv")
    skill_check("correctness", "TSS enrichment of insertions", if (tsse < 5) "NOTE" else "PASS",
                sprintf("%.2f (mean per bp at TSS +/-50 bp over the outer 100 bp of +/-2 kb); ENCODE ATAC, GRCh38: < 5 concerning, 5-7 acceptable, > 7 ideal (different TSS set and normalisation: a reference, not a threshold)", tsse))
    skill_check("correctness", "fraction of fragments in peaks (FRiP)", if (frip < 0.3) "NOTE" else "PASS",
                sprintf("%.3f of %d fragments overlap a peak of %s; ENCODE ATAC: > 0.3", frip, length(fr), nm[1]))
}
ov <- range_overlap(inp$peaks[[1]], a, inp$txdb)

## ---- 3. report ---------------------------------------------------------------------------
prom <- range_promoter_pct(anno)
dist <- vapply(nm, function(s) stats::median(abs(as.data.frame(anno[[s]])$distanceToTSS)), 0)
sec <- list(
    "Input" = range_input_section(inp, anno, a),
    "Assay label" = range_assay_section(qc),
    "Accessibility" = c(paste0("- ", nm, ": ", sprintf("%.1f", prom), "% in promoters (+/-", tss,
                               " bp); median |distance to TSS| ", dist, " bp"),
                        "- tables/annotation_stat.tsv, figures/annotation_bar.pdf, distance_to_tss.pdf, tss_profile.pdf"),
    "Fragments" = if (is.null(fq)) "not given" else
        c(paste0("- ", fq$fragments, " fragments; TSS enrichment ", fq$tss_enrichment, "; FRiP ", fq$frip,
                 "; median length ", fq$fragment_length_median, " bp, ", round(100 * fq$fraction_below_147bp, 1),
                 "% below 147 bp"), "- tables/fragments_qc.tsv, figures/tss_insertions.pdf"),
    "Overlap" = range_overlap_section(ov),
    "Reading" = paste("Accessible regions show where chromatin is open, not which factor binds there or",
                      "whether an element is active; a footprint or ChIP is needed for binding.")
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(anno[[1]], model = a$model, context = a$context, assay = a$assay,
                       question = a$question, overlap = ov)
skill_flag_phrases(res, c(
    "assay label called verified" = "\\b(assay|label)\\b[^.]{0,60}\\b(verified|confirmed|validated)\\b",
    "binding from accessibility" = "\\b(binds?|bound by|occupied by|binding sites?)\\b",
    "activity from accessibility" = "\\b(active enhancers?|activates|regulates|drives|controls)\\b",
    "nearest gene called the target" = "\\btarget genes?\\b|\\bis the target\\b"))
skill_report("range-atac report", sec)
saveRDS(list(anno = anno, fragments_qc = fq, overlap = ov, assay_qc = qc), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-atac done: ", skill_path("report", "report.md"))
