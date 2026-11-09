## range-tf: narrow peaks (TF ChIP-seq, CUT&RUN, CUT&Tag for TFs) -> input checks -> annotation,
## profile, comparison, overlap with reference sets in a background -> correctness checks ->
## interpret() last -> provenance. See ../SKILL.md. The steps shared with range-histone and
## range-atac are in ../../shared/range-steps.R.
##
## Rscript run_range_tf.R --example --out <dir> [--model provider:model]
## Rscript run_range_tf.R --peaks q.bed --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     --annodb org.Hs.eg.db --assay TF --context "..." [--question "..."] \
##     [--compare a.bed,b.bed --names q,a,b] [--targets dir_or_files --background open.bed] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
source(system.file("skills", "shared", "range-steps.R", package = "epiSeeker"))

a <- skill_args(list(assay = "TF", tss = "3000"))
if (isTRUE(a$example)) {
    ## the package's own example: AR ChIP-seq peaks of LNCaP-ARmo cells at 0, 1 and 100 nM DHT and
    ## CBX6/CBX7 peaks (GEO GSM1174480-2, GSM1295076-7; lifted from hg19 to hg38 for the package).
    ## Background: DNase-seq peaks of LNCaP after 1 nM R1881 (ENCODE ENCFF570MWX), chr1-chr2 only,
    ## so the overlap test runs on chr1-chr2 (inst/extdata/encode/SOURCE.txt)
    f <- getSampleFiles()
    a$peaks <- f$ARmo_1nM
    a$compare <- paste(f$ARmo_0M, f$ARmo_100nM, sep = ",")
    a$names <- "ARmo_1nM,ARmo_0M,ARmo_100nM"
    a$targets <- paste(unlist(f[c("ARmo_0M", "ARmo_100nM", "CBX6_BF", "CBX7_BF")]), collapse = ",")
    a$background <- system.file("extdata", "encode", "LNCaP_R1881_DNase_ENCFF570MWX_chr1-2.bed.gz", package = "epiSeeker")
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$context <- "ChIP-seq in a prostate cancer cell line treated with 1 nM androgen"
    a$question <- "Which transcription factor was profiled?"
}
skill_start("range-tf", a, a$out)
tss <- as.numeric(a$tss)

## ---- 1. input checks -------------------------------------------------------------------
inp <- range_inputs(a)
nm <- inp$names
## assay: TF peaks are narrow. Reference only, not a calibrated threshold: the median peak width
## of 20 ENCODE TF ChIP-seq sets (P1 pilot) was 184-556 bp.
for (s in nm) {
    w <- inp$width[[s]]
    skill_check("input", paste("peak width of", s, "fits a narrow-peak assay"), if (w > 556) "WARN" else "PASS",
                paste0("median ", w, " bp (20 ENCODE TF sets: 184-556 bp)",
                       if (w > 556) "; wider than any of them: broad mark? use range-histone" else ""))
}
## the declared assay against ENCODE peak sets (seqAssayQC; partial: it cannot tell a broad mark
## declared as TF, which the width check above can sometimes)
qc <- range_assay_qc(inp$peaks, inp$txdb, a$assay)

## ---- 2. analysis with epiSeeker functions (correctness checks of each annotation inside) ----
anno <- range_annotate(inp$peaks, inp$txdb, a$annodb, tss, allow_suspect = a$allow_suspect)
anno <- range_keep_assay_qc(anno, qc)
range_profile(inp$peaks, inp$txdb, tss, tss, "tss_profile.pdf")
ov <- range_overlap(inp$peaks[[1]], a, inp$txdb)

## ---- 3. report (numbers from the functions above) --------------------------------------
prom <- range_promoter_pct(anno)
sec <- list(
    "Input" = range_input_section(inp, anno, a),
    "Assay label" = range_assay_section(qc),
    "Annotation" = c(paste0("- ", nm, ": ", sprintf("%.1f", prom), "% in promoters (+/-", tss, " bp)"),
                     "- tables/annotation_stat.tsv, figures/annotation_bar.pdf, distance_to_tss.pdf, tss_profile.pdf"),
    "Overlap" = range_overlap_section(ov)
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(anno[[1]], model = a$model, context = a$context, assay = a$assay,
                       question = a$question, overlap = ov)
skill_flag_phrases(res, c(
    "assay label called verified" = "\\b(assay|label)\\b[^.]{0,60}\\b(verified|confirmed|validated)\\b",
    "causal claim from binding" = "\\b(regulates|activates|represses|drives|controls)\\b",
    "nearest gene called the target" = "\\btarget genes?\\b|\\bis the target\\b",
    "overlap as identity proof" = "\\b(proves?|confirms?|demonstrates?)\\b"))
skill_report("range-tf report", sec)
saveRDS(list(anno = anno, overlap = ov, assay_qc = qc), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-tf done: ", skill_path("report", "report.md"))
