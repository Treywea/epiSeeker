## range-histone: broad histone marks (H3K27me3, H3K36me3, H3K9me3, H3K79me2, ...) -> input
## checks -> peaks merged into domains -> domain summaries, genes under domains, annotation,
## profiles over wide windows and gene bodies, optional overlap in a background -> correctness
## checks -> interpret() last -> provenance. See ../SKILL.md. Shared steps: ../../shared/range-steps.R.
##
## Rscript run_range_histone.R --example --out <dir> [--model provider:model]
## Rscript run_range_histone.R --peaks mark.broadPeak --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     --annodb org.Hs.eg.db --assay H3K27me3 [--merge_gap 5000] [--window 10000] \
##     [--compare b.bed --names a,b] [--targets refs --background universe.bed] [--context "..."] \
##     [--question "..."] [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
source(system.file("skills", "shared", "range-steps.R", package = "epiSeeker"))

a <- skill_args(list(assay = "broad histone mark", tss = "3000", merge_gap = "5000", window = "10000"))
if (isTRUE(a$example)) {
    ## H3K27me3 replicated peaks of adult prostate gland tissue (ENCODE ENCSR690CSD, ENCFF198QRW),
    ## chr1-chr2 only (inst/extdata/encode/SOURCE.txt). Peaks of the ENCODE histone pipeline are
    ## called in narrow mode, so merging them into domains is the first step.
    a$peaks <- system.file("extdata", "encode", "prostate_H3K27me3_ENCFF198QRW_chr1-2.bed.gz", package = "epiSeeker")
    a$names <- "prostate_H3K27me3"
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$assay <- "H3K27me3 (broad histone mark)"
    a$context <- "H3K27me3 ChIP-seq of adult human prostate tissue, chromosomes 1 and 2 only"
    a$question <- "Which genes and kinds of genomic regions does this mark cover?"
}
skill_start("range-histone", a, a$out)
tss <- as.numeric(a$tss)
gap <- as.numeric(a$merge_gap)
win <- as.numeric(a$window)

## ---- 1. input checks -------------------------------------------------------------------
inp <- range_inputs(a)
nm <- inp$names
## broad marks: peaks kilobases wide. Peaks as narrow as TF peaks (<= 556 bp, the widest median
## of 20 ENCODE TF sets) come from a narrow-mode caller; merging into domains handles them.
for (s in nm) {
    w <- inp$width[[s]]
    skill_check("input", paste("peak width of", s, "for a broad mark"), if (w <= 556) "NOTE" else "PASS",
                paste0("median ", w, " bp", if (w <= 556)
                    paste0("; as narrow as TF peaks (20 ENCODE TF sets: 184-556 bp): narrow-mode caller? ",
                           "domains are built by merging peaks within ", gap, " bp") else ""))
}

## ---- 2. domains and their summaries ---------------------------------------------------------
dom <- lapply(inp$peaks, function(p) GenomicRanges::reduce(GenomicRanges::granges(p), min.gapwidth = gap + 1L,
                                                            ignore.strand = TRUE))
dsum <- do.call(rbind, lapply(nm, function(s) {
    dw <- GenomicRanges::width(dom[[s]])
    q <- stats::quantile(dw, c(0.25, 0.5, 0.75), names = FALSE)
    data.frame(set = s, peaks = length(inp$peaks[[s]]), domains = length(dw), merge_gap = gap,
               domain_width_q1 = q[1], domain_width_median = q[2], domain_width_q3 = q[3],
               domain_width_max = max(dw), domains_over_10kb = sum(dw >= 1e4), mb_covered = round(sum(dw) / 1e6, 3))
}))
skill_table(dsum, "domains_summary.tsv")
for (s in nm) skill_table(as.data.frame(dom[[s]])[, 1:4], paste0("domains_", s, ".tsv"))
skill_check("analysis", "peaks merged into domains", "PASS",
            paste0(dsum$set, ": ", dsum$peaks, " peaks -> ", dsum$domains, " domains (median ",
                   dsum$domain_width_median, " bp)", collapse = "; "))

## annotation of the domains (correctness checks inside), genes under domains
anno <- range_annotate(dom, inp$txdb, a$annodb, tss, allow_suspect = a$allow_suspect, what = "domains of ")
genes <- lapply(nm, function(s) skill_capture(seq2gene(dom[[s]], tssRegion = c(-tss, tss), flankDistance = 0,
                                                       TxDb = inp$txdb), paste("seq2gene", s)))
names(genes) <- nm
gtab <- do.call(rbind, lapply(nm, function(s) if (length(genes[[s]])) data.frame(set = s, gene_id = sort(genes[[s]]))))
if (!is.null(a$annodb) && !is.null(gtab)) {
    org <- skill_load(a$annodb)
    gtab$symbol <- unname(suppressMessages(AnnotationDbi::mapIds(org, gtab$gene_id, "SYMBOL", "ENTREZID")))
}
if (!is.null(gtab)) skill_table(gtab, "genes_under_domains.tsv")

## profiles: wide window around TSSs, and scaled gene bodies with 20% flanks
range_profile(dom, inp$txdb, win, win, "tss_profile_wide.pdf")
range_profile(dom, inp$txdb, rel(0.2), rel(0.2), "gene_body_profile.pdf", type = "body", nbin = 100,
              step = "getTagMatrix gene body")

ov <- range_overlap(dom[[1]], a, inp$txdb)

## ---- 3. report ---------------------------------------------------------------------------
prom <- range_promoter_pct(anno)
sec <- list(
    "Input" = range_input_section(inp, anno, a),
    "Domains" = c(paste0("- ", dsum$set, ": ", dsum$peaks, " peaks merged within ", gap, " bp into ", dsum$domains,
                         " domains; width median ", dsum$domain_width_median, " bp (IQR ", dsum$domain_width_q1, "-",
                         dsum$domain_width_q3, ", max ", dsum$domain_width_max, "), ", dsum$domains_over_10kb,
                         " of >= 10 kb; ", dsum$mb_covered, " Mb covered"),
                  paste0("- ", nm, ": ", lengths(genes), " genes with a domain over the promoter (+/-", tss,
                         " bp) or the gene body (seq2gene(); tables/genes_under_domains.tsv)"),
                  paste0("- ", nm, ": ", sprintf("%.1f", prom), "% of domains annotated to a promoter")),
    "Figures" = "annotation_bar, distance_to_tss (domains), tss_profile_wide (+/-window), gene_body_profile",
    "Overlap" = range_overlap_section(ov),
    "Reading" = paste("Domains are merged peaks; their genes are genes under a domain, not genes the mark",
                      "regulates. No motif or binding-site statement follows from a histone mark.")
)

## ---- 4. interpret() last -----------------------------------------------------------------
ctx <- paste(c(if (!is.null(a$context)) sub("([^.])$", "\\1.", a$context),
               paste0("Regions are broad domains merged from peaks within ", gap, " bp, not binding sites.")), collapse = " ")
res <- skill_interpret(anno[[1]], model = a$model, context = ctx, assay = a$assay, question = a$question, overlap = ov)
skill_flag_phrases(res, c(
    "TF-style claim from a broad mark" = "\\b(motifs?|binding sites?|binds?|bound by)\\b",
    "causal claim from a mark" = "\\b(regulates|silences|represses|activates|drives|controls)\\b",
    "nearest gene called the target" = "\\btarget genes?\\b|\\bis the target\\b"))
skill_report("range-histone report", sec)
saveRDS(list(domains = dom, anno = anno, genes = genes, overlap = ov), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-histone done: ", skill_path("report", "report.md"))
