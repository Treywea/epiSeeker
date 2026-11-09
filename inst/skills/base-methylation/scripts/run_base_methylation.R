## base-methylation: single-base modification levels (bisulfite / EM-seq, nanopore or PacBio 5mC,
## Fiber-seq m6A) -> input checks -> bmData, clusters annotated, profile in motif context,
## per-read views -> correctness checks -> interpret() last -> provenance. See ../SKILL.md.
##
## Rscript run_base_methylation.R --example --out <dir> [--model provider:model]        (bulk)
## Rscript run_base_methylation.R --example_reads --out <dir> [--model provider:model]  (per read)
## Rscript run_base_methylation.R --input f1[,f2] --format bedmethyl|modcalls|table|rds \
##     --genome hg38 [--names s1,s2] [--mod m] [--type 5mC] [--txdb TxDb.Hsapiens.UCSC.hg38.knownGene] \
##     [--annodb org.Hs.eg.db] [--bsgenome BSgenome.Hsapiens.UCSC.hg38 --motif CG --region chr22:1-2] \
##     [--assay WGBS] [--context "..."] [--question "..."] [--min_cov 5] --out <dir> [--model ...]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

a <- skill_args(list(min_cov = "5", gap = "500", type = "5mC", motif = "CG"))
reads <- NULL
if (isTRUE(a$example) || isTRUE(a$example_reads)) {
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$bsgenome <- "BSgenome.Hsapiens.UCSC.hg38"
    a$region <- "chr22:10525991-10526342"
}
if (isTRUE(a$example)) {
    ## the bundled single-base table: 50 CpGs on chr22 of one sample (GEO GSE223153), the data of
    ## data(demo_bmdata), read here from its text form to show the table path
    a$input <- system.file("extdata", "demo_bisseq.txt", package = "epiSeeker")
    a$format <- "table"
    a$names <- "acinar_methyl"
    a$assay <- "5mC levels from bisulfite-type sequencing (bundled example)"
}
if (isTRUE(a$example_reads)) {
    ## a modkit `extract full` table simulated from the bundled levels: 40 reads across the region,
    ## each CpG modified with the probability given by demo_bisseq.txt (seed 1).
    ## It shows the per-read path; it is not real nanopore data.
    lv <- utils::read.delim(system.file("extdata", "demo_bisseq.txt", package = "epiSeeker"))
    set.seed(1)
    win <- list(c(10525900, 10526450))
    rows <- list()
    for (w in seq_along(win)) for (r in seq_len(40)) {
        st <- win[[w]][1] + sample(0:60, 1); en <- win[[w]][2] - sample(0:60, 1)
        s <- lv[lv$pos >= st & lv$pos <= en, ]
        q <- ifelse(stats::runif(nrow(s)) < s$Methylation, stats::runif(nrow(s), 0.8, 1), stats::runif(nrow(s), 0, 0.2))
        rows[[length(rows) + 1]] <- data.frame(
            read_id = sprintf("read%d_%02d", w, r), forward_read_position = s$pos - st, ref_position = s$pos - 1,
            chrom = s$chr, mod_strand = "+", ref_strand = "+", ref_mod_strand = "+", fw_soft_clipped_start = 0,
            fw_soft_clipped_end = 0, alignment_start = st - 1, alignment_end = en, read_length = en - st + 1,
            mod_qual = round(q, 3), mod_code = "m", base_qual = 30, ref_kmer = "ACGTA", query_kmer = "ACGTA",
            canonical_base = "C", modified_primary_base = "C", inferred = "false", flag = 0)
    }
    a$out_input <- file.path(a$out, "provenance", "example_modkit_extract.tsv")
    dir.create(dirname(a$out_input), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(do.call(rbind, rows), a$out_input, sep = "\t", quote = FALSE, row.names = FALSE)
    a$input <- a$out_input
    a$format <- "modcalls"
    a$names <- "simulated_reads"
    a$assay <- "5mC per-read calls (simulated from the bundled levels)"
}
skill_start("base-methylation", a, a$out)
min_cov <- as.numeric(a$min_cov)

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$input) || is.null(a$format)) skill_check("input", "input and format given", "STOP",
                                                     "give --input <file(s)> and --format bedmethyl|modcalls|table|rds")
if (is.null(a$genome) || isTRUE(a$genome)) skill_check("input", "genome build declared", "STOP", "give --genome (e.g. hg38)")
txdb <- skill_load(a$txdb)
if (!is.null(txdb)) skill_check_genome(a$genome, txdb)
files <- skill_split(a$input)
skill_input(files, "input")
nm <- skill_split(a$names)
if (is.null(nm)) nm <- sub("\\..*$", "", basename(files))
bm <- switch(a$format,
    bedmethyl = {
        if (length(files) > 1) skill_check("input", "one bedMethyl file per run", "STOP", "run once per sample")
        skill_capture(readBedMethyl(files, mod = a$mod, sampleName = nm), "readBedMethyl")
    },
    modcalls = {
        if (length(files) > 1) skill_check("input", "one per-read table per run", "STOP", "run once per sample")
        reads <- skill_capture(readModCalls(files), "readModCalls")
        src <- S4Vectors::metadata(reads)$modSource
        skill_check("input", "per-read calls placed on the reference", if (src$n_unplaced > 0) "NOTE" else "PASS",
                    paste0(src$n_unplaced, " calls without a reference position dropped; format ", src$format))
        skill_check("input", paste("call type", a$type, "present"), if (a$type %in% reads$type) "PASS" else "STOP",
                    paste("types:", paste(unique(reads$type), collapse = ", ")))
        skill_capture(summarizeModCalls(reads, type = a$type, sampleName = nm), "summarizeModCalls")
    },
    table = {
        ## columns chr, pos, coverage, level -- as inst/extdata/demo_bisseq.txt
        dfs <- lapply(files, function(f) utils::read.delim(f))
        if (any(vapply(dfs, function(d) ncol(d) < 4, TRUE))) {
            skill_check("input", "table columns chr, pos, coverage, level", "STOP", "need four columns")
        }
        skill_capture(makeBmDataFromData(if (length(dfs) == 1) dfs[[1]] else dfs, sampleNames = nm), "makeBmDataFromData")
    },
    rds = readRDS(files),
    skill_check("input", "known --format", "STOP", a$format))
if (!is(bm, "bmData")) skill_check("input", "input became a bmData", "STOP", class(bm)[1])

gr <- SummarizedExperiment::rowRanges(bm)
an <- SummarizedExperiment::assayNames(bm)
skill_check("input", "sites read", "PASS", paste0(length(gr), " sites, ", ncol(bm), " sample(s), assays ", paste(an, collapse = ", ")))
skill_check_recorded_genome(gr, a$genome, what = "sites")
if (!is.null(txdb)) {
    st <- GenomeInfoDb::seqlevelsStyle(gr); tt <- GenomeInfoDb::seqlevelsStyle(txdb)
    skill_check("input", "seqlevels style of sites", if (length(intersect(st, tt))) "PASS" else "STOP",
                paste0(paste(st, collapse = "/"), " vs TxDb ", paste(tt, collapse = "/")))
}
## level scale: interpret() and the plots expect a fraction in [0, 1]; percentages (0-100) are
## the commonest silent error here
rng <- vapply(an, function(x) { v <- as.matrix(SummarizedExperiment::assay(bm, x)); v <- v[is.finite(v)]
                                if (length(v)) c(min(v), max(v)) else c(NA, NA) }, numeric(2))
lev <- an[!is.na(rng[1, ]) & rng[1, ] >= 0 & rng[2, ] <= 1]
if (!length(lev)) {
    skill_check("input", "a level assay in [0, 1]", "STOP", paste0(
        paste0(an, " ", rng[1, ], "-", rng[2, ], collapse = "; "),
        if (any(rng[2, ] > 1 & rng[2, ] <= 100, na.rm = TRUE)) "; looks like percent: divide by 100" else ""))
}
lev <- lev[1]
skill_check("input", "a level assay in [0, 1]", "PASS", lev)
covn <- setdiff(an, lev)
covn <- covn[vapply(covn, function(x) { v <- as.matrix(SummarizedExperiment::assay(bm, x)); v <- v[is.finite(v)]
                                         length(v) && all(v >= 0 & v == round(v)) }, TRUE)]
if (!length(covn)) {
    skill_check("input", "coverage assay", "WARN", "none: levels of 1 read count as much as levels of 100")
} else {
    cv <- as.matrix(SummarizedExperiment::assay(bm, covn[1]))
    low <- mean(cv < min_cov, na.rm = TRUE)
    skill_check("input", paste("coverage >=", min_cov), if (low > 0.5) "WARN" else if (low > 0) "NOTE" else "PASS",
                sprintf("%.0f%% of site x sample values below %d (left out of levels by interpret())", 100 * low, min_cov))
}
if (!is.null(reads)) {
    fmt <- S4Vectors::metadata(reads)$modSource$format
    skill_check("input", "how coverage of per-read calls is counted", "NOTE", if (identical(fmt, "fibertools"))
        paste("fibertools lists only calls above its threshold: coverage = reads spanning the site; for 5mC",
              "only reads of the strand that gave calls there (PacBio calls CpGs on the read's own strand)") else
        "the table lists every tested base: coverage = calls at the site (see ?summarizeModCalls)")
}

## ---- 2. analysis with epiSeeker functions ----------------------------------------------
gap <- as.numeric(a$gap)
lv_tab <- cbind(as.data.frame(gr)[, 1:2], do.call(cbind, lapply(an, function(x)
    stats::setNames(as.data.frame(as.matrix(SummarizedExperiment::assay(bm, x))), paste0(x, "_", colnames(bm))))))
skill_table(lv_tab, "sites.tsv")
anno <- NULL
if (!is.null(txdb)) {
    cl <- GenomicRanges::reduce(GenomicRanges::granges(gr), min.gapwidth = gap + 1L, ignore.strand = TRUE)
    anno <- skill_capture(annotateSeq(cl, TxDb = txdb, annoDb = a$annodb, verbose = FALSE), "annotateSeq clusters")
    skill_check_anno(anno, paste0("site clusters (within ", gap, " bp)"), allow_suspect = a$allow_suspect)
    skill_table(as.data.frame(anno), "clusters_annotation.tsv")
    skill_figure(plotAnnoBar(anno), "clusters_annotation_bar.pdf")
}
if (!is.null(a$region)) {
    rg <- GenomicRanges::GRanges(a$region)
    if (!is.null(a$bsgenome)) {
        bs <- skill_load(a$bsgenome)
        g <- unique(stats::na.omit(GenomeInfoDb::genome(bs)))
        skill_check("input", "BSgenome matches --genome", if (a$genome %in% g) "PASS" else "STOP",
                    paste0(a$bsgenome, ": ", paste(g, collapse = ",")))
        m <- skill_capture(getBmMatrix(region = data.frame(chr = as.character(GenomicRanges::seqnames(rg)),
                                                          start = GenomicRanges::start(rg), end = GenomicRanges::end(rg)),
                                       BSgenome = bs, input = bm, base = "C", motif = skill_split(a$motif)), "getBmMatrix")
        skill_figure(plotBmProf(m), "profile_region.pdf", width = 10, height = 6)
    }
    if (!is.null(reads)) {
        skill_figure(plotModReads(reads, rg, types = a$type), "reads_region.pdf", width = 8, height = 6)
        ## two halves of the region: modified together on the same molecules?
        mid <- GenomicRanges::start(rg) + GenomicRanges::width(rg) %/% 2
        h1 <- GenomicRanges::GRanges(GenomicRanges::seqnames(rg), IRanges::IRanges(GenomicRanges::start(rg), mid))
        h2 <- GenomicRanges::GRanges(GenomicRanges::seqnames(rg), IRanges::IRanges(mid + 1, GenomicRanges::end(rg)))
        co <- modCooccurrence(reads, h1, h2, type = a$type)
        skill_table(data.frame(region1 = as.character(h1), region2 = as.character(h2), n_reads = co$n_reads,
                               odds_ratio = co$odds_ratio, p_value = co$p_value), "cooccurrence.tsv")
        skill_check("correctness", "reads span both halves for the co-occurrence test",
                    if (is.na(co$p_value)) "NOTE" else "PASS",
                    paste0(co$n_reads, " reads span both",
                           if (is.na(co$p_value) && co$n_reads >= 2) "; not testable: a margin of the 2x2 table is empty (e.g. every read modified in one half)"
                           else if (is.na(co$p_value)) "; not testable: fewer than 2 reads cover both halves" else ""))
    }
}

## ---- 3. report ---------------------------------------------------------------------------
sec <- list(
    "Input" = c(paste0("- ", ncol(bm), " sample(s): ", paste(colnames(bm), collapse = ", "), "; ", length(gr),
                       " sites; level assay `", lev, "`", if (length(covn)) paste0(", coverage `", covn[1], "`") else ""),
                paste0("- genome ", a$genome, if (!is.null(txdb)) paste0("; TxDb ", annoSource(anno)$TxDb$source) else "")),
    "Analysis" = c(if (!is.null(anno)) paste0("- ", length(anno@anno), " clusters of sites (gap ", gap, " bp) annotated; tables/clusters_annotation.tsv"),
                   if (!is.null(a$region)) paste0("- region ", a$region, ": figures/profile_region.pdf",
                                                  if (!is.null(reads)) ", reads_region.pdf, tables/cooccurrence.tsv" else ""))
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(bm, model = a$model, context = a$context, assay = a$assay, question = a$question,
                       value = lev, coverage = if (length(covn)) covn[1] else NULL, TxDb = txdb,
                       annoDb = a$annodb, gap = gap, min_cov = min_cov)
skill_flag_phrases(res, c(
    "silencing from levels alone" = "\\b(silenc(ed|es|ing)|represses|switch(es|ed)? off)\\b",
    "difference without replicates called significant" = "\\bsignificant(ly)?\\b"))
skill_report("base-methylation report", sec)
saveRDS(list(bmData = bm, clusters = anno), skill_path("provenance", "objects.rds"))
skill_finish()
message("base-methylation done: ", skill_path("report", "report.md"))
