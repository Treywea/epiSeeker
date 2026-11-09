## interaction-loops: chromatin loops (Hi-C / HiChIP / ChIA-PET BEDPE) and optional peaks ->
## input checks -> loop annotation, peak-to-gene links -> correctness checks -> interpret() last
## -> provenance. See ../SKILL.md. A loop partner is a candidate gene, never "the target".
##
## Rscript run_interaction_loops.R --example --out <dir> [--model provider:model]
## Rscript run_interaction_loops.R --loops loops.bedpe --genome hg38 \
##     --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db [--peaks peaks.bed] \
##     [--protein_coding] [--strict] [--context "..."] [--question "..."] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

## what the paper's Fig 4 found; sent with the context and written into the report
caveat <- paste(
    "Loop partners are candidate genes, not established targets: on a CRISPRi enhancer-gene",
    "benchmark (K562), genes linked to elements by real loops were regulated no more often than",
    "genes linked by distance-matched rewired loops (precision 0.220 vs 0.219); loops added",
    "recall over the nearest gene, not specificity.")

a <- skill_args(list())
if (isTRUE(a$example)) {
    ## toy loops from the bundled peaks: each AR peak (ARmo_1nM) joined to the nearest CBX7 peak
    ## 20 kb - 1 Mb away on the same chromosome, anchors 5 kb. They show the mechanics only and
    ## carry no biology.
    f <- getSampleFiles()
    c1 <- GenomicRanges::resize(GenomicRanges::granges(readPeakFile(f$ARmo_1nM)), 1, fix = "center")
    c2 <- GenomicRanges::resize(GenomicRanges::granges(readPeakFile(f$CBX7_BF)), 1, fix = "center")
    h <- GenomicRanges::findOverlaps(c1, c2, maxgap = 1e6, ignore.strand = TRUE)
    d <- GenomicRanges::distance(c1[S4Vectors::queryHits(h)], c2[S4Vectors::subjectHits(h)])
    h <- h[d >= 2e4]; d <- d[d >= 2e4]
    o <- order(S4Vectors::queryHits(h), d); h <- h[o]
    h <- h[!duplicated(S4Vectors::queryHits(h))]
    a1 <- GenomicRanges::resize(c1[S4Vectors::queryHits(h)], 5000, fix = "center")
    a2 <- GenomicRanges::resize(c2[S4Vectors::subjectHits(h)], 5000, fix = "center")
    a$loops <- file.path(a$out, "provenance", "example_loops.bedpe")
    dir.create(dirname(a$loops), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(data.frame(chrom1 = as.character(GenomicRanges::seqnames(a1)), start1 = GenomicRanges::start(a1) - 1,
                                  end1 = GenomicRanges::end(a1), chrom2 = as.character(GenomicRanges::seqnames(a2)),
                                  start2 = GenomicRanges::start(a2) - 1, end2 = GenomicRanges::end(a2)),
                       a$loops, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    a$peaks <- f$ARmo_1nM
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$protein_coding <- TRUE
    a$strict <- TRUE
    a$context <- "Toy loops built from the package's example peaks (AR and CBX7 ChIP-seq), for demonstration"
}
skill_start("interaction-loops", a, a$out)

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$loops)) skill_check("input", "loops given", "STOP", "give --loops <BEDPE>")
if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP", "give --txdb <TxDb package>")
txdb <- skill_load(a$txdb)
skill_check_genome(a$genome, txdb)
skill_input(a$loops, "loops")
gi <- skill_capture(readLoops(a$loops), "readLoops")
src <- S4Vectors::metadata(gi)$loopSource
skill_check("input", "loop rows with coordinates", if (src$dropped > 0) "NOTE" else "PASS",
            paste0(src$loops, " loops read, ", src$dropped, " rows without coordinates dropped"))
reg <- InteractionSet::regions(gi)
skill_check_recorded_genome(reg, a$genome, what = "loop anchors")
st <- GenomeInfoDb::seqlevelsStyle(reg); tt <- GenomeInfoDb::seqlevelsStyle(txdb)
if (!length(intersect(st, tt))) {
    if (!isTRUE(a$fix_seqlevels)) {
        skill_check("input", "seqlevels style of loop anchors", "STOP",
                    paste0(paste(st, collapse = "/"), " vs TxDb ", paste(tt, collapse = "/"), "; --fix_seqlevels converts (recorded)"))
    }
    GenomeInfoDb::seqlevelsStyle(gi) <- tt[1]
    skill_check("input", "seqlevels style of loop anchors", "NOTE", paste("converted to", tt[1]))
} else skill_check("input", "seqlevels style of loop anchors", "PASS", paste(st, collapse = "/"))
peak <- NULL
if (!is.null(a$peaks)) {
    skill_input(a$peaks, "peaks")
    peak <- readPeakFile(a$peaks)
    skill_check_recorded_genome(peak, a$genome, what = "peaks")
    peak <- skill_check_seqlevels(peak, txdb, fix = a$fix_seqlevels, what = "peaks")
}
gene_ids <- NULL
if (isTRUE(a$protein_coding)) {
    if (is.null(a$annodb)) skill_check("input", "annoDb for --protein_coding", "STOP", "give --annodb (an org.*.eg.db)")
    org <- skill_load(a$annodb)
    gt <- suppressMessages(AnnotationDbi::select(org, AnnotationDbi::keys(org), "GENETYPE"))
    gene_ids <- unique(gt$ENTREZID[gt$GENETYPE %in% "protein-coding"])
    skill_note("gene_filter", paste(length(gene_ids), "protein-coding genes of", a$annodb))
}

## ---- 2. analysis with epiSeeker functions ----------------------------------------------
la <- skill_capture(annotateLoops(gi, TxDb = txdb, annoDb = a$annodb, gene_ids = gene_ids, verbose = FALSE),
                    "annotateLoops")
skill_check_anno(la@anchors, "loop anchors", allow_suspect = a$allow_suspect)
lp <- as.data.frame(la)
skill_table(lp, "loops_annotated.tsv")
skill_table(la@loopStat, "loop_classes.tsv")
inter <- sum(is.na(lp$span))
short <- sum(!is.na(lp$span) & lp$span < 2 * stats::median(GenomicRanges::width(InteractionSet::regions(la@loops))))
skill_check("correctness", "loops between chromosomes / shorter than two anchors", if (inter + short > 0) "NOTE" else "PASS",
            paste0(inter, " between chromosomes, ", short, " with anchors closer than two anchor widths"))
skill_figure(plotLoopAnno(la, type = "class"), "loop_classes.pdf")
skill_figure(plotLoopAnno(la, type = "span"), "loop_span.pdf")

links <- NULL
if (!is.null(peak)) {
    pa <- skill_capture(annotateSeq(peak, TxDb = txdb, annoDb = a$annodb, verbose = FALSE), "annotateSeq peaks")
    skill_check_anno(pa, "peaks", allow_suspect = a$allow_suspect)
    same <- identical(annoSource(pa)$TxDb, annoSource(la@anchors)$TxDb) && identical(annoSource(pa)$genome, annoSource(la@anchors)$genome)
    skill_check("correctness", "peaks and loops annotated with the same TxDb", if (same) "PASS" else "STOP",
                paste(annoSource(pa)$TxDb$source, "/", annoSource(la@anchors)$TxDb$source))
    strict <- if (isTRUE(a$strict)) list(TxDb = txdb, anchor_flank = 1000, distal_distance = 3000, exclude_own_anchor = TRUE) else list()
    links <- skill_capture(do.call(linkPeaksToGenes, c(list(peak = pa, loops = la, gene_ids = gene_ids,
                                                            annoDb = a$annodb), strict)), "linkPeaksToGenes")
    skill_note("link_definition", if (isTRUE(a$strict)) "strict: anchor_flank 1000, distal_distance 3000, exclude_own_anchor (paper Fig 4)" else "annotateLoops() promoter anchors")
    skill_table(links, "peak_links.tsv")
    skill_figure(plotLoopLinks(links), "peak_links.pdf")
}

## ---- 3. report: candidates, never targets ---------------------------------------------------
lt <- if (is.null(links)) NULL else table(links$link)
sec <- list(
    "Input" = c(paste0("- ", length(la@loops), " loops, ", length(InteractionSet::regions(la@loops)), " anchors; genome ", a$genome,
                       "; TxDb ", annoSource(la@anchors)$TxDb$source),
                if (!is.null(gene_ids)) paste0("- genes counted: ", .skill$notes$gene_filter)),
    "Loop classes" = paste0("- ", la@loopStat$loop_class, ": ", la@loopStat$n, " (", round(la@loopStat$Frequency, 1), "%)"),
    "Peak links" = if (is.null(lt)) "not run (no --peaks)" else
        c(paste0("- ", names(lt), ": ", as.integer(lt)),
          "", "Genes in `loop_genes` of tables/peak_links.tsv are **loop-linked candidate genes**.", "", caveat)
)

## ---- 4. interpret() last -----------------------------------------------------------------
ctx <- paste(c(if (!is.null(a$context)) sub("([^.])$", "\\1.", a$context), caveat), collapse = " ")
res <- skill_interpret(la, model = a$model, context = ctx, question = a$question, links = links)
skill_flag_phrases(res, c(
    "loop partner called the target" = "\\btarget genes?\\b|\\bis the target\\b|\\btargets\\b",
    "regulation from a loop" = "\\b(regulates|controls|drives|activates|represses)\\b",
    "loop as proof" = "\\b(proves?|confirms?|demonstrates?)\\b"))
skill_report("interaction-loops report", sec)
saveRDS(list(loops = la, links = links), skill_path("provenance", "objects.rds"))
skill_finish()
message("interaction-loops done: ", skill_path("report", "report.md"))
