## range-singlecell: clusters of single cells -> region set per group via the adapters
## (fragmentsToGroups / peaksToGroups) -> input checks -> annotateGroups, comparison ->
## correctness checks -> interpret() last -> provenance. See ../SKILL.md. Spatial spots go to
## range-spatial, which adds the depth / row / column artefact check. Shared steps:
## ../../shared/group-steps.R.
##
## Rscript run_range_singlecell.R --example --out <dir> [--model provider:model]
## Rscript run_range_singlecell.R --fragments frags.tsv.gz --groups barcode_group.csv \
##     | --peaks_table markers.tsv  --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     [--annodb org.Hs.eg.db] [--assay ATAC] [--tile 500] [--min_fragments 5] [--fold 2] \
##     [--context "..."] [--question "..."] [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
source(system.file("skills", "shared", "group-steps.R", package = "epiSeeker"))

a <- skill_args(list(tile = "500", min_fragments = "5", fold = "2", assay = "ATAC"))
if (isTRUE(a$example)) {
    ## toy scATAC fragments simulated from the bundled peaks (see group_simulate_cells())
    sim <- group_simulate_cells(file.path(a$out, "provenance"))
    a$fragments <- sim$fragments
    a$groups <- sim$groups
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$context <- "Simulated pseudo-single-cell ATAC fragments built from the package's example peaks"
}
skill_start("range-singlecell", a, a$out)

## ---- 1. input checks -------------------------------------------------------------------
if (!is.null(a$positions)) {
    skill_check("input", "spatial data go to range-spatial", "STOP",
                "--positions given: run range-spatial, which checks the domains against depth, chip rows and columns")
}
if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP", "give --txdb <TxDb package>")
txdb <- skill_load(a$txdb)
skill_check_genome(a$genome, txdb)
gin <- group_inputs(a, txdb)
group_size_checks(gin$regions, gin$members)

## ---- 2. analysis with epiSeeker functions (correctness checks inside) ----------------------
x <- group_annotate(gin$regions, txdb, a$annodb, "single-cell", gin$members, allow_suspect = a$allow_suspect)

## ---- 3. report ---------------------------------------------------------------------------
sec <- list(
    "Groups" = group_section(x),
    "Note" = paste("Each group is a pseudo-bulk region set (tiles enriched in the group, or the group's peaks);",
                   "differences between groups are differences between these sets, not between single cells.")
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(x, model = a$model, context = a$context, assay = a$assay, question = a$question)
skill_flag_phrases(res, c(
    "single-cell claim from a pseudo-bulk" = "\\b(every|each|all) (single )?(cells?|spots?)\\b",
    "causal claim from accessibility" = "\\b(regulates|drives|controls|activates|represses)\\b"))
skill_report("range-singlecell report", sec)
saveRDS(list(groups = gin$regions, anno = x), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-singlecell done: ", skill_path("report", "report.md"))
