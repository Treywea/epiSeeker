## range-singlecell: clusters of single cells or spatial domains -> region set per group via the
## adapters (fragmentsToGroups / peaksToGroups) -> input checks -> annotateGroups, comparison ->
## correctness checks -> interpret() last -> provenance. See ../SKILL.md.
##
## Rscript run_range_singlecell.R --example --out <dir> [--model provider:model]          (scATAC)
## Rscript run_range_singlecell.R --example_spatial --out <dir> [--model provider:model]  (spatial)
## Rscript run_range_singlecell.R --fragments frags.tsv.gz --groups barcode_group.csv \
##     | --peaks_table markers.tsv  --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     [--annodb org.Hs.eg.db] [--assay ATAC] [--tile 500] [--min_fragments 5] [--fold 2] \
##     [--positions tissue_positions_list.csv] [--context "..."] [--question "..."] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

a <- skill_args(list(tile = "500", min_fragments = "5", fold = "2", assay = "ATAC"))
if (isTRUE(a$example) || isTRUE(a$example_spatial)) {
    ## toy fragments simulated from the bundled peak sets (seed 1): cells of group "AR_like" put 70%
    ## of 300 fragments in AR peaks (ARmo_1nM), "PcG_like" in CBX7 peaks, the small group
    ## "AR_0M_like" (20 cells) in ARmo_0M peaks; the rest fall anywhere on chr1-22. Fragment
    ## barcodes carry "-1", the group table does not. It shows the mechanics; it is not real data.
    dir.create(file.path(a$out, "provenance"), recursive = TRUE, showWarnings = FALSE)
    f <- getSampleFiles()
    sl <- GenomeInfoDb::seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)[paste0("chr", 1:22)]
    spec <- list(AR_like = list(f$ARmo_1nM, 60), PcG_like = list(f$CBX7_BF, 60), AR_0M_like = list(f$ARmo_0M, 20))
    set.seed(1)
    fr <- list(); grp <- list(); k <- 0
    for (g in names(spec)) {
        pk <- readPeakFile(spec[[g]][[1]])
        for (i in seq_len(spec[[g]][[2]])) {
            k <- k + 1; bc <- sprintf("CELL%04d", k)
            n_in <- 210; n_out <- 90
            j <- sample(length(pk), n_in, replace = TRUE)
            s_in <- GenomicRanges::start(pk)[j] + floor(stats::runif(n_in) * GenomicRanges::width(pk)[j])
            ch <- sample(names(sl), n_out, replace = TRUE, prob = sl)
            s_out <- floor(stats::runif(n_out) * (sl[ch] - 1000))
            len <- sample(150:300, n_in + n_out, replace = TRUE)
            fr[[k]] <- data.frame(chr = c(as.character(GenomicRanges::seqnames(pk))[j], ch),
                                  start = c(s_in, s_out) - 1, end = c(s_in, s_out) - 1 + len,
                                  barcode = paste0(bc, "-1"), count = 1L)
            grp[[k]] <- data.frame(barcode = bc, group = g)
        }
    }
    fr <- do.call(rbind, fr); fr <- fr[order(fr$chr, fr$start), ]
    a$fragments <- file.path(a$out, "provenance", "example_fragments.tsv.gz")
    con <- gzfile(a$fragments, "w"); utils::write.table(fr, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(con)
    grp <- do.call(rbind, grp)
    a$groups <- file.path(a$out, "provenance", "example_groups.csv")
    utils::write.csv(grp, a$groups, row.names = FALSE, quote = FALSE)
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$context <- "Simulated pseudo-single-cell ATAC fragments built from the package's example peaks"
    if (isTRUE(a$example_spatial)) {
        ## the same cells as spots on a 14 x 10 grid: AR_like on the left, PcG_like on the right,
        ## AR_0M_like in the bottom-right corner
        xy <- expand.grid(col = 0:13, row = 0:9)
        ord <- c(which(xy$col < 7)[1:60], which(xy$col >= 7 & !(xy$col >= 11 & xy$row >= 5))[1:60],
                 which(xy$col >= 11 & xy$row >= 5)[1:20])
        pos <- data.frame(barcode = grp$barcode, in_tissue = 1, row = xy$row[ord], col = xy$col[ord],
                          pxl_row = 100 * xy$row[ord], pxl_col = 100 * xy$col[ord])
        a$positions <- file.path(a$out, "provenance", "example_tissue_positions_list.csv")
        utils::write.table(pos, a$positions, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
        a$context <- "Simulated spatial ATAC spots built from the package's example peaks"
    }
}
skill_start("range-singlecell", a, a$out)
resolution <- if (!is.null(a$positions)) "spatial" else "single-cell"

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP", "give --txdb <TxDb package>")
txdb <- skill_load(a$txdb)
skill_check_genome(a$genome, txdb)
if (!is.null(a$fragments)) {
    if (is.null(a$groups)) skill_check("input", "groups given", "STOP", "give --groups <csv: barcode, group>")
    skill_input(c(a$fragments, a$groups), c("fragments", "groups"))
    fr <- readFragments(a$fragments)
    skill_check_recorded_genome(fr, a$genome, what = "fragments")
    fr <- skill_check_seqlevels(fr, txdb, fix = a$fix_seqlevels, what = "fragments")
    groups <- utils::read.csv(a$groups, colClasses = "character")
    tile <- as.numeric(a$tile)
    if (grepl("H3K27me3|H3K9me3|H3K36me3|H3K79me2|broad", a$assay, ignore.case = TRUE) && tile < 5000) {
        skill_check("input", "tile width for a broad mark", "NOTE",
                    paste0("tile ", tile, " bp; ?fragmentsToGroups suggests e.g. 5000 for a broad histone mark"))
    }
    rs <- skill_capture(fragmentsToGroups(fr, groups, tile = tile, min_fragments = as.numeric(a$min_fragments),
                                          specific = TRUE, fold = as.numeric(a$fold)), "fragmentsToGroups")
    m <- S4Vectors::metadata(rs)$match
    skill_check("input", "barcodes of the group table found in the fragments",
                if (m$matched / m$total < 0.9) "WARN" else "PASS",
                sprintf("%d of %d (%s); the adapter evaluation asked for >= 90%%", m$matched, m$total, m$how))
    members <- S4Vectors::metadata(rs)$members
} else if (!is.null(a$peaks_table)) {
    skill_input(a$peaks_table, "peaks_table")
    df <- utils::read.delim(a$peaks_table, check.names = FALSE)
    rs <- skill_capture(peaksToGroups(df), "peaksToGroups")
    members <- NULL
    groups <- NULL
    skill_check("input", "cells per group", "NOTE", "not known from a peak table; give them to judge the groups")
} else {
    skill_check("input", "fragments or peak table given", "STOP", "give --fragments + --groups, or --peaks_table")
}
n_reg <- lengths(rs)
for (g in names(rs)) {
    small <- (!is.null(members) && members[[g]] < 50) || n_reg[[g]] < 200
    skill_check("input", paste("size of group", g), if (small) "NOTE" else "PASS",
                paste0(if (!is.null(members)) paste0(members[[g]], " cells/spots, ") else "", n_reg[[g]], " regions",
                       if (small) "; small: fewer than 50 members (adapter evaluation rule) or 200 regions (build check unreliable)" else ""))
}
coords <- NULL
if (!is.null(a$positions)) {
    if (is.null(groups)) skill_check("input", "groups for the spots", "STOP", "--positions needs --fragments and --groups")
    skill_input(a$positions, "positions")
    pos <- readSpotPositions(a$positions)
    tis <- pos[pos$in_tissue == 1, ]
    gmap <- stats::setNames(groups$group, groups$barcode)
    tis$group <- unname(gmap[tis$barcode])
    found <- mean(!is.na(tis$group))
    skill_check("input", "in-tissue spots with a group", if (found < 0.9) "WARN" else "PASS",
                sprintf("%.0f%% of %d", 100 * found, nrow(tis)))
    coords <- tis[!is.na(tis$group), c("barcode", "x", "y", "group")]
}

## ---- 2. analysis with epiSeeker functions ----------------------------------------------
x <- skill_capture(annotateGroups(rs, TxDb = txdb, annoDb = a$annodb, resolution = resolution,
                                  members = members, coords = coords), "annotateGroups")
for (g in names(x)) skill_check_anno(x[[g]], paste("group", g), allow_suspect = a$allow_suspect)
gi <- groupInfo(x)
skill_table(gi, "groups.tsv")
stat <- do.call(rbind, lapply(names(x), function(g) cbind(group = g, getAnnoStat(x[[g]]))))
skill_table(stat, "annotation_stat.tsv")
skill_figure(plotAnnoBar(x), "annotation_bar.pdf")
skill_figure(plotDistToTSS(x), "distance_to_tss.pdf")
if (!is.null(coords)) skill_figure(plotSpatialGroups(coords), "spatial_groups.pdf")

## ---- 3. report ---------------------------------------------------------------------------
prom <- vapply(names(x), function(g) { st <- getAnnoStat(x[[g]]); sum(st$Frequency[grepl("^Promoter", st$Feature)]) }, 0)
sec <- list(
    "Groups" = c("| group | resolution | members | regions | % promoter |", "|---|---|---|---|---|",
                 paste0("| ", gi$group, " | ", gi$resolution, " | ", gi$n_members, " | ", gi$n_regions, " | ",
                        sprintf("%.1f", prom[gi$group]), " |")),
    "Note" = paste("Each group is a pseudo-bulk region set (tiles enriched in the group, or the group's peaks);",
                   "differences between groups are differences between these sets, not between single cells.")
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- skill_interpret(x, model = a$model, context = a$context, assay = a$assay, question = a$question)
skill_flag_phrases(res, c(
    "single-cell claim from a pseudo-bulk" = "\\b(every|each|all) (single )?(cells?|spots?)\\b",
    "causal claim from accessibility" = "\\b(regulates|drives|controls|activates|represses)\\b"))
skill_report("range-singlecell report", sec)
saveRDS(list(groups = rs, anno = x), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-singlecell done: ", skill_path("report", "report.md"))
