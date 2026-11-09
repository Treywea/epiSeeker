## Steps shared by the skills for groups of cells or spots (range-singlecell, range-spatial).
## Source skill-helpers.R first. They order the epiSeeker adapter calls and record the checks;
## the regions of each group come from fragmentsToGroups() / peaksToGroups().

## ---- the simulated example of both skills -------------------------------------------------

## Toy fragments simulated from the bundled peak sets (seed 1): cells of group "AR_like" put 70%
## of their fragments in AR peaks (ARmo_1nM), "PcG_like" in CBX7 peaks, the small group
## "AR_0M_like" (20 cells) in ARmo_0M peaks; the rest fall anywhere on chr1-22. Fragment
## barcodes carry "-1", the group table does not. 300 fragments per cell, or, with
## vary_depth = TRUE, 100-300 of them kept per cell independently of its group (log-uniform),
## so that depth varies between spots as it does in real data. It shows the mechanics; not real data.
group_simulate_cells <- function(dir, vary_depth = FALSE) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
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
    if (vary_depth) {
        ## keep a random share of each cell's fragments: depth independent of the group
        set.seed(2)
        n_keep <- round(exp(stats::runif(k, log(100), log(300))))
        fr <- lapply(seq_len(k), function(i) fr[[i]][sort(sample(nrow(fr[[i]]), n_keep[i])), ])
    }
    fr <- do.call(rbind, fr); fr <- fr[order(fr$chr, fr$start), ]
    out <- list(fragments = file.path(dir, "example_fragments.tsv.gz"), groups = file.path(dir, "example_groups.csv"))
    con <- gzfile(out$fragments, "w"); utils::write.table(fr, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(con)
    utils::write.csv(do.call(rbind, grp), out$groups, row.names = FALSE, quote = FALSE)
    out
}

## ---- inputs --------------------------------------------------------------------------------

## fragments + barcode -> group table -> fragmentsToGroups(), with the barcode check; or a
## marker-peak table -> peaksToGroups(). Returns the regions per group, members, the table.
group_inputs <- function(a, txdb) {
    if (!is.null(a$fragments)) {
        if (is.null(a$groups)) skill_check("input", "groups given", "STOP", "give --groups <csv: barcode, group>")
        skill_input(c(a$fragments, a$groups), c("fragments", "groups"))
        fr <- skill_capture(readFragments(a$fragments), "readFragments")
        skill_check_recorded_genome(fr, a$genome, what = "fragments")
        fr <- skill_check_seqlevels(fr, txdb, fix = a$fix_seqlevels, what = "fragments")
        groups <- utils::read.csv(a$groups, colClasses = "character")
        ## the column names fragmentsToGroups() accepts, renamed to barcode / group
        bcol <- intersect(c("barcode", "cell", "spot"), names(groups))[1]
        gcol <- intersect(c("group", "cluster", "domain", "cell_type", "celltype"), names(groups))[1]
        if (is.na(bcol) || is.na(gcol)) {
            skill_check("input", "group table has a barcode and a group column", "STOP",
                        "columns barcode|cell|spot and group|cluster|domain|cell_type")
        }
        groups <- data.frame(barcode = groups[[bcol]], group = groups[[gcol]], stringsAsFactors = FALSE)
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
        return(list(regions = rs, members = S4Vectors::metadata(rs)$members, groups = groups, fragments = fr))
    }
    if (!is.null(a$peaks_table)) {
        skill_input(a$peaks_table, "peaks_table")
        df <- utils::read.delim(a$peaks_table, check.names = FALSE)
        rs <- skill_capture(peaksToGroups(df), "peaksToGroups")
        skill_check("input", "cells per group", "NOTE", "not known from a peak table; give them to judge the groups")
        return(list(regions = rs, members = NULL, groups = NULL, fragments = NULL))
    }
    skill_check("input", "fragments or peak table given", "STOP", "give --fragments + --groups, or --peaks_table")
}

## a group of few cells is a noisy pseudo-bulk; below 200 regions the build check is unreliable
group_size_checks <- function(rs, members) {
    n_reg <- lengths(rs)
    for (g in names(rs)) {
        small <- (!is.null(members) && members[[g]] < 50) || n_reg[[g]] < 200
        skill_check("input", paste("size of group", g), if (small) "NOTE" else "PASS",
                    paste0(if (!is.null(members)) paste0(members[[g]], " cells/spots, ") else "", n_reg[[g]], " regions",
                           if (small) "; small: fewer than 50 members (adapter evaluation rule) or 200 regions (build check unreliable)" else ""))
    }
}

## ---- annotation ------------------------------------------------------------------------------

## annotateGroups() -> csAnnoList; annoChecks() of every group; group table, feature table and
## the two standard figures
group_annotate <- function(rs, txdb, annodb, resolution, members, coords = NULL, allow_suspect = FALSE) {
    x <- skill_capture(annotateGroups(rs, TxDb = txdb, annoDb = annodb, resolution = resolution,
                                      members = members, coords = coords), "annotateGroups")
    for (g in names(x)) skill_check_anno(x[[g]], paste("group", g), allow_suspect = allow_suspect)
    skill_table(groupInfo(x), "groups.tsv")
    stat <- do.call(rbind, lapply(names(x), function(g) cbind(group = g, getAnnoStat(x[[g]]))))
    skill_table(stat, "annotation_stat.tsv")
    skill_figure(plotAnnoBar(x), "annotation_bar.pdf")
    skill_figure(plotDistToTSS(x), "distance_to_tss.pdf")
    x
}

## the Groups section of the report
group_section <- function(x) {
    gi <- groupInfo(x)
    prom <- vapply(names(x), function(g) { st <- getAnnoStat(x[[g]]); sum(st$Frequency[grepl("^Promoter", st$Feature)]) }, 0)
    c("| group | resolution | members | regions | % promoter |", "|---|---|---|---|---|",
      paste0("| ", gi$group, " | ", gi$resolution, " | ", gi$n_members, " | ", gi$n_regions, " | ",
             sprintf("%.1f", prom[gi$group]), " |"))
}
