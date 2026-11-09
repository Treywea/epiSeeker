## range-spatial: spatial domains (spatial ATAC / CUT&Tag spots with a barcode -> domain table)
## -> input checks -> artefact check (do sequencing depth, chip rows or columns explain the
## domains?) -> region set per domain via fragmentsToGroups() -> annotateGroups(), comparison ->
## correctness checks -> interpret() last -> provenance. See ../SKILL.md. Shared steps:
## ../../shared/group-steps.R. The artefact check is the one of the paper's Fig 2 spatial panel.
##
## Rscript run_range_spatial.R --example --out <dir> [--model provider:model]
## Rscript run_range_spatial.R --fragments frags.tsv.gz --groups barcode_domain.csv \
##     --positions tissue_positions_list.csv --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
##     [--annodb org.Hs.eg.db] [--assay ATAC] [--tile 500] [--min_fragments 5] [--fold 2] \
##     [--spot_scores barcode_scores.csv] [--allow_depth_confound] [--context "..."] [--question "..."] \
##     [--fix_seqlevels] [--allow_suspect] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))
source(system.file("skills", "shared", "group-steps.R", package = "epiSeeker"))

a <- skill_args(list(tile = "500", min_fragments = "5", fold = "2", assay = "ATAC"))
if (isTRUE(a$example)) {
    ## the simulated cells of range-singlecell as spots on a 14 x 10 grid, depth varying per spot
    ## independently of the domain; domains concentric: AR_0M_like the 20 spots nearest the
    ## centre, AR_like the next 60, PcG_like the outer 60 (group_simulate_cells())
    sim <- group_simulate_cells(file.path(a$out, "provenance"), vary_depth = TRUE)
    a$fragments <- sim$fragments
    a$groups <- sim$groups
    grp <- utils::read.csv(sim$groups, colClasses = "character")
    xy <- expand.grid(col = 0:13, row = 0:9)
    ord <- order((xy$col - 6.5)^2 + (xy$row - 4.5)^2, xy$row, xy$col)
    slot <- c(ord[1:20], ord[21:80], ord[81:140])                      # centre, ring, outside
    who <- c(which(grp$group == "AR_0M_like"), which(grp$group == "AR_like"), which(grp$group == "PcG_like"))
    pos <- data.frame(barcode = grp$barcode[who], in_tissue = 1, row = xy$row[slot], col = xy$col[slot],
                      pxl_row = 100 * xy$row[slot], pxl_col = 100 * xy$col[slot])
    a$positions <- file.path(a$out, "provenance", "example_tissue_positions_list.csv")
    utils::write.table(pos, a$positions, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$context <- "Simulated spatial ATAC spots built from the package's example peaks"
}
skill_start("range-spatial", a, a$out)

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$txdb)) skill_check("input", "TxDb named", "STOP", "give --txdb <TxDb package>")
if (is.null(a$positions) || is.null(a$fragments)) {
    skill_check("input", "fragments, groups and spot positions given", "STOP",
                "give --fragments, --groups and --positions (the depth / row / column check needs all three)")
}
txdb <- skill_load(a$txdb)
skill_check_genome(a$genome, txdb)
gin <- group_inputs(a, txdb)
group_size_checks(gin$regions, gin$members)
skill_input(a$positions, "positions")
pos <- readSpotPositions(a$positions)
tis <- pos[pos$in_tissue == 1, ]
## barcodes of the fragments carry a "-1"-like suffix more often than the tables do
strip <- function(b) sub("-[0-9]+$", "", b)
key <- function(b, ref) if (mean(b %in% ref) >= 0.5) b else strip(b)
gmap <- stats::setNames(gin$groups$group, gin$groups$barcode)
tis$group <- unname(gmap[key(tis$barcode, names(gmap))])
found <- mean(!is.na(tis$group))
skill_check("input", "in-tissue spots with a group", if (found < 0.9) "WARN" else "PASS",
            sprintf("%.0f%% of %d", 100 * found, nrow(tis)))

## ---- artefact check: domains vs sequencing depth, chip rows, chip columns ---------------
## fragments per spot; the spots of the check are in tissue, with a group and >= 1 fragment
dep <- table(gin$fragments$barcode)
dk <- names(dep)
if (mean(dk %in% tis$barcode) < 0.5) dk <- strip(dk)
depth <- tapply(as.numeric(dep), dk, sum)
tis$depth <- unname(depth[tis$barcode])
sp <- tis[!is.na(tis$group) & !is.na(tis$depth) & tis$depth > 0, ]
skill_table(sp[, c("barcode", "row", "col", "group", "depth")], "spots.tsv")
coords <- sp[, c("barcode", "x", "y", "group")]
skill_figure(plotSpatialGroups(coords), "spatial_groups.pdf")
skill_figure(plotSpatialGroups(coords, value = stats::setNames(log10(sp$depth), sp$barcode)), "spatial_depth_log10.pdf")

## cross-validated balanced accuracy of LDA (5 folds x 5 repeats, folds stratified by domain),
## exactly as the Fig 2 check
bal_acc <- function(y, p) mean(tapply(p == y, y, mean))
cv_lda <- function(X, y, reps = 5, kf = 5) {
    y <- factor(y); X <- as.matrix(X)
    mean(vapply(seq_len(reps), function(r) {
        set.seed(r); folds <- integer(length(y))
        for (lv in levels(y)) { i <- which(y == lv); folds[i] <- sample(rep_len(seq_len(kf), length(i))) }
        pred <- factor(rep(NA, length(y)), levels = levels(y))
        for (f in seq_len(kf)) {
            tr <- folds != f
            v <- apply(X[tr, , drop = FALSE], 2, stats::var) > 0
            fit <- suppressWarnings(MASS::lda(X[tr, v, drop = FALSE], y[tr], tol = 1e-8))
            pred[!tr] <- stats::predict(fit, X[!tr, v, drop = FALSE])$class
        }
        bal_acc(y, pred)
    }, 0))
}
onehot <- function(f) { m <- stats::model.matrix(~ factor(f)); m[, -1, drop = FALSE] }
nspot <- table(sp$group)
cvd <- sp[sp$group %in% names(nspot)[nspot >= 5], ]
if (any(nspot < 5)) {
    skill_check("artefact", "domains with >= 5 spots for the cross-validation", "NOTE",
                paste("left out of the check:", paste(names(nspot)[nspot < 5], collapse = ", ")))
}
if (length(unique(cvd$group)) < 2) skill_check("artefact", "two or more domains to check", "STOP", "fewer than 2 domains with >= 5 spots")
k <- length(unique(cvd$group)); chance <- 1 / k
kap <- function(ba) (ba - chance) / (1 - chance)
ba <- c(depth = cv_lda(log10(cvd$depth), cvd$group), row = cv_lda(onehot(cvd$row), cvd$group),
        column = cv_lda(onehot(cvd$col), cvd$group))
if (!is.null(a$spot_scores)) {
    skill_input(a$spot_scores, "spot_scores")
    sc <- utils::read.csv(a$spot_scores, check.names = FALSE)
    m <- match(key(cvd$barcode, sc[[1]]), sc[[1]])
    if (anyNA(m)) skill_check("artefact", "spot scores for every spot of the check", "STOP", paste(sum(is.na(m)), "spots without scores"))
    ba <- c(ba, biology = cv_lda(as.matrix(sc[m, -1, drop = FALSE]), cvd$group))
}
r2 <- summary(stats::lm(log10(cvd$depth) ~ factor(cvd$group)))$r.squared
art <- data.frame(predictor = c(names(ba), "R2: log10 depth ~ domain"), balanced_accuracy = c(round(ba, 4), NA),
                  chance = c(rep(round(chance, 4), length(ba)), NA), kappa = c(round(kap(ba), 4), round(r2, 4)),
                  spots = nrow(cvd), domains = k)
skill_table(art, "artefact_check.tsv")
skill_figure(ggplot2::ggplot(art[!is.na(art$balanced_accuracy), ], ggplot2::aes(x = .data$predictor, y = .data$balanced_accuracy)) +
                 ggplot2::geom_col(fill = "grey45") + ggplot2::geom_hline(yintercept = chance, linetype = 2) +
                 ggplot2::coord_flip() + ggplot2::theme_minimal() +
                 ggplot2::labs(x = NULL, y = "CV balanced accuracy of LDA predicting the domain (dashed: chance)"),
             "artefact_check.pdf", width = 6, height = 3)

kd <- kap(ba[["depth"]])
dtxt <- sprintf("depth predicts the domain with balanced accuracy %.3f (chance %.3f, kappa %.2f); R^2 of log10 depth ~ domain %.2f",
                ba[["depth"]], chance, kd, r2)
if (kd >= 0.4 || r2 >= 0.4) {
    skill_check("artefact", "domains explained by sequencing depth", if (isTRUE(a$allow_depth_confound)) "WARN" else "STOP",
                paste0(dtxt, "; thresholds kappa 0.4 / R^2 0.4. The domains follow depth: they cannot be read as anatomy",
                       if (isTRUE(a$allow_depth_confound)) " (--allow_depth_confound: run continued)" else
                           ". Re-cluster with depth removed (e.g. drop depth-correlated components, downsample), or --allow_depth_confound to go on"))
    depth_status <- "confounded"
} else if (kd >= 0.15 || r2 >= 0.15) {
    skill_check("artefact", "domains explained by sequencing depth", "WARN",
                paste0(dtxt, "; partly entangled with depth (thresholds 0.15): read domain differences with care"))
    depth_status <- "partly entangled"
} else {
    skill_check("artefact", "domains explained by sequencing depth", "PASS", dtxt)
    depth_status <- "not explained by depth"
}
for (ax in c("row", "column")) {
    ka <- kap(ba[[ax]])
    skill_check("artefact", paste0("domains follow chip ", ax, "s"), if (ka >= 0.4) "WARN" else "PASS",
                sprintf("balanced accuracy %.3f from the %s alone (kappa %.2f)%s", ba[[ax]], ax, ka,
                        if (ka >= 0.4) paste0("; a channel artefact (DBiT-type chips) or banded anatomy: compare with the histology") else ""))
}
stripe <- vapply(split(sp, sp$group), function(d) if (nrow(d) < 10) NA_real_ else max(max(table(d$row)), max(table(d$col))) / nrow(d), 0)
st_bad <- names(stripe)[!is.na(stripe) & stripe >= 0.5]
skill_check("artefact", "domains confined to one row or column (stripes)", if (length(st_bad)) "WARN" else "PASS",
            if (length(st_bad)) paste0(paste0(st_bad, " (", round(100 * stripe[st_bad]), "% of spots in one row/column)"), collapse = "; ")
            else sprintf("largest share of a domain in one row or column: %.0f%%", 100 * max(stripe, na.rm = TRUE)))
if ("biology" %in% names(ba)) {
    best <- max(ba[c("depth", "row", "column")])
    skill_check("artefact", "spot scores beat depth, rows and columns by >= 0.05", if (ba[["biology"]] >= best + 0.05) "PASS" else "WARN",
                sprintf("biology %.3f vs best of depth/row/column %.3f (Fig 2 rule SP3)", ba[["biology"]], best))
}

## ---- 2. analysis with epiSeeker functions (correctness checks inside) ----------------------
x <- group_annotate(gin$regions, txdb, a$annodb, "spatial", gin$members, coords = coords, allow_suspect = a$allow_suspect)

## ---- 3. report ---------------------------------------------------------------------------
sec <- list(
    "Domains" = group_section(x),
    "Artefact check" = c("| predictor | balanced accuracy | chance | kappa (R^2 for the last row) |", "|---|---|---|---|",
                         paste0("| ", art$predictor, " | ", art$balanced_accuracy, " | ", art$chance, " | ", art$kappa, " |"),
                         "", paste0("Depth: domains ", depth_status, ". ", nrow(cvd), " spots, ", k, " domains; tables/artefact_check.tsv, ",
                                    "figures/spatial_groups.pdf, spatial_depth_log10.pdf, artefact_check.pdf.")),
    "Note" = paste("Each domain is a pseudo-bulk region set of its spots. Rows and columns are those of the chip;",
                   if (depth_status == "confounded") "the domains follow sequencing depth, so no anatomical reading is supported." else
                       "a domain is an anatomical structure only if the histology or independent markers say so.")
)

## ---- 4. interpret() last -----------------------------------------------------------------
ctx <- paste(c(if (!is.null(a$context)) sub("([^.])$", "\\1.", a$context),
               paste0("Artefact check: ", dtxt, " (domains ", depth_status, ").")), collapse = " ")
res <- skill_interpret(x, model = a$model, context = ctx, assay = a$assay, question = a$question)
skill_flag_phrases(res, c(
    "single-spot claim from a pseudo-bulk" = "\\b(every|each|all) (single )?(cells?|spots?)\\b",
    "causal claim from accessibility" = "\\b(regulates|drives|controls|activates|represses)\\b"))
skill_report("range-spatial report", sec)
saveRDS(list(groups = gin$regions, anno = x, artefact = art, spots = sp), skill_path("provenance", "objects.rds"))
skill_finish()
message("range-spatial done: ", skill_path("report", "report.md"))
