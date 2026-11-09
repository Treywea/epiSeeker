## base-singlemolecule: per-read tables of fibertools (Fiber-seq: m6A, nucleosomes, MSPs) or modkit
## (5mC / 6mA calls of any modBAM) -> input checks -> reads, per-element openness, pairwise
## co-occurrence of elements -> the per-read openness confound (a per-read covariate, its link to
## element openness, co-occurrence re-fitted with it) -> correctness checks -> interpret() last ->
## provenance. See ../SKILL.md. Raw co-occurrence is confounded by how open each read is overall
## (checks/20261107_adapter-single-molecule in the epiSeeker-LLM records): read the adjusted column.
##
## Rscript run_base_singlemolecule.R --example --out <dir> [--model provider:model]
## Rscript run_base_singlemolecule.R --example_control --out <dir>   (equal labelling efficiency)
## Rscript run_base_singlemolecule.R --input reads.ft_all.tsv.gz --genome hg38 --elements peaks.bed \
##     [--format auto|fibertools|modkit] [--type msp] [--summary_type m6A] [--cutoff 0.5] \
##     [--min_dist 1000] [--max_dist 10000] [--region chr1:1-2] [--txdb TxDb...] [--annodb org...] \
##     [--assay "Fiber-seq"] [--context "..."] [--question "..."] --out <dir> [--model provider:model]

source(system.file("skills", "shared", "skill-helpers.R", package = "epiSeeker"))

a <- skill_args(list(format = "auto", cutoff = "0.5", min_dist = "1000", max_dist = "10000"))
if (isTRUE(a$example) || isTRUE(a$example_control)) {
    ## A fibertools `ft extract --all` table simulated over chr1:1,410,000-1,440,000 (seed 1):
    ## 200 reads of 12-25 kb; the elements are the LNCaP DNase peaks of the window
    ## (inst/extdata/encode). Each read has a labelling efficiency z ~ N(0, 1): an element is open
    ## (covered by an MSP) with probability plogis(logit(base) + 1.5 z), and linker MSPs outside
    ## the elements occur at a rate exp(0.8 z) -- so, given the read, elements open independently,
    ## but reads labelled well look open everywhere. --example_control sets z = 0 for every read.
    ## Nucleosomes fill the rest; m6A calls sit on a fixed grid of A positions inside MSPs.
    ## It shows the per-read path and the confound; it is not real Fiber-seq data.
    dir.create(file.path(a$out, "provenance"), recursive = TRUE, showWarnings = FALSE)
    w0 <- 1410000; w1 <- 1440000
    dhs <- readPeakFile(system.file("extdata", "encode", "LNCaP_R1881_DNase_ENCFF570MWX_chr1-2.bed.gz", package = "epiSeeker"))
    el <- dhs[as.character(GenomicRanges::seqnames(dhs)) == "chr1" & GenomicRanges::start(dhs) >= w0 & GenomicRanges::end(dhs) <= w1]
    el <- GenomicRanges::reduce(GenomicRanges::granges(el), min.gapwidth = 101L)
    es <- GenomicRanges::start(el); ee <- GenomicRanges::end(el)
    a$elements <- file.path(a$out, "provenance", "example_elements.bed")
    utils::write.table(data.frame("chr1", es - 1, ee), a$elements, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    set.seed(1)
    base <- stats::runif(length(el), 0.2, 0.6)
    agrid <- sort(sample(w0:w1, (w1 - w0) %/% 8))
    nr <- 200
    z <- if (isTRUE(a$example_control)) rep(0, nr) else stats::rnorm(nr)
    rows <- vector("list", nr)
    for (r in seq_len(nr)) {
        len <- sample(12000:25000, 1)
        st <- sample((w0 - 3000):(w1 + 3000 - len), 1); en <- st + len
        inr <- which(es >= st & ee <= en)
        op <- inr[stats::runif(length(inr)) < stats::plogis(stats::qlogis(base[inr]) + 1.5 * z[r])]
        ms <- es[op] - 10; ml <- ee[op] - es[op] + 21                      # MSPs over open elements
        nlk <- stats::rpois(1, len / 1000 * 1.5 * exp(0.8 * z[r]))        # linker MSPs elsewhere
        if (nlk) {
            ls <- sort(sample(st:(en - 200), nlk))
            ll <- sample(60:150, nlk, replace = TRUE)
            ok <- !IRanges::overlapsAny(IRanges::IRanges(ls, width = ll), IRanges::IRanges(es - 200, ee + 200))
            ms <- c(ms, ls[ok]); ml <- c(ml, ll[ok])
        }
        o <- order(ms); ms <- ms[o]; ml <- ml[o]
        keep <- c(TRUE, ms[-1] > (ms + ml)[-length(ms)] + 20)              # no overlapping MSPs
        ms <- ms[keep]; ml <- ml[keep]
        ## nucleosomes (147 bp, 20-60 bp linkers) in the gaps between MSPs
        gs <- c(st, ms + ml); ge <- c(ms, en)
        ns <- unlist(lapply(seq_along(gs), function(i) if (ge[i] - gs[i] > 200) seq(gs[i] + 20, ge[i] - 170, by = 147 + sample(20:60, 1))))
        m6 <- agrid[agrid >= st & agrid <= en & IRanges::overlapsAny(IRanges::IRanges(agrid, width = 1)[agrid >= st & agrid <= en],
                                                                     IRanges::IRanges(ms, width = ml))]
        m6 <- m6[stats::runif(length(m6)) < 0.7]
        cs <- function(v) if (length(v)) paste(v, collapse = ",") else "."
        rows[[r]] <- data.frame(ct = "chr1", st = st - 1, en = en, fiber = sprintf("fiber%03d", r), strand = "+",
                                ref_m6a = cs(m6 - 1), m6a_qual = cs(sample(200:255, length(m6), replace = TRUE)),
                                ref_nuc_starts = cs(ns - 1), ref_nuc_lengths = cs(rep(147, length(ns))),
                                ref_msp_starts = cs(ms - 1), ref_msp_lengths = cs(ml), check.names = FALSE)
    }
    tab <- do.call(rbind, rows)
    names(tab)[1] <- "#ct"
    a$input <- file.path(a$out, "provenance", if (isTRUE(a$example_control)) "example_control.ft_all.tsv.gz" else "example.ft_all.tsv.gz")
    con <- gzfile(a$input, "w"); utils::write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
    a$format <- "fibertools"
    a$type <- "msp"
    a$genome <- "hg38"
    a$txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    a$annodb <- "org.Hs.eg.db"
    a$assay <- "Fiber-seq m6A (simulated)"
    a$context <- paste("Simulated Fiber-seq reads over LNCaP DNase peaks on chr1",
                       if (isTRUE(a$example_control)) "(every read labelled equally)" else "(labelling efficiency varies between reads)")
}
skill_start("base-singlemolecule", a, a$out)
cutoff <- as.numeric(a$cutoff)

## ---- 1. input checks -------------------------------------------------------------------
if (is.null(a$input)) skill_check("input", "per-read table given", "STOP", "give --input <ft extract --all / modkit extract table>")
if (is.null(a$genome) || isTRUE(a$genome)) skill_check("input", "genome build declared", "STOP", "give --genome (e.g. hg38)")
txdb <- skill_load(a$txdb)
if (!is.null(txdb)) skill_check_genome(a$genome, txdb)
skill_input(a$input, "reads")
reads <- skill_capture(readModCalls(a$input, format = a$format), "readModCalls")
src <- S4Vectors::metadata(reads)$modSource
rd <- S4Vectors::metadata(reads)$reads
skill_check("input", "per-read calls placed on the reference", if (src$n_unplaced > 0) "NOTE" else "PASS",
            paste0(nrow(rd), " reads; ", src$n_unplaced, " calls without a reference position dropped; format ", src$format))
types <- unique(reads$type)
if (is.null(a$type)) a$type <- if ("msp" %in% types) "msp" else if ("5mC" %in% types) "5mC" else types[1]
skill_check("input", paste("call type", a$type, "present"), if (a$type %in% types) "PASS" else "STOP",
            paste("types:", paste(types, collapse = ", ")))
skill_check_recorded_genome(reads, a$genome, what = "reads")
if (!is.null(txdb)) skill_check_seqlevels(GenomicRanges::GRanges(unique(rd$chr), IRanges::IRanges(1, 1)), txdb, what = "reads")
el <- NULL
if (!is.null(a$elements)) {
    skill_input(a$elements, "elements")
    el <- GenomicRanges::granges(readPeakFile(a$elements))
    st <- GenomeInfoDb::seqlevelsStyle(el); rs <- GenomeInfoDb::seqlevelsStyle(reads)
    skill_check("input", "seqlevels style of elements matches the reads", if (length(intersect(st, rs))) "PASS" else "STOP",
                paste0(paste(st, collapse = "/"), " vs reads ", paste(rs, collapse = "/")))
    rgr <- GenomicRanges::GRanges(rd$chr, IRanges::IRanges(rd$start, rd$end))
    cov_el <- IRanges::overlapsAny(el, rgr, type = "within")
    skill_check("input", "elements spanned by at least one read", if (mean(cov_el) < 0.5) "WARN" else "PASS",
                sprintf("%d of %d", sum(cov_el), length(el)))
    el <- el[cov_el]
    el <- el[order(as.character(GenomicRanges::seqnames(el)), GenomicRanges::start(el))]
    names(el) <- paste0(GenomicRanges::seqnames(el), ":", GenomicRanges::start(el), "-", GenomicRanges::end(el))
}

## ---- 2. analysis with epiSeeker functions ----------------------------------------------
st_type <- if (!is.null(a$summary_type)) a$summary_type else if ("5mC" %in% types) "5mC" else if ("m6A" %in% types) "m6A" else NULL
bm <- NULL
if (!is.null(st_type)) {
    bm <- skill_capture(summarizeModCalls(reads, type = st_type, sampleName = "reads"), "summarizeModCalls")
    skill_check("input", "how coverage of per-read calls is counted", "NOTE", if (identical(src$format, "fibertools"))
        "fibertools lists only calls above its threshold: coverage = reads spanning the site (see ?summarizeModCalls)" else
        "the table lists every tested base: coverage = calls at the site (see ?summarizeModCalls)")
}
anno <- NULL; open_tab <- NULL; co <- NULL; rho_med <- NA; read_cov <- NULL
if (!is.null(el) && length(el)) {
    if (!is.null(txdb)) {
        anno <- skill_capture(annotateSeq(el, TxDb = txdb, annoDb = a$annodb, verbose = FALSE), "annotateSeq elements")
        skill_check_anno(anno, "elements", allow_suspect = a$allow_suspect)
        skill_table(as.data.frame(anno), "elements_annotation.tsv")
    }
    m <- skill_capture(modReadMatrix(reads, el, type = a$type), "modReadMatrix")
    open_tab <- data.frame(element = names(el), reads_spanning = colSums(!is.na(m)),
                           fraction_open = round(colMeans(m >= cutoff, na.rm = TRUE), 4), row.names = NULL)
    skill_table(open_tab, "element_openness.tsv")

    ## per-read covariate: how open / modified the read is outside the elements
    elr <- GenomicRanges::reduce(el)
    bp_out <- function(gr) {
        h <- GenomicRanges::findOverlaps(gr, elr, ignore.strand = TRUE)
        w_in <- GenomicRanges::width(GenomicRanges::pintersect(gr[S4Vectors::queryHits(h)], elr[S4Vectors::subjectHits(h)], ignore.strand = TRUE))
        bin <- tapply(w_in, factor(S4Vectors::queryHits(h), levels = seq_along(gr)), sum)
        bin[is.na(bin)] <- 0
        GenomicRanges::width(gr) - as.numeric(bin)
    }
    rgr <- GenomicRanges::GRanges(rd$chr, IRanges::IRanges(rd$start, rd$end))
    rout <- stats::setNames(bp_out(rgr), rd$read_id)
    y <- reads[reads$type == a$type]
    if (a$type %in% c("msp", "nuc")) {
        yo <- tapply(bp_out(GenomicRanges::granges(y)), y$read_id, sum)
        cv <- as.numeric(yo[rd$read_id]) / rout
        cv_def <- paste0("fraction of the read outside the elements covered by ", a$type)
    } else {
        yo <- y[!IRanges::overlapsAny(y, elr, ignore.strand = TRUE)]
        if (identical(src$format, "modkit")) {
            cv <- as.numeric(tapply(yo$modified, yo$read_id, mean)[rd$read_id])
            cv_def <- paste0("fraction of ", a$type, " calls modified outside the elements")
        } else {
            cv <- as.numeric(tapply(yo$modified, yo$read_id, sum)[rd$read_id]) / rout * 1000
            cv_def <- paste0(a$type, " calls per kb outside the elements")
        }
    }
    cv[is.na(cv)] <- 0
    names(cv) <- rd$read_id
    read_cov <- data.frame(read_id = rd$read_id, chr = rd$chr, start = rd$start, end = rd$end, covariate = round(cv, 5))
    skill_table(read_cov, "reads.tsv")
    q <- stats::quantile(cv, c(0.25, 0.5, 0.75), names = FALSE)
    skill_note("per_read_covariate", cv_def)

    ## does a read's overall openness track whether it is open at each element?
    rho <- vapply(seq_along(el), function(j) {
        ok <- !is.na(m[, j]); o <- m[ok, j] >= cutoff
        if (sum(ok) < 10 || length(unique(o)) < 2) return(NA_real_)
        suppressWarnings(stats::cor(as.numeric(o), cv[rownames(m)[ok]], method = "spearman"))
    }, 0)
    rho_med <- stats::median(rho, na.rm = TRUE)
    skill_check("correctness", "per-read openness tracks element openness", if (!is.na(rho_med) && rho_med >= 0.1) "WARN" else "PASS",
                sprintf("covariate = %s (median %.3g, IQR %.3g-%.3g); median Spearman(open at element, covariate) %.2f over %d elements%s",
                        cv_def, q[2], q[1], q[3], rho_med, sum(!is.na(rho)),
                        if (!is.na(rho_med) && rho_med >= 0.1) "; reads open everywhere make every pair look co-open: raw co-occurrence is confounded, read the adjusted odds ratios" else ""))

    ## co-occurrence of element pairs min_dist-max_dist apart: raw (modCooccurrence) and adjusted
    ## for the per-read covariate (glm; modCooccurrence() has no such option yet)
    cen <- (GenomicRanges::start(el) + GenomicRanges::end(el)) / 2
    chr <- as.character(GenomicRanges::seqnames(el))
    pr <- which(outer(seq_along(el), seq_along(el), "<") & outer(chr, chr, "==") &
                    abs(outer(cen, cen, "-")) >= as.numeric(a$min_dist) & abs(outer(cen, cen, "-")) <= as.numeric(a$max_dist), arr.ind = TRUE)
    if (nrow(pr)) {
        co <- do.call(rbind, lapply(seq_len(nrow(pr)), function(k) {
            i <- pr[k, 1]; j <- pr[k, 2]
            r <- modCooccurrence(reads, el[i], el[j], type = a$type, cutoff = cutoff)
            ok <- stats::complete.cases(m[, c(i, j)])
            adj <- c(NA_real_, NA_real_)
            if (!is.na(r$p_value)) {
                d <- data.frame(o1 = m[ok, i] >= cutoff, o2 = m[ok, j] >= cutoff, cv = cv[rownames(m)[ok]])
                f <- tryCatch(suppressWarnings(stats::glm(o2 ~ o1 + cv, family = stats::binomial, data = d)), error = function(e) NULL)
                if (!is.null(f) && f$converged && "o1TRUE" %in% rownames(summary(f)$coefficients)) {
                    cf <- summary(f)$coefficients
                    adj <- c(exp(cf["o1TRUE", 1]), cf["o1TRUE", 4])
                }
            }
            data.frame(element1 = names(el)[i], element2 = names(el)[j], distance = round(abs(cen[j] - cen[i])),
                       n_reads = r$n_reads, odds_ratio_unadjusted = signif(r$odds_ratio, 4), p_unadjusted = signif(r$p_value, 4),
                       odds_ratio_adjusted = signif(adj[1], 4), p_adjusted = signif(adj[2], 4))
        }))
        skill_table(co, "cooccurrence.tsv")
        testable <- !is.na(co$p_unadjusted)
        pos_raw <- sum(co$odds_ratio_unadjusted > 1 & co$p_unadjusted < 0.05, na.rm = TRUE)
        pos_adj <- sum(co$odds_ratio_adjusted > 1 & co$p_adjusted < 0.05, na.rm = TRUE)
        n_adj <- sum(!is.na(co$p_adjusted))
        skill_check("correctness", "co-occurrence adjusted for per-read openness",
                    if (n_adj == sum(testable)) "PASS" else "WARN",
                    sprintf("%d pairs %s-%s bp apart, %d testable; OR > 1 and p < 0.05: %d unadjusted, %d of %d adjusted%s",
                            nrow(co), a$min_dist, a$max_dist, sum(testable), pos_raw, pos_adj, n_adj,
                            if (n_adj < sum(testable)) "; pairs without an adjusted fit stay confounded by per-read openness" else ""))
        skill_figure(ggplot2::ggplot(co[testable, ], ggplot2::aes(x = log2(.data$odds_ratio_unadjusted), y = log2(.data$odds_ratio_adjusted))) +
                         ggplot2::geom_abline(linetype = 2) + ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
                         ggplot2::geom_point(alpha = 0.6) + ggplot2::theme_minimal() +
                         ggplot2::labs(x = "log2 odds ratio, unadjusted (modCooccurrence)", y = "log2 odds ratio, adjusted for per-read openness"),
                     "cooccurrence_adjusted_vs_unadjusted.pdf", width = 5, height = 5)
    } else {
        skill_check("correctness", "element pairs for co-occurrence", "NOTE", paste0("none ", a$min_dist, "-", a$max_dist, " bp apart"))
    }
}
## reads over a region: --region, else the span of the elements on their first chromosome (<= 150 kb)
rg <- if (!is.null(a$region)) GenomicRanges::GRanges(a$region) else if (!is.null(el) && length(el)) {
    e1 <- el[as.character(GenomicRanges::seqnames(el)) == as.character(GenomicRanges::seqnames(el))[1]]
    GenomicRanges::GRanges(GenomicRanges::seqnames(e1)[1], IRanges::IRanges(min(GenomicRanges::start(e1)) - 1000,
                                                                             min(max(GenomicRanges::end(e1)), min(GenomicRanges::start(e1)) + 150000) + 1000))
} else NULL
if (!is.null(rg)) {
    skill_figure(plotModReads(reads, rg, types = intersect(c("nuc", "msp", "m6A", "5mC"), types), peaks = el),
                 "reads_region.pdf", width = 10, height = 7)
}

## ---- 3. report ---------------------------------------------------------------------------
sec <- list(
    "Input" = c(paste0("- ", nrow(rd), " reads (", src$format, "); call types ", paste(types, collapse = ", "),
                       "; genome ", a$genome),
                if (!is.null(el)) paste0("- ", length(el), " elements spanned by reads (", a$elements, ")"),
                if (!is.null(bm)) paste0("- ", st_type, " summarised to ", nrow(bm), " sites (summarizeModCalls)")),
    "Elements" = if (is.null(open_tab)) "no --elements given" else
        paste0("- open (", a$type, " >= ", cutoff, ") on ", round(100 * stats::median(open_tab$fraction_open), 1),
               "% of spanning reads at the median element; tables/element_openness.tsv"),
    "Co-occurrence" = if (is.null(co)) "not run" else
        c(paste0("- per-read covariate: ", .skill$notes$per_read_covariate, "; median Spearman with element openness ",
                 round(rho_med, 2)),
          paste0("- ", nrow(co), " pairs; OR > 1 and p < 0.05: ", sum(co$odds_ratio_unadjusted > 1 & co$p_unadjusted < 0.05, na.rm = TRUE),
                 " unadjusted, ", sum(co$odds_ratio_adjusted > 1 & co$p_adjusted < 0.05, na.rm = TRUE), " adjusted (tables/cooccurrence.tsv)"),
          "", paste("**Unadjusted co-occurrence is confounded by per-read openness**: a read that is labelled",
                    "or opened more everywhere makes every pair of elements look co-open. Only the adjusted odds",
                    "ratios (or a comparison with distance-matched control pairs) speak to elements opening together."))
)

## ---- 4. interpret() last -----------------------------------------------------------------
res <- NULL
if (!is.null(bm)) {
    res <- skill_interpret(bm, model = a$model, context = a$context, assay = a$assay, question = a$question,
                           TxDb = txdb, annoDb = a$annodb, min_cov = 5)
} else {
    skill_check("interpret", "interpret() run", "NOTE", "no base-level calls to summarise; interpret() not run")
}
skill_flag_phrases(res, c(
    "coordination from co-occurrence" = "\\b(coordinat(ed|ion|ely)|cooperat\\w*|co-regulat\\w*|in concert|act(s|ing)? together)\\b",
    "causal claim from openness" = "\\b(regulates|drives|controls|activates|represses)\\b"))
skill_report("base-singlemolecule report", sec)
saveRDS(list(bmData = bm, elements = el, openness = open_tab, cooccurrence = co, reads = read_cov), skill_path("provenance", "objects.rds"))
skill_finish()
message("base-singlemolecule done: ", skill_path("report", "report.md"))
