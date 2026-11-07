library(epiSeeker)

context("test interpret() for bmData and loopAnno, and the strict options of linkPeaksToGenes()")

## ---- fixtures ---------------------------------------------------------------------------------

## two samples, five sites in two clusters (100-300 and 2000-2100 on chr1), one site under
## coverage 5 in sample b
bm_fixture <- function() {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 200, 300, 2000, 2100), width = 1))
    lev <- cbind(a = c(0.9, 0.8, 1.0, 0.1, 0.0), b = c(0.2, 0.9, 0.4, 0.1, 0.6))
    cov <- cbind(a = c(10, 10, 10, 10, 10), b = c(10, 10, 2, 10, 10))
    suppressWarnings(bmData(value1 = cov, value2 = lev, gr = gr, sampleNames = c("a", "b"),
                            valueNames = c("Cov", "Methylation")))
}

## same genes and loops as test-loops.R, plus a non-coding gene G4 (TSS 150,000 +)
lp_genes <- function() {
    GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10000, 100000, 300000, 150000),
                                                    c(20000, 110000, 310000, 152000)),
                           strand = c("+", "+", "-", "+"), gene_id = c("G1", "G2", "G3", "G4"))
}
lp_bedpe <- function() {
    data.frame(c1 = c("chr1", "chr1", "chr1", "chr1", "chr1"),
               s1 = c(40000, 150000, 8000, 40000, 115000), e1 = c(45000, 155000, 12000, 45000, 120000),
               c2 = c("chr1", "chr1", "chr1", "chr2", "chr1"),
               s2 = c(98000, 200000, 308000, 1000, 98000), e2 = c(103000, 205000, 312000, 5000, 103000),
               score = c(5, 4, 3, 2, 1))
}
lp_anno <- function() {
    suppressWarnings(suppressMessages(annotateLoops(lp_bedpe(), TxDb = lp_genes(), verbose = FALSE)))
}

mock_schema_answer <- function(env) {
    local_mocked_bindings(
        generate_object = function(...) {
            a <- list(...); env$seen <- a
            list(object = lapply(a$schema$properties, epiSeeker:::.na_of), finish_reason = "stop",
                 raw_text = "{}")
        },
        .package = "aisdk", .env = parent.frame()
    )
}

## ---- bmData ------------------------------------------------------------------------------------

test_that("the bmData digest reports clusters, differences, levels and coverage exactly", {
    x <- bm_fixture()
    ev <- epiSeeker:::.bm_digest(x, assay = "WGBS")
    expect_equal(names(ev$sections), c("overview", "regions", "differences", "levels", "coverage",
                                       "genomic_context"))
    expect_true(all(nchar(ev$sha256) == 64))
    expect_match(ev$sections$overview, "Sites: 5 single bases in 2 clusters", fixed = TRUE)
    expect_match(ev$sections$overview, "Level assay: Methylation; coverage assay: Cov", fixed = TRUE)
    expect_match(ev$sections$regions, "Clusters (2 of 2, most sites first)", fixed = TRUE)
    ## cluster 1: a = mean(0.9, 0.8, 1.0); b = mean(0.2, 0.9), site 300 has coverage 2
    cl <- ev$evidence$clusters
    expect_equal(cl$sites, c(3L, 2L))
    expect_equal(cl$mean_level_a, c(mean(c(0.9, 0.8, 1.0)), mean(c(0.1, 0))))
    expect_equal(cl$mean_level_b, c(mean(c(0.2, 0.9)), mean(c(0.1, 0.6))))
    expect_match(ev$sections$regions, "chr1:100-300\t3\t0.90\t0.55", fixed = TRUE)
    ## differences: 4 sites covered in both; |b - a| >= 0.2 at 100 (0.7) and 2100 (0.6)
    expect_match(ev$sections$differences, "sites covered in both: 4", fixed = TRUE)
    expect_match(ev$sections$differences, "sites with |difference| >= 0.2: 2 of 4", fixed = TRUE)
    expect_match(ev$sections$differences, "chr1:100 0.90 -> 0.20", fixed = TRUE)
    ## levels of a: 0.9, 0.8, 1.0, 0.1, 0.0 -> < 0.2: 40%, > 0.8: 40% (0.9, 1.0), 0.2-0.8: 20%
    expect_match(ev$sections$levels, "a: median 0.80; < 0.2: 40%; 0.2-0.8: 20%; > 0.8: 40% (5 sites; 0 without",
                 fixed = TRUE)
    expect_match(ev$sections$levels, "b: median 0.40; < 0.2: 25%; 0.2-0.8: 50%; > 0.8: 25% (4 sites; 1 without",
                 fixed = TRUE)
    expect_match(ev$sections$coverage, "b: median 10, IQR 10-10; below 5: 20%", fixed = TRUE)
    expect_match(ev$sections$genomic_context, "not annotated", fixed = TRUE)
    ## the same object gives the same digest
    expect_identical(ev$text, epiSeeker:::.bm_digest(x, assay = "WGBS")$text)
})

test_that("the bmData digest finds the level assay, accepts an override and stops without one", {
    x <- bm_fixture()
    expect_equal(epiSeeker:::.bm_assays(x, NULL, NULL), list(value = "Methylation", coverage = "Cov"))
    expect_error(epiSeeker:::.bm_digest(x, value = "nope"), "not an assay")
    y <- suppressWarnings(bmData(value1 = cbind(a = c(5, 7)), gr = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(1, 2), width = 1)), sampleNames = "a", valueNames = "Cov"))
    expect_error(epiSeeker:::.bm_digest(y), "value =")
    ## one sample: no comparison
    data(demo_bmdata, package = "epiSeeker")
    ev <- epiSeeker:::.bm_digest(demo_bmdata)
    expect_match(ev$sections$differences, "one sample; no comparison", fixed = TRUE)
    expect_equal(ev$totals$n_sites, 50)
})

test_that("the bmData digest annotates the clusters when TxDb is given", {
    x <- bm_fixture()
    tx <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 5000), strand = "+", gene_id = "G1")
    ev <- suppressWarnings(suppressMessages(epiSeeker:::.bm_digest(x, TxDb = tx)))
    expect_match(ev$sections$regions, "\tfeature\tnearest_gene", fixed = TRUE)
    expect_equal(ev$evidence$clusters$gene, c("G1", "G1"))
    expect_match(ev$sections$genomic_context, "% of 2 clusters", fixed = TRUE)
})

test_that("interpret() reads a bmData through the shared skeleton", {
    skip_if_not_installed("aisdk")
    env <- new.env(); mock_schema_answer(env)
    res <- clusterProfiler::interpret(bm_fixture(), assay = "WGBS", model = "openai:m")
    expect_s3_class(res, "seq_interpretation")
    expect_equal(names(env$seen$schema$properties),
                 c("answer", "candidates", "methylation_profile", "sample_differences",
                   "genomic_context", "data_quality", "narrative"))
    expect_identical(env$seen$prompt, attr(res, "evidence")$text)
    expect_match(env$seen$system, "^You interpret methylation levels at single bases")
    expect_false("temperature" %in% names(env$seen))
    p <- clusterProfiler::provenance(res)
    expect_true(p$ok)
    expect_equal(p$entry_point, "interpret,bmData")
    expect_true(methods::existsMethod(clusterProfiler::interpret, "bmData"))
})

## ---- loopAnno ----------------------------------------------------------------------------------

test_that("the loopAnno digest reports gene hubs, classes and spans", {
    la <- lp_anno()
    ev <- epiSeeker:::.loop_digest(la)
    expect_equal(names(ev$sections), c("overview", "gene_hubs", "peak_links", "loop_classes",
                                       "anchor_distribution", "span"))
    expect_true(all(nchar(ev$sha256) == 64))
    expect_match(ev$sections$overview, "Loops: 5 between", fixed = TRUE)
    expect_match(ev$sections$overview, "Loops between chromosomes: 1", fixed = TRUE)
    ## G2 sits on L1 and L5 (other end distal both times); G1 and G3 on L3 (other end promoter);
    ## G4 (TSS 150,000) on anchor C of L2, whose other end is distal
    hub <- ev$evidence$gene_hubs
    expect_equal(hub$gene, c("G2", "G1", "G3", "G4"))
    expect_equal(hub$loops, c(2L, 1L, 1L, 1L))
    expect_equal(hub$to_promoter, c(0L, 1L, 1L, 0L))
    expect_match(ev$sections$gene_hubs, "(4 of 4, most loops first)", fixed = TRUE)
    expect_match(ev$sections$gene_hubs, "G2\t2\t0\t2", fixed = TRUE)
    expect_match(ev$sections$peak_links, "no peaks given", fixed = TRUE)
    expect_match(ev$sections$loop_classes, "promoter-promoter: 1 (20%)", fixed = TRUE)
})

test_that("the loopAnno digest reports peak links when they are given", {
    g <- lp_genes()
    la <- lp_anno()
    pk <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(42000, 10500, 60000, 152500, 117000),
                                                         c(42500, 10800, 60500, 152800, 117300)))
    pa <- suppressWarnings(suppressMessages(annotateSeq(pk, TxDb = g, verbose = FALSE)))
    links <- linkPeaksToGenes(pa, la)
    ev <- epiSeeker:::.loop_digest(la, links = links)
    nd <- sum(as.character(links$link) != "promoter")
    expect_match(ev$sections$peak_links, paste0("(", nd, " of 5 peaks are outside promoters)"), fixed = TRUE)
    expect_match(ev$sections$peak_links, "Genes reached only through loops (1 of 1", fixed = TRUE)
    expect_match(ev$sections$peak_links, "G2\t1\tG1", fixed = TRUE)
    expect_error(epiSeeker:::.loop_digest(la, links = data.frame(x = 1)), "linkPeaksToGenes")
})

test_that("interpret() reads a loopAnno through the shared skeleton", {
    skip_if_not_installed("aisdk")
    env <- new.env(); mock_schema_answer(env)
    res <- clusterProfiler::interpret(lp_anno(), model = "openai:m")
    expect_equal(names(env$seen$schema$properties),
                 c("answer", "candidates", "loop_profile", "gene_hubs", "peak_links", "narrative"))
    expect_identical(env$seen$prompt, attr(res, "evidence")$text)
    expect_match(env$seen$system, "^You interpret chromatin loops")
    expect_equal(clusterProfiler::provenance(res)$entry_point, "interpret,loopAnno")
    ## an unfinished answer is a failure, not a result
    local_mocked_bindings(generate_object = function(...) list(object = list(answer = "x"),
                                                               finish_reason = "length", raw_text = "{"),
                          .package = "aisdk")
    expect_warning(res2 <- clusterProfiler::interpret(lp_anno(), model = "openai:m"), "no complete answer")
    expect_false(clusterProfiler::provenance(res2)$ok)
})

test_that("the region methods still send the system prompt they always sent", {
    expect_match(epiSeeker:::.seq_system_prompt(),
                 "^You interpret genomic regions \\(peaks from ChIP-seq, ATAC-seq, CUT&Tag or similar\\)\nfrom")
    expect_identical(epiSeeker:::.seq_system_prompt(), epiSeeker:::.seq_system_prompt(NULL))
})

## ---- linkPeaksToGenes strict options ----------------------------------------------------------

## independent re-implementation of checks/20261107_fig4-loops/fig4.R (loop_tables, loop_genes,
## nearest_pc and the classes of distal peaks) on plain GRanges
fig4_classes <- function(gi, peaks, tss, flank = 1000, distal = 3000) {
    a1 <- InteractionSet::anchors(gi, type = "first"); a2 <- InteractionSet::anchors(gi, type = "second")
    reg <- unique(c(GenomicRanges::granges(a1), GenomicRanges::granges(a2)))
    i1 <- match(a1, reg); i2 <- match(a2, reg)
    h <- GenomicRanges::findOverlaps(reg + flank, tss, ignore.strand = TRUE)
    rg <- split(tss$gene_id[S4Vectors::subjectHits(h)], factor(S4Vectors::queryHits(h), levels = seq_along(reg)))
    partners <- split(c(i2, i1), factor(c(i1, i2), levels = seq_along(reg)))
    ctr <- GenomicRanges::resize(peaks, 1, fix = "center")
    d <- GenomicRanges::distanceToNearest(ctr, tss, ignore.strand = TRUE)
    near <- tss$gene_id[S4Vectors::subjectHits(d)]
    dist <- S4Vectors::mcols(d)$distance
    vapply(seq_along(peaks), function(i) {
        if (dist[i] <= distal) return("promoter")
        own <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(peaks[i], reg, ignore.strand = TRUE))
        if (!length(own)) return("no loop")
        pg <- setdiff(unique(unlist(rg[unique(unlist(partners[own]))])), unique(unlist(rg[own])))
        if (!length(pg)) return("loop, no gene")
        if (near[i] %in% pg) "loop includes nearest" else "loop, other genes"
    }, character(1))
}

test_that("the strict options of linkPeaksToGenes reproduce the Fig 4 definitions", {
    g <- lp_genes()
    la <- lp_anno()
    gi <- la@loops
    pk <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(42000, 10500, 60000, 152500, 117000, 205000, 9000),
                                                         c(42500, 10800, 60500, 152800, 117300, 205300, 9200)))
    pa <- suppressWarnings(suppressMessages(annotateSeq(pk, TxDb = g, verbose = FALSE)))
    pc <- c("G1", "G2", "G3")                         # G4 is "non-coding"
    strict <- linkPeaksToGenes(pa, la, TxDb = g, gene_ids = pc, anchor_flank = 1000,
                               distal_distance = 3000, exclude_own_anchor = TRUE)
    ## only cis loops between non-overlapping anchors enter fig4; the trans loop L4 is left out there
    cis <- !is.na(S4Vectors::mcols(gi)$span)
    tss <- epiSeeker:::.tss_of(g, pc)
    ref <- fig4_classes(gi[cis], pa@anno, tss)
    la_cis <- la; la_cis@loops <- gi[cis]
    strict_cis <- linkPeaksToGenes(pa, la_cis, TxDb = g, gene_ids = pc, anchor_flank = 1000,
                                   distal_distance = 3000, exclude_own_anchor = TRUE)
    expect_equal(as.character(strict_cis$link), ref)
    ## G4 (non-coding) never appears; peak 4 is near G4 only and is distal under the strict rules
    expect_false(any(grepl("G4", strict$loop_genes)))
    expect_false(any(strict$nearest_gene %in% "G4"))
    expect_equal(as.character(strict$link[4]), ref[4])
    ## defaults are unchanged
    expect_identical(linkPeaksToGenes(pa, la), linkPeaksToGenes(pa, la, exclude_own_anchor = FALSE))
    expect_error(linkPeaksToGenes(pa, la, anchor_flank = 1000), "need `TxDb`")
    expect_null(strict$loop_symbols)                    # no annoDb: ids only
})
