# epiSeeker 1.1.4

+ Fixes from rerunning the 2025 case studies (2026-11-09):
  - `enrichPeakOverlap()` ignores strand in every method and column. The observed overlap was
    counted with strand (`intersect()` pairs `*` only with `*`) while the null shuffled unstranded
    sets, so a target set with `+`/`-` strand (genes, transposable elements) gave `N_OL = 0` and
    `p = 1` without a warning; `method = "fisher"` reported `N_OL = 0` the same way.
    `attr(, "params")$strand` is `"ignored"`.
  - `getMotifMatrix()` with more than one region placed the hits at wrong coordinates (the
    region starts were recycled over the rows of the hits). Each hit is now placed in its own
    region (new column `region`); hits of `-` strand regions are mapped back to the genome;
    `attr(, "range")` spans the regions, `attr(, "regions")` keeps them. `plotMotifProf()` stops
    for hits on several chromosomes.
  - The genome-build check of `annotateSeq()` also runs when the annotation is a `GRanges`
    (custom TE or enhancer sets), against its `seqlengths()`, and `getTagMatrix()` runs it for
    every set of regions, against the chromosome lengths of the annotation the regions were made
    from (kept by `getBioRegion()`). Chromosome lengths recorded on the peaks that differ from the
    annotation's count as evidence too. A `GRanges` without lengths cannot be checked: a message
    says so and `annoChecks()` records it (`seqlengths_known = FALSE`). `annoChecks()` accepts
    tag matrices. hg38 peaks against TEs in T2T-CHM13 coordinates (the 2025 case) now warn.
  - `getTagMatrix(type = "body")` warns when regions with peaks are left out because they are
    not longer than the number of bins (most Alu elements with `nbin = 800`), with the number and
    fraction, and when rows have an empty body (flanks in bp); recorded in
    `annoChecks(<tag matrix>)$short_regions`. Which regions are left out is unchanged. Flanks
    given in bp with a fractional number of bins (`nbin * 0.15`, e.g. `nbin = 50`) were columns
    of zeros; they get `round(nbin * 0.15)` bins.
  - `getTagMatrix(type = "body")` no longer expands every region to one value per base: the bins
    are summed on the run-length coverage, and regions without peaks are skipped. Same matrix
    (identical on the human genes with 5,000-25,000 peaks); on those inputs 1.3-1.6 GB and
    1.6-3.9 s instead of 15.7-16.3 GB and 44-91 s.

+ Agent skills in `inst/skills/` (`system.file("skills", package = "epiSeeker")`): `range-tf`,
  `base-methylation`, `interaction-loops` and `range-singlecell`, each a `SKILL.md` (Agent Skills
  format) with one runnable script and a worked example on the package's data. All follow one
  structure: input checks (genome build, seqlevels, assay), epiSeeker functions only, the built-in
  correctness checks (`annoChecks()` etc.), `interpret()` last, provenance. The scripts stop on a
  failed check, record every warning, and call a model only when one is given. A test checks the
  format and that every function the skills name exists (2026-11-09).

+ `enrichPeakOverlap(method = "fisher")` tests every target set against the `background` as a
  universe, without permutations: the units are the background regions (overlapping ones merged,
  adjacent ones such as tiles kept apart) plus every query region outside them (so the table stays
  consistent, where LOLA's raw mode fails with "Negative b entry"), and the one-sided
  Cochran-Mantel-Haenszel test over `strata` bins of unit width (default 50) gives the p-value.
  Stratifying by width keeps the test calibrated on universes of unequal region sizes: with peaks
  placed at random within a background of promoters and enhancers, 5.4% of the tests had p < 0.05
  (24.7% without strata, on development data); it has less power than the permutation test for
  weak overlaps. With a background the result also has `odds_ratio`, the Mantel-Haenszel odds
  ratio of that table; `attr(, "params")$rank_by` names the statistic to rank by, which the
  `interpret()` digest uses and states. It stays `jaccard`: `odds_ratio` put the profiled factor
  first for 13 of 14 ENCODE peak sets against a universe of the 60 reference sets (`jaccard` 11),
  but in simulation it put the true set first less often at 2% shared peaks (0.78 against 1.00;
  related and small random sets came first) and more often at 1% (0.40 against 0.02).
  The permutation test with `background` is about four times faster (positions are drawn from
  precomputed cumulative widths and all target sets are shuffled and counted in one pass per
  permutation; same statistic, strand ignored), and `N_OL` and the effect sizes are computed from
  overlap pairs. The default (`method = "permutation"`, `background = NULL`) gives the same
  numbers as before. (2026-11-09, Mon)
+ `annotateLoops()` takes `gene_ids` to restrict the genes recorded at the anchors (for example
  to protein-coding genes); an anchor is then a promoter anchor only when it holds the promoter
  of one of them, so pseudogenes and uncharacterised loci no longer make an anchor a promoter
  anchor or appear among the gene hubs of `interpret(<loopAnno>)`, whose digest states the
  restriction. The default is unchanged. (2026-11-08, Sun)
+ Single-molecule modification data can be read: `readModCalls()` takes the per-read tables of
  `ft extract --all` (fibertools, Fiber-seq: m6A, 5mC, nucleosomes, MSPs and FIREs) and of
  `modkit extract full` / `extract calls` (any modBAM) and returns one range per call, with the
  read, its haplotype and the source in the columns and `metadata()`; coordinates are made
  1-based and calls without a reference position are counted, not kept. `summarizeModCalls()`
  turns them into a `bmData` (coverage and fraction modified), and `readBedMethyl()` reads a
  `modkit pileup` file, so `getBmMatrix()` and `plotBmProf()` work on both unchanged.
  `modReadMatrix()` gives one value per read and region, `modCooccurrence()` tests whether two
  regions open on the same molecules (Fisher's test on the reads spanning both) and
  `plotModReads()` draws the reads over a region. BAM/CRAM modification tags are not decoded
  in the package; the tables come from modkit or fibertools. (2026-11-07, Sat)
+ Single-cell and spatial chromatin data enter as groups of a `csAnnoList`. `readFragments()`
  reads a fragments file; `fragmentsToGroups()` sums the fragments of the cells or spots of each
  cluster or spatial domain in tiles and keeps the group-specific tiles, matching barcodes with
  or without the `-1` suffix of fragments files; `peaksToGroups()` takes per-group peaks or the
  marker tables of ArchR, Signac and 10x; `annotateGroups()` annotates every group and records
  the resolution, the number of cells or spots and spot coordinates in `metadata()`, read with
  `groupInfo()`. `readSpotPositions()` and `plotSpatialGroups()` read and plot spot positions.
  The clustering itself stays with the upstream tools. `interpret()` states the resolution and
  the cells or spots behind every group. No new class and no new dependency. (2026-11-07, Sat)
+ `interpret()` reads a `bmData` and a `loopAnno` (experimental: neither has been compared
  with a baseline that uses no model), so every group -- interval, base and pair -- has a
  method on the shared skeleton. The `bmData` digest merges sites into clusters and reports
  their mean level per sample (with `TxDb`, their feature and nearest gene), the differences
  between samples, the distribution of levels and the coverage; the level and coverage assays
  are found by their values or named with `value =` / `coverage =`. The `loopAnno` digest
  reports the genes at promoter anchors ranked by loops, the peak links of
  `linkPeaksToGenes()` when given (`links =`), the loop classes, the anchor distribution and
  the span per class. `linkPeaksToGenes()` gains `TxDb`, `gene_ids`, `anchor_flank`,
  `distal_distance` and `exclude_own_anchor`, which reproduce the definitions used to score
  loops against a CRISPRi enhancer-gene benchmark (protein-coding genes; genes of an anchor =
  transcript starts within the anchor +/- 1 kb; nearest gene from transcript starts; distal =
  more than 3 kb from one); the defaults are unchanged. (2026-11-07, Sat)
+ `enrichPeakOverlap()` without `background` says, once per session, that the p-values of a
  whole-genome null are anti-conservative for peaks concentrated in open chromatin (on 20
  ENCODE TF ChIP-seq peak sets against 60 reference sets, 1181 of the 1200 pairs were
  significant) and suggests `background =` and ranking by `jaccard`, which put the profiled
  factor first for 11 of 14 sets against 3-4 for `fold_enrichment`. The default is
  unchanged; `verbose = FALSE` or `options(epiSeeker.null_note = FALSE)` turns the note off.
  (2026-11-07, Sat)
+ The genome-build check of `annotateSeq()` records `n_peaks` and `reliable` in
  `annoChecks()`: below 200 peaks it is unreliable (only peaks near chromosome ends can fall
  outside them; the mismatched sets it missed had 50-160 peaks), and `show()` says so. No
  warning is raised for that alone, so small example and test sets stay quiet.
  `interpret()` now puts these checks in the overview of its evidence digest (a suspected
  build mismatch, dropped peaks, an unreliable check) and tells the model to state them
  first; for an overlap table it states the null the p-values came from, warns that a
  whole-genome null makes most reference sets significant, and ranks by `jaccard`, giving
  `fold_enrichment` as an effect size. (2026-11-07, Sat)
+ `annotateSeq()` checks two silent errors and keeps what it finds, read with
  `annoChecks()`. It warns when the peaks look like coordinates of another genome build than
  `TxDb`, that is when peaks extend beyond the chromosome ends of `TxDb`. An annotation
  against the wrong build looks normal (hg19 example peaks annotated with the hg38 `TxDb` of
  the vignette had 15% of their peaks within 1 kb of a TSS instead of 58%). On 116 ENCODE
  and example files annotated with their own build and another one, the check raised no
  false alarm on the 110 matched pairs and found 132 of the 136 mismatched ones; the four it
  missed were sets of 50-150 peaks. The warning on dropped peaks now says how many were
  dropped on which sequences and why (not in `TxDb`, named differently, no gene on the
  sequence). The checks take under 1% of the time of `annotateSeq()` and never change the
  annotation. (2026-11-05, Thu)
+ `enrichPeakOverlap()` gains `background`: the regions target peaks are shuffled into under
  the null (open chromatin, mappable regions, exons). Shuffling over whole chromosomes ignores
  that peaks of real assays sit in a small part of the genome, so almost every target set came
  out significant (1,193 of 1,200 query-target pairs on ENCODE data). With a background each
  target set gets its own null, and the result reports `expected` (mean overlap under the null)
  and `fold_enrichment`, which, unlike the p-value, can rank target sets. Fixed with
  `pool = FALSE`: the adjusted p-values are corrected across all target sets (each call saw one
  p-value), every target set keeps its name (all rows were `targetPeak1`), and lifted targets are
  no longer lifted a second time. Fixed in `shuffle()`: chromosome lengths are taken by name;
  when the seqlevels were not in alphabetical order (as `readPeakFile()` returns them for a file
  that is not sorted by chromosome, or chr1, chr2, ..., chr10) peaks were placed along the length
  of another chromosome; for the CBX7 example peaks with seqlevels chr1, chr2, ..., 272 of
  1,649 shuffled peaks fell past the end of their chromosome. This changes the p-values for such
  input. Defaults are unchanged. With `pool = FALSE` and target files, `tSample` is now the file name,
  as with `pool = TRUE`. Two p-value errors are fixed for every setting: a null value equal to the
  observed overlap now counts against it (`>=`, not `>`; for small target sets, whose overlap
  is a small integer, 20% of unrelated sets were at p < 0.05, and 3.5% are now), and the
  mirrored direction of `symmetric = TRUE` shuffles the whole query set in each permutation
  and compares overlap counts, where it shuffled a single query peak and compared its 0/1
  overlap with the observed ratio (two sets without any overlap came out at p = 0.04). Both
  change existing p-values. (2026-11-05, Thu)
+ Chromatin loops, the pair group next to regions (`csAnno`) and single bases (`bmData`).
  `readLoops()` reads BEDPE files (also gzipped, with a Juicer/HiCCUPS `#` header, a plain
  header or none), data frames, `GInteractions` and `Pairs`; loops are InteractionSet's
  `GInteractions`, the class HiCExperiment and GenomicInteractions use, so InteractionSet is
  now imported. `annotateLoops()` annotates both anchors with `annotateSeq()`, marks an
  anchor as a promoter anchor when it overlaps the promoter window (`tssRegion`) of a
  transcript of `TxDb`, records the genes (and, with `annoDb`, symbols) of those promoters,
  and classes each loop as promoter-promoter, promoter-distal or distal-distal; the result
  is a `loopAnno`, whose anchors are an ordinary `csAnno`. `linkPeaksToGenes()` gives each
  peak of a `csAnno` the genes at the other end of the loops it sits in and compares them
  with its nearest gene; `plotLoopAnno()` and `plotLoopLinks()` plot the loop classes, the
  loop spans and the links. On K562 (ENCODE HiCCUPS loops and ATAC-seq peaks), 22% of the
  87,875 peaks outside promoters sit in a loop anchor, and for 10,842 of them the loops
  point to genes other than the nearest one. (2026-11-05, Thu)
+ `interpret()` reads a `csAnno`: a method on the `interpret()` generic of clusterProfiler, which
  epiSeeker now imports. The annotation (and, optionally, an `enrichPeakOverlap()` result) is
  turned into an evidence digest without a model -- top regions with their nearest genes,
  promoter genes, overlap with reference samples ranked by effect size rather than by
  p-value, genomic distribution, distance to TSS and region width -- in which every section
  states how many items it sends and how many there are. One structured call answers an
  optional `question` first (`answer`, ranked `candidates`), then describes the regions.
  A call counts as successful only when the model finished and returned an object; the
  execution record is read with `clusterProfiler::provenance()`. Needs the S4 generic of
  clusterProfiler (Treywea/clusterProfiler feat/interpret-s4-generic) and the `aisdk`
  package. A `GRanges` is annotated with `annotateSeq()` first (`TxDb` must be given, so
  regions are never annotated against a default genome); a `csAnnoList` is compared group
  by group (experimental: distributions, widths, distances to TSS, shared and group-specific
  promoter genes, top regions); a `peakOverlapResult` alone is refused with a pointer to
  `interpret(<csAnno>, overlap = )`. The digest states the genome and annotation source and
  how the overlap p-values were obtained. (2026-11-03, Tue)
+ Result objects record more of how they were made, following the result classes of
  clusterProfiler. `annotateSeq()` stores the genome, the source, organism and creation time
  of `TxDb` and the name and version of `annoDb` in `metadata()` of the annotated regions,
  read with `annoSource()` and shown by `show()`; objects made earlier fall back to the
  genome stored with the regions. The `csAnno` class is unchanged: a new slot would leave
  every `csAnno` saved by an earlier version without it, and `saveRDS()` then fails on such
  an object ("no slot of name ..."). `enrichPeakOverlap()` adds three effect sizes that need
  no permutation -- `jaccard`, `frac_query` and `frac_target` -- because its p-value cannot
  rank target sets (its smallest value, 1/(nShuffle+1), is shared by most of them); the
  result is a `peakOverlapResult`, still a data.frame, with the arguments of the test in
  `attr(, "params")`. With the same seed the existing columns are identical to before, for
  `pool = TRUE` and `pool = FALSE`. A new class `csAnnoList` (made with `csAnnoList()`) holds
  named `csAnno` objects for comparison; it is a list, so `plotAnnoBar()`,
  `plotDistToTSS()` and `lapply()` accept it as before. (2026-11-03, Tue)
+ The example peak files are hg38, like the `TxDb` the vignette and the tests annotate them
  with. The five GEO files of `getSampleFiles()` (GSE40740 CBX6/CBX7, GSE48308 ARmo) and
  `sample_peaks.txt` were hg19 peaks annotated against hg38 genes: CBX6_BF had 15% of its
  peaks within 1 kb of a TSS instead of 58%, and its Reactome enrichment (kidney
  development, gastrulation) was gone. They are lifted with the UCSC hg19ToHg38 chain,
  keeping the peaks that map to one range of the same width (98.5-100% per file,
  `data-raw/lift_sample_data_hg38.R`); `tagMatrix` and `peakAnnoList` are rebuilt with
  `TxDb.Hsapiens.UCSC.hg38.knownGene`, and the outputs pasted in the vignette are updated.
  Annotated as published (hg19) and as shipped (hg38), the same peaks fall in the same
  feature category for 96-98% of them, and CBX6 gives the same top Reactome pathways.
  (2026-11-03, Tue)
+ `annotateSeq(annoDb = "org.Hs.eg.db")` adds the gene annotation (SYMBOL,
  GENENAME, ENSEMBL) again under R >= 4.6. The package was loaded with
  `requireNamespace(annoDb, character.only = TRUE)`; `character.only` is an
  argument of `require()`, which R 4.6's `loadNamespace()` rejects, so
  `requireNamespace()` returned `FALSE`, the database was never looked up, and
  every call ended with the misleading warning "ID type not matched". The
  database is now taken from the package namespace, so it need not be attached
  either, and an `OrgDb` object is accepted as well as its name. The same fix
  applies to the BSgenome lookup of `getBmMatrix()`, which also no longer drops
  genomes whose organism starts with B or S (`Scerevisiae`): its pattern
  `"^[^BSgenome]"` was a character class, not "does not start with BSgenome".
  (2026-10-23, Fri)
+ `annotateSeq()` keeps `annotation`, `geneChr/geneStart/geneEnd`, `geneId` and
  `transcriptId` consistent with each other at transcript level (issue #252 of
  ChIPseeker). The transcript id of an exon/intron/UTR hit was taken from
  `names(genomicRegion)[subjectIndex]`, indexing the unlisted ranges with the
  names of the `GRangesList`; this mostly returned NA and occasionally a wrong
  transcript, so a peak could be reported with the metadata of an unrelated
  transcript (even on another chromosome). The ids are now expanded before
  indexing, and the same alignment is applied at `level = "gene"`, where the gene
  of the overlapping exon/intron/UTR is resolved with `TXID2EG()` and used
  instead of the nearest gene. (2026-10-01, Thu)
+ `seq2gene()` no longer fails with `$ operator is invalid for atomic vectors`
  when none of the queried regions overlaps an exon
  (`getGenomicAnnotation.internal()` returns `NA` then); host-gene extraction is
  skipped and the nearest/flanking genes are still reported (issue #248 of
  ChIPseeker). (2026-10-01, Thu)
+ `downloadGEObedFiles()`/`downloadGSMbedFiles()` rewrite the `ftp://` urls of
  `gsminfo$supplementary_file` to `https://` before downloading and report the
  underlying download error when a file cannot be fetched (issue #254 of
  ChIPseeker). (2026-10-01, Thu)
+ Documented the distance definition of `flank_gene_distances` reported by
  `annotateSeq(..., addFlankGeneInfo = TRUE)`: a distance of 0 means that the
  peak overlaps the feature range - at `level = "transcript"` the feature is the
  whole transcript, which is why most entries can be 0 - while non-overlapping
  peaks get the signed distance to the feature TSS (issue #235 of ChIPseeker).
  (2026-10-01, Thu)
+ `annotateSeq()` no longer fails with `Error: invalid subscript` when *every*
  peak is dropped, which typically happens because none of the seqlevels of the
  peaks matches `TxDb` (`chr1` vs `NC_000001.11`): `.get_distance_to_gene_end()`
  no longer hands the `SortedByQueryHits` object returned by `follow()` for an
  empty query to `features[]`, and `annotateSeq()` reports the seqlevels mismatch
  together with how to align the chromosome names (`seqlevelsStyle()`).
  (2026-10-01, Thu)
+ `enrichPeakOverlap()`/`enrichAnnoOverlap()` accept a single `GRanges` as
  `targetPeak`: it was passed on unwrapped while the overlap code works on a
  list of target peak sets, so the call failed with "GRanges objects don't
  support [[, as.list(), lapply()". The documentation of `enrichPeakOverlap()`
  and `enrichAnnoOverlap()` now also states the direction of the test (the
  observed ratio is the fraction of *target* peaks covered by the query peaks and
  the target is the shuffled set, so `N_OL` is direction free while the p-value
  is not) and that `N_OL` of `enrichAnnoOverlap()` counts genes and can therefore
  exceed the number of input peaks (issue #84 of ChIPseeker).
  (2026-10-01, Thu)
+ An empty input is now handled: `annotateSeq(GRanges(), ...)` returns an empty
  `csAnno` instead of failing with "Error: invalid subscript". The root cause
  (follow()/precede() returning a SortedByQueryHits object instead of indices for
  an empty query) is guarded in `getNearestFeatureIndicesAndDistances()`, which
  also makes `seq2gene()` return `character(0)` for an empty input.
  (2026-10-01, Thu)
+ `enrichPeakOverlap()` gained an opt-in `symmetric` argument. With
  `symmetric = TRUE` the mirrored direction (query and target exchanged) is
  computed as well and the two one-sided permutation p-values are combined as
  `min(1, 2*min(p, p_rev))` (Hedges), so the result no longer depends on the order
  of the arguments. The default stays `FALSE` (one-sided) and all existing numbers
  are unchanged; the mirrored test doubles the number of permutations and raises
  the smallest reportable p-value to `2/(nShuffle+1)`. `enrichAnnoOverlap()` needs
  no such argument, its hypergeometric p-value is already exactly symmetric (issue
  #84 of ChIPseeker). (2026-10-01, Thu)
+ `upsetplot(vennpie = TRUE)` can be drawn again with ggplot2 >= 4.0. The vennpie
  sub-view is embedded as an `annotation_custom()` layer, which ggplot2 only
  supports below `coord_cartesian()`. The sub-view no longer uses `coord_fixed()`
  and the parent plot is converted with `ggplotify` as well, because the ggupset
  coordinate system (`CoordCombMatrix`) rejects such layers. No replacement for
  the removed `coord_fixed()` is added: ggplotify rasterises the grob with the
  aspect ratio of the device, and forcing a square panel distorted the venn
  diagram further (measured anisotropy 1.43 without versus 1.50 with
  `theme(aspect.ratio = 1)`, against 1.41 for the undistorted drawing). Before,
  the plot failed with `` `annotation_custom()` only works with
  `coord_cartesian()` ``. Note that ggupset 0.4.1 itself still maps the
  deprecated `size` aesthetic, so drawing the plot may print
  "Please use `linewidth` instead"; the warning comes from ggupset
  (https://github.com/const-ae/ggupset/issues) and can be ignored.
  (2026-10-01, Thu)

# epiSeeker 1.1.3

+ `annotateSeq()` now reports `geneChr` and `geneStrand` as characters instead of
  factor codes. `as.data.frame()` returns 'seqnames'/'strand' as factors and
  assigning a factor into `mcols()` dropped the class, so both columns came out
  as integers (e.g. `geneStrand` = 1/2 instead of +/-, and a wrong `geneChr`
  whenever the seqlevels were not in numeric order). (2026-09-15, Tue)
+ `annotateSeq()` now warns when peaks are dropped because no feature of `TxDb`
  can be found for them (e.g. peaks on contigs/scaffolds without genes, or a
  seqlevels style mismatch). Previously they disappeared silently. (2026-09-15, Tue)
+ `annotateSeq(..., sameStrand = TRUE)` no longer assigns a peak to a feature on
  the opposite strand. Overlap detection in `getNearestFeatureIndicesAndDistances()`
  was calling `findOverlaps()` with `unstrand(features)`, which bypassed `sameStrand`
  and overrode the strand-aware nearest-feature result. Peaks with ambiguous
  strand (`*`) are unaffected and still match features on any strand.
  (2026-09-15, Tue)
+ `plotAnnoBar()` no longer uses the deprecated `ggplot2::aes_string()`. It follows
  the tidy evaluation idiom already used by `plotDistToTSS()`. (2026-09-15, Tue)
+ dropped unused `@importFrom` directives with no call site (`aplot::xlim2`,
  `ggplot2::geom_segment`, `ggplot2::geom_text`, `ggplot2::scale_fill_hue`,
  `utils::getFromNamespace`) and a duplicate `ggplot2::scale_fill_brewer`. The
  `grid::unit()` import, previously declared in `upsetplot()` which never used it,
  now sits in `plotBmProf()` where `unit()` is actually called, keeping `grid` a
  genuine import. (2026-09-17, Thu)

# epiSeeker 1.1.2

+ documentation build migrated to roxygen2 8.1.0 (re-generated Rd pages and
  NAMESPACE; fixed `csAnno` class `@aliases` to be a single line); the exported
  API is unchanged
+ fixed bug in `getNearestFeatureIndicesAndDistances()` where results for
  `overlap=="all"` were silently overridden by the `overlap=="TSS"` branch,
  causing the two modes to behave identically (2026-09-07, Mon)

# epiSeeker 1.1.1

+ make csAnno subset robust to GRanges metadata columns (2027-07-07, Tue)

# epiSeeker 1.0.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# epiSeeker 0.99.12

+ add more runnable examples and increase the coverage of tests (44.18%) (2025-11-30, Sun)
+ use new demo data for base modification (`demo_bmdata`) and move some dependency packages from 'Imports' to 'Suggests' (2025-11-10, Mon)
+ consistently use 'hg38' for demo with new demo file, a small subset derived from GSM6418464 (2025-11-07, Fri)
+ fixed Notes reported by `BiocCheck` and add 'GeneRegulation' to 'biocViews' (2025-11-07, Fri)
+ update test for the changes of TxDb.Hsapiens.UCSC.hg19.knownGene (2025-11-06, Thu)
+ new cache mechanism from 'yulab.utils' (2025-10-15, Wed)
+ fixed R check (2025-10-05, Sun)
+ remove dependent of `Vennerable` package (2025-09-27, Sat)
+ epiSeeker inherits from ChIPseeker to support analysis of multi-omics epigenomic data (2025-09-19, Fri)
    - epiSeeker inherits functions from ChIPseeker to plot peak profile, plot coverage, peak annotaions, visualization of annotation, data mining and statistic for peak overla
    - epiSeeker merges multiple functions to plot profile (getTagMatrix + plotPeakProf/plotPeakHeatmap), and add featrues of dealing missing data and different statistic methods to combind data
    - epiSeeker provides new features for plotCov to plot cluster tree and peak co-accessibility
    - epiSeeekr provides interactive function for user to explore data
    - epiSeeker provides functions to plot gene structure, base modification and motif
    - see more on <https://github.com/YuLab-SMU/ChIPseeker/blob/devel/NEWS.md>
