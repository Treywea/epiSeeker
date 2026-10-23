# epiSeeker 1.1.4

+ `annotateSeq(annoDb = "org.Hs.eg.db")` adds the gene annotation (SYMBOL,
  GENENAME, ENSEMBL) again under R >= 4.6. The package was loaded with
  `requireNamespace(annoDb, character.only = TRUE)`; `character.only` is an
  argument of `require()`, which R 4.6's `loadNamespace()` rejects, so
  `requireNamespace()` returned `FALSE`, the database was never looked up, and
  every call ended with the misleading warning "ID type not matched". The
  database is now taken from the package namespace, so it need not be attached
  either, and an `OrgDb` object is accepted as well as its name. The same fix
  applies to the BSgenome lookup of `getBmMatrix()`. (2026-10-23, Fri)
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
