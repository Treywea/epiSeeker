---
name: range-atac
description: Analyse chromatin accessibility peak sets (ATAC-seq, DNase-seq) and, if given, their fragments with the epiSeeker R package - accessibility at promoters, distance to TSS, TSS profile, comparison of sets, overlap with reference sets in a background, and from fragments the TSS enrichment, fraction of fragments in peaks and fragment lengths - with checks for the silent errors (wrong genome build, seqlevels of peaks or fragments that do not match, no background), then interpret() and provenance. Use when a user has ATAC/DNase peaks or a fragments file and asks where chromatin is open, near which genes, how conditions differ, or whether the library looks right. Single-cell clusters go to range-singlecell, spatial domains to range-spatial.
---

# range-atac: accessibility with epiSeeker

Accessible regions say where chromatin is open, not which factor binds there or whether an element
is active. The script `scripts/run_range_atac.R` runs the steps below; prefer it. The steps it
shares with `range-tf` and `range-histone` are in `../shared/range-steps.R`.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`) and equal to `GenomeInfoDb::genome(TxDb)`; TxDb named.
2. **Seqlevels** of the peaks and of the fragments share a style with the TxDb (fragments files
   often carry contigs the TxDb lacks; those are reported, not an error).
3. **Peak width.** ATAC / DNase peak calls are narrow (the LNCaP DNase example: median 180 bp); a
   median above 1 kb suggests a broad-mode caller or merged regions (NOTE).

## 2. Analysis (epiSeeker functions only)
- `readPeakFile()`; `annotateSeq(peaks, TxDb =, annoDb =, tssRegion = c(-3000, 3000))`;
  several sets with `csAnnoList()`; `plotAnnoBar()`, `plotDistToTSS()`.
- TSS profile of the peaks: `getTagMatrix(peaks, upstream = 3000, downstream = 3000,
  type = "start_site", by = "gene", TxDb =)`, `plotPeakProf()`.
- With `--fragments`: `readFragments()`; Tn5 insertion sites = both fragment ends; their profile
  with `getTagMatrix(insertions, upstream = 2000, downstream = 2000, ...)`. TSS enrichment = mean
  insertions per bp within +/-50 bp of the TSS over the outer 100 bp at each end. FRiP = share of
  fragments overlapping a peak. Fragment length median and share < 147 bp (nucleosome-free).
- Overlap with reference sets (e.g. TF ChIP peaks): as in `range-tf`, background required.

## 3. Correctness checks
- `annoChecks()` of every set: build suspect -> **stop**; dropped peaks -> report.
- TSS enrichment next to ENCODE's ATAC standard for GRCh38 (< 5 concerning, 5-7 acceptable, > 7
  ideal) and FRiP next to ENCODE's 0.3, as NOTEs: the TSS set and normalisation differ from
  ENCODE's, so they are references, not thresholds.
- Every warning of an epiSeeker call as a WARN row.

## 4. interpret() last
`clusterProfiler::interpret(anno, context =, assay = "ATAC", question =, overlap =, model =)`;
skipped without a model. A call counts only when `clusterProfiler::provenance(res)$ok` is TRUE.

## 5. Provenance and outputs
`figures/` (annotation_bar, distance_to_tss, tss_profile, tss_insertions), `tables/`
(annotation_<set>, annotation_stat, fragments_qc, overlap), `report/report.md`, `provenance/`
(run.txt, checks.tsv, sessionInfo.txt, interpret_provenance.txt, objects.rds; the simulated
fragments of `--example_fragments`).

## Do not say
- "factor X binds here" or "binding sites" from accessibility alone; that needs ChIP or a footprint.
- "active enhancer", "activates", "regulates" from openness alone.
- "the target gene" for the nearest gene.
- A QC verdict from TSS enrichment or FRiP of the simulated fragments of the example.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/range-atac/scripts/run_range_atac.R", package = "epiSeeker"))')
Rscript "$S" --example --out atac_example                       # ENCODE LNCaP DNase, chr1-2
Rscript "$S" --example_fragments --out atac_fragments_example   # + simulated fragments
Rscript "$S" --peaks atac_peaks.narrowPeak --fragments fragments.tsv.gz --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db --assay ATAC --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
