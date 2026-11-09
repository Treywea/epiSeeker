---
name: range-tf
description: Annotate, profile, compare and test the overlap of narrow peak sets (TF ChIP-seq, CUT&RUN, CUT&Tag for TFs, ATAC summits) with the epiSeeker R package, catching the silent errors that make such analyses look fine but be wrong (wrong genome build, seqlevels that do not match, a whole-genome null for overlap tests), then interpret the result with interpret() and keep provenance. Use when a user asks what a peak set binds near, which genes or programs it points to, how several peak sets differ, which factor or reference sets it resembles, or to "run epiSeeker/ChIPseeker" on peaks.
---

# range-tf: narrow peaks with epiSeeker

Every number comes from an epiSeeker function; your job is to call them in the right order, stop
when a check fires, and report what the functions returned. The script
`scripts/run_range_tf.R` does all of the steps below; prefer running it over re-typing them.
Broad histone marks (H3K27me3, H3K36me3, H3K9me3) are not covered: their peaks are kilobases wide and
"nearest gene" means little for them.

## 1. Input checks (stop on failure)
1. **Genome build.** Ask for it (or read it from the pipeline that made the peaks) and pass it as
   `--genome`. It must equal `GenomeInfoDb::genome(TxDb)`. Bed files almost never record their
   build, so nothing else will notice hg19 peaks annotated with an hg38 TxDb.
2. **TxDb named explicitly.** `annotateSeq()` with `TxDb = NULL` silently falls back to
   `TxDb.Hsapiens.UCSC.hg38.knownGene` (mouse data would be annotated against human). Always pass
   `TxDb =` and `annoDb =`.
3. **Seqlevels.** `GenomeInfoDb::seqlevelsStyle(peaks)` must share a style with the TxDb ("1" vs
   "chr1" drops every peak). Convert only on purpose (`--fix_seqlevels`, recorded).
4. **Assay.** TF peaks are narrow: the median peak width of 20 ENCODE TF ChIP-seq sets was
   184-556 bp. A median far above that suggests a broad mark; say so instead of going on.
   Promoter fraction does not tell TFs apart (25-92% in the same 20 sets).

## 2. Analysis (epiSeeker functions only)
- `readPeakFile()` each file; `annotateSeq(peak, TxDb =, annoDb =, tssRegion = c(-3000, 3000))`.
- Several sets: `csAnnoList(...)`, then `plotAnnoBar()`, `plotDistToTSS()`.
- Profile around TSSs: `getTagMatrix(list_of_peaks, upstream = 3000, downstream = 3000,
  type = "start_site", by = "gene", TxDb =)`, then `plotPeakProf()`.
- Overlap with reference sets: `enrichPeakOverlap(query, targets, TxDb =, method = "fisher",
  background = <universe>)`.
  - **The background is required.** Use open chromatin (DNase/ATAC) of the same cell type, or the
    union of a large reference collection (dozens of sets). Without one, peaks are shuffled over
    whole chromosomes and nearly everything is "significant" (false-positive rate 0.996 in the
    package's evaluation; 1181 of 1200 ENCODE pairs significant).
  - A universe made of little more than the sets you test is too small: every unit is hit by
    something, so p-values sit near 1 (the worked example shows this). Rank, do not test, in that case.
  - **Rank by the effect size** named in `attr(result, "params")$rank_by` (`jaccard`), never by
    p-value: the smallest p-value is shared by most sets. Report `odds_ratio` and `p.adjust` next to it.
  - The test is directional (query in target); `N_OL` is the same both ways.

## 3. Correctness checks (built in; read them before using any result)
- `annoChecks(anno)$build`: `suspect = TRUE` means peaks run past the chromosome ends of the TxDb,
  i.e. another build. **Stop.** On 116 files it raised no false alarm and found 132 of 136
  mismatched pairs; below 200 peaks (`reliable = FALSE`) it can miss a mismatch, so say the check is
  unreliable for small sets.
- `annoChecks(anno)$dropped`: peaks not annotated, with the reason per sequence ("named
  differently", "not in TxDb", "no feature"). Report the number; more than a few percent needs
  an explanation.
- `annoSource(anno)`: genome, TxDb source and date, annoDb version: copy it into the report.
- Every warning raised by an epiSeeker call is a finding. The script records each as a WARN row.

## 4. interpret() last
`clusterProfiler::interpret(anno, context =, assay = "TF", question =, overlap = <the
enrichPeakOverlap result>, model =)`. It builds a digest without a model (top regions, promoter
genes, overlap ranked by effect size, distribution, distances, widths, and the annoChecks), makes
one structured call and records it in `clusterProfiler::provenance(res)`.
- When the question is the identity of the factor, do not put the answer in `context`.
- A call counts only when `clusterProfiler::provenance(res)$ok` is TRUE; otherwise the fields are NA.
- Pass no model and the script skips this step and says so. Never call a model the user did not
  configure.

## 5. Provenance and outputs
`--out` gets `figures/`, `tables/` (annotation, annotation_stat, overlap), `report/` (report.md
with the check table; interpretation.txt and evidence_digest.txt when interpret() ran) and
`provenance/` (run.txt: arguments, input md5s, package versions, notes, warnings, output md5s;
checks.tsv; sessionInfo.txt; interpret_provenance.txt; objects.rds).

## Do not say
- "X regulates / activates gene Y" from binding near Y; binding is not regulation.
- "the target gene" for the nearest gene; it is the nearest gene.
- "the overlap proves / confirms the factor"; promoter-rich reference sets overlap any promoter-rich query.
- Any p-value from an overlap test without a background, or a ranking by p-value.
- Anything about genes or distributions when the build check fired.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/range-tf/scripts/run_range_tf.R", package = "epiSeeker"))')
Rscript "$S" --example --out range_tf_example            # bundled AR ChIP-seq peaks, no model
Rscript "$S" --peaks peaks.narrowPeak --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
  --annodb org.Hs.eg.db --context "ChIP-seq, HepG2" --targets refs/ --background dnase.bed \
  --out out --model openai:gpt-4o-mini
```
Exit status 1 with "STOP at ..." means a check fired; read `report/report.md` and
`provenance/checks.tsv`, fix the input, run again. Worked example: `references/worked-example.md`.
