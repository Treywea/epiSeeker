---
name: base-methylation
description: Analyse single-base modification data with the epiSeeker R package - bisulfite or EM-seq methylation tables, modkit pileup bedMethyl files, and single-molecule per-read tables (modkit extract for nanopore/PacBio 5mC, fibertools for Fiber-seq m6A) summarised to a bmData - with checks for the silent errors of such data (levels in percent instead of fractions, low coverage, another genome build, strand handling of per-read calls), then interpret() and provenance. Use when a user asks what methylation or base modification looks like in their data, around which genes, how samples differ, or what the reads of one locus show.
---

# base-methylation: single-base levels and per-read calls with epiSeeker

Every number comes from an epiSeeker function. The script `scripts/run_base_methylation.R`
runs the steps below for one input; prefer it over re-typing them. Not covered yet: calling
differentially methylated regions (DMRs) - the package has no DMR caller, so do not hand-roll one
and present it as an epiSeeker result.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`); when a TxDb is used it must match
   `GenomeInfoDb::genome(TxDb)`, and a BSgenome for motif context must be of the same build.
2. **Seqlevels** of the sites share a style with the TxDb.
3. **Level scale.** Levels must be fractions in [0, 1]. Percentages (0-100, as many tools write)
   are the commonest silent error: `interpret()` would pick the wrong assay or stop, plots would be
   off by 100. Divide by 100 first.
4. **Coverage.** A coverage assay should be present; sites below `min_cov` (default 5) are left out
   of levels and differences by `interpret()`. Report the fraction below it.
5. **Per-read tables.** `readModCalls()` reads `modkit extract full` / `extract calls` and
   `ft extract --all`; calls without a reference position are dropped and counted
   (`metadata(x)$modSource$n_unplaced`). `summarizeModCalls(x, type = "5mC")` counts coverage as
   calls at the site when the tool lists every tested base (modkit), and as reads spanning the site
   when it lists only calls above a threshold (fibertools; for 5mC only reads of the strand that
   gave calls there, since PacBio calls CpGs on the read's own strand). Record which applied.

## 2. Analysis (epiSeeker functions only)
- Build the `bmData`: `makeBmDataFromData(<chr, pos, coverage, level table>, sampleNames =)`,
  `readBedMethyl(file, mod = "m")`, or `summarizeModCalls(readModCalls(file), type =)`; a saved
  bmData (`data(demo_bmdata)`) works directly.
- Where the sites are: merge sites within 500 bp into clusters and `annotateSeq(clusters, TxDb =,
  annoDb =)`; `plotAnnoBar()`.
- A locus in sequence context: `getBmMatrix(region = data.frame(chr, start, end), input = bm,
  BSgenome =, base = "C", motif = "CG")`, then `plotBmProf()`.
- Per read (single-molecule data only): `plotModReads(reads, region, types = "5mC")`;
  `modReadMatrix(reads, regions, type =)`; `modCooccurrence(reads, region1, region2, type =)`
  for whether two regions are modified on the same molecules (needs reads spanning both).

## 3. Correctness checks
- `annoChecks()` of the cluster annotation (build, dropped clusters). With fewer than 200 clusters
  the build check is unreliable; say so.
- The level-scale and coverage checks above, recorded as rows of `provenance/checks.tsv`.
- `modCooccurrence()` returns `NA` when fewer than 2 reads span both regions: report "not
  testable", not "no association".
- Every warning of an epiSeeker call is recorded as a WARN row.

## 4. interpret() last
`clusterProfiler::interpret(bm, assay =, context =, question =, value = <level assay>,
coverage = <coverage assay>, TxDb =, annoDb =, min_cov = 5, model =)`. The digest lists
clusters with mean levels per sample (and feature / nearest gene with a TxDb), the sites with the
largest differences between samples, the level distribution and the coverage. The method is
**experimental**: it has not been compared with a baseline that uses no model. Without a model the
script skips the step and says so.

## 5. Provenance and outputs
`figures/` (clusters_annotation_bar, profile_region, reads_region), `tables/` (sites,
clusters_annotation, cooccurrence), `report/report.md`, `provenance/` (run.txt with input md5s,
versions, warnings, output md5s; checks.tsv; sessionInfo.txt; interpret_provenance.txt;
objects.rds; the simulated input of the per-read example).

## Do not say
- "the gene is silenced / repressed" from methylation levels alone; levels are not expression.
- "significant(ly) different" without a test on replicates; the digest's differences are
  descriptive (it counts sites with |difference| >= 0.2).
- Anything about co-occurrence on molecules from bulk levels; only per-read data can show it.
- A DMR call; the package does not make them yet.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/base-methylation/scripts/run_base_methylation.R", package = "epiSeeker"))')
Rscript "$S" --example --out bm_example              # bundled 50 CpGs on chr22
Rscript "$S" --example_reads --out bm_reads_example  # per-read path, simulated from them
Rscript "$S" --input sample.bed --format bedmethyl --mod m --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db --assay "ONT 5mC" --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
