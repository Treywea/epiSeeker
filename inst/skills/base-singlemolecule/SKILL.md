---
name: base-singlemolecule
description: Analyse single-molecule chromatin data from per-read tables - fibertools ft extract --all for Fiber-seq (m6A, nucleosomes, MSPs) or modkit extract for nanopore / PacBio 5mC and 6mA - with the epiSeeker R package - per-element openness on each read, pairwise co-occurrence of elements on the same molecules, read plots and a bmData summary - with checks for the silent errors (unplaced calls, missing call types, element and read seqlevels that do not match) and, above all, the per-read openness confound - co-occurrence between elements is inflated by how open or how well labelled each read is overall, so the skill adjusts for it - then interpret() and provenance. Use when a user asks whether elements open together on the same molecules, what reads of a locus look like, or how open each element is per molecule.
---

# base-singlemolecule: per-read chromatin with epiSeeker

Per-read data can show what bulk data cannot: whether two elements are open **on the same
molecule**. But reads differ in how open (or how well m6A-labelled) they are everywhere, and that
alone makes any two elements look co-open. On real Fiber-seq (COLO829BL, MYC window) 7 of 44 FIRE
element pairs were "co-open" (OR > 1, p < 0.05), 1 of 44 after adjusting for each read's overall
MSP fraction, and 23% of distance-matched non-FIRE control pairs were "co-open" too. The script
`scripts/run_base_singlemolecule.R` reports co-occurrence adjusted for this; prefer it.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`); with a TxDb, equal to `GenomeInfoDb::genome(TxDb)`.
2. **Per-read table**: `readModCalls(file, format =)` reads `ft extract --all` and `modkit extract
   full|calls`; calls without a reference position are dropped and counted
   (`metadata(x)$modSource$n_unplaced`). The call type analysed (`--type`, default `msp` for
   fibertools) must be present.
3. **Elements** (`--elements`, e.g. FIRE peaks, DNase / ATAC peaks): their seqlevels style must
   match the reads; elements no read spans are dropped (WARN when fewer than half are spanned).

## 2. Analysis (epiSeeker functions only)
- `summarizeModCalls(reads, type =)` -> `bmData` (5mC if present, else m6A); coverage counted as
  calls (modkit) or spanning reads (fibertools).
- Elements: `annotateSeq()` with a TxDb; per read and element `modReadMatrix(reads, elements,
  type =)`; an element is open on a read when the value is >= 0.5 (`tables/element_openness.tsv`).
- Per-read covariate (`tables/reads.tsv`): for `msp` / `nuc`, the fraction of the read **outside
  the elements** covered by that type; for base types, the fraction (modkit) or per-kb count
  (fibertools) of modified calls outside the elements.
- Pairs of elements 1-10 kb apart: `modCooccurrence(reads, e1, e2, type =, adjust = "read",
  method = "logistic", covariate = <the per-read covariate, named by read id>)` - the unadjusted
  Fisher test (`odds_ratio`, `p_value`) and, on the same reads, `glm(open2 ~ open1 + covariate,
  binomial)` (`adj_odds_ratio`, `adj_p_value`) - `tables/cooccurrence.tsv`.
- **Distance-matched control pairs**: the adjustment is only as good as its covariate (a noisy one
  leaves part of the confound, a strongly predictive one over-adjusts). Before calling pairs
  co-open, run the same test on control pairs at the same distances away from the elements and
  report their positive rate next to that of the pairs.
- `plotModReads(reads, region, peaks = elements)`.

## 3. Correctness checks
- **Per-read openness tracks element openness**: median over elements of Spearman(open at the
  element, covariate) >= 0.1 -> WARN; then unadjusted co-occurrence is confounded and only the
  adjusted odds ratios speak to elements opening together.
- Pairs without an adjusted fit stay confounded (WARN). `modCooccurrence()` returns `NA` when
  fewer than 2 reads span both elements or a margin is empty: "not testable", not "no association".
- `annoChecks()` of the elements (fewer than 200: build check unreliable); every warning of an
  epiSeeker call as a WARN row.

## 4. interpret() last
`clusterProfiler::interpret(bm, assay =, context =, question =, TxDb =, annoDb =, model =)` on the
summarised levels; the co-occurrence table is not part of its digest. Experimental method;
skipped without a model.

## 5. Provenance and outputs
`figures/` (reads_region, cooccurrence_adjusted_vs_unadjusted), `tables/` (element_openness,
reads, cooccurrence, elements_annotation), `report/report.md`, `provenance/` (run.txt,
checks.tsv, sessionInfo.txt, interpret_provenance.txt, objects.rds; the simulated table).

## Do not say
- "elements A and B open together / are coordinated / cooperate" from unadjusted co-occurrence,
  or from the adjusted test without distance-matched control pairs.
- Any co-occurrence statement from bulk levels; only per-read data can show it.
- "the element regulates ..." from openness on reads.
- Anything biological about the example: its reads are simulated to show the confound.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/base-singlemolecule/scripts/run_base_singlemolecule.R", package = "epiSeeker"))')
Rscript "$S" --example --out sm_example                  # simulated Fiber-seq, uneven labelling
Rscript "$S" --example_control --out sm_control          # the same with equal labelling
Rscript "$S" --input sample.ft_all.tsv.gz --elements fire_peaks.bed --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db --assay Fiber-seq --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
