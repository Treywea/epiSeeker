---
name: range-histone
description: Analyse broad histone-mark peak sets (H3K27me3, H3K36me3, H3K9me3, H3K79me2 and other marks that cover kilobases) with the epiSeeker R package - merge peaks into domains, summarise domains (number, width, bp covered, genes under them), annotate them, profile them over wide windows and gene bodies, optionally test their overlap with reference sets in a background - with checks for the silent errors (wrong genome build, seqlevels that do not match, narrow-mode peak calls of a broad mark), then interpret() and provenance. Use when a user has broad-mark ChIP-seq or CUT&Tag peaks or domains and asks which genes or regions the mark covers or how two such sets differ. Not for TF peaks (range-tf) or ATAC (range-atac).
---

# range-histone: broad histone marks with epiSeeker

A broad mark covers domains of kilobases, so the unit of analysis is the **domain**, not the peak:
"nearest gene of a peak", motifs and binding sites, which make sense for a TF, do not here. The
script `scripts/run_range_histone.R` runs the steps below; prefer it. The steps it shares with
`range-tf` and `range-atac` are in `../shared/range-steps.R`.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`) and equal to `GenomeInfoDb::genome(TxDb)`; TxDb named
   (`annotateSeq()` with `TxDb = NULL` falls back to hg38).
2. **Seqlevels** of the peaks share a style with the TxDb ("1" vs "chr1" drops every peak);
   convert only on purpose (`--fix_seqlevels`, recorded).
3. **Peak width.** A median width inside the range of TF peaks (184-556 bp in 20 ENCODE TF sets)
   means the peaks were called in narrow mode (as the ENCODE histone pipeline does): say so; the
   domains below are built by merging. Broad-mode calls (`.broadPeak`, SICER, epic2) need less
   merging.

## 2. Analysis (epiSeeker functions only)
- `readPeakFile()`; domains = peaks merged within `--merge_gap` bp (default 5000;
  `GenomicRanges::reduce(min.gapwidth =)`). Report the gap: domain counts depend on it.
- Domain summary: number, width quartiles, number >= 10 kb, Mb covered (`tables/domains_summary.tsv`).
- `annotateSeq(domains, TxDb =, annoDb =)`; `plotAnnoBar()`, `plotDistToTSS()` of the domains.
- Genes under domains: `seq2gene(domains, tssRegion =, flankDistance = 0, TxDb =)` - genes whose
  promoter or body a domain overlaps (`tables/genes_under_domains.tsv`).
- Wide profiles: `getTagMatrix(domains, upstream = 10000, downstream = 10000, type = "start_site",
  by = "gene", TxDb =)`; gene bodies scaled: `getTagMatrix(type = "body", nbin = 100,
  upstream = rel(0.2), downstream = rel(0.2))`; both with `plotPeakProf()`. Promoter-centred marks
  (H3K27me3 at CpG-island promoters, H3K4me3) and gene-body marks (H3K36me3, H3K79me2) differ here.
- Overlap with reference sets: as in `range-tf` (`enrichPeakOverlap(method = "fisher",
  background =)`, background required, ranked by `jaccard`), with the domains as query.

## 3. Correctness checks
- `annoChecks()` of the domains: build suspect -> **stop**; dropped domains -> report.
- The width check above; every warning of an epiSeeker call as a WARN row.
- A gene under a domain is a gene the mark covers. Expression is not measured here.

## 4. interpret() last
`clusterProfiler::interpret(<csAnno of the domains>, context =, assay = "H3K27me3", question =,
model =)`. The script adds to `context` that the regions are domains merged from peaks, not binding
sites. Without a model the step is skipped and the report says so. A call counts only when
`clusterProfiler::provenance(res)$ok` is TRUE.

## 5. Provenance and outputs
`figures/` (annotation_bar, distance_to_tss, tss_profile_wide, gene_body_profile), `tables/`
(domains_summary, domains_<set>, annotation_<set>, annotation_stat, genes_under_domains, overlap),
`report/report.md`, `provenance/` (run.txt, checks.tsv, sessionInfo.txt, interpret_provenance.txt,
objects.rds).

## Do not say
- Motifs, binding sites, "X binds": a histone mark is not a factor (flagged as "TF-style claim").
- "the mark represses / silences / activates gene Y": coverage is not regulation.
- "the target gene" of a domain; it is a gene under (or nearest to) the domain.
- Domain counts or widths without the merge gap they depend on.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/range-histone/scripts/run_range_histone.R", package = "epiSeeker"))')
Rscript "$S" --example --out histone_example     # ENCODE prostate H3K27me3, chr1-2, no model
Rscript "$S" --peaks H3K36me3.broadPeak --genome hg38 --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
  --annodb org.Hs.eg.db --assay H3K36me3 --merge_gap 5000 --out out --model openai:gpt-4o-mini
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
