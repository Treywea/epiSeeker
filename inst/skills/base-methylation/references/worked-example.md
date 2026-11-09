# base-methylation worked examples (package data)

```sh
S=$(Rscript -e 'cat(system.file("skills/base-methylation/scripts/run_base_methylation.R", package = "epiSeeker"))')
Rscript "$S" --example --out bm_example              [--model provider:model]
Rscript "$S" --example_reads --out bm_reads_example  [--model provider:model]
```

## 1. Bulk levels (`--example`)
**Data.** `inst/extdata/demo_bisseq.txt`: 50 CpGs on chr22:10525991-10526342 of one sample
(`acinar_methyl`, GEO GSE223153), columns chr, pos, Cov, Methylation - the data of
`data(demo_bmdata)`. Declared hg38; TxDb hg38 knownGene; BSgenome hg38 for the motif profile.

**Input checks.** Level assay `Methylation` in [0, 1] (PASS), coverage assay `Cov` (3-38 reads):
34% of sites below 5 (NOTE: `interpret()` leaves them out of levels). Seqlevels UCSC (PASS).

**Analysis.** The 50 sites form 1 cluster (gap 500 bp), annotated as an intron of an
uncharacterised transcript (`tables/clusters_annotation.tsv`); `getBmMatrix()` + `plotBmProf()`
draw the CG profile of the region (`figures/profile_region.pdf`).

**Correctness.** One cluster: the build check is unreliable (NOTE, fewer than 200 regions); nothing
was dropped.

**interpret()** sends: 50 sites in 1 cluster, mean level 0.63 (sites with coverage >= 5), feature
Intron, then the level distribution and coverage (`report/evidence_digest.txt`). One sample: no
differences between samples to report.

**What can be said.** "The 50 CpGs of this region (one intronic cluster) have a mean level of 0.63
among sites with at least 5 reads; a third of the sites have fewer." **Not:** "the gene is
silenced", or any comparison - there is one sample.

## 2. Per-read calls (`--example_reads`)
**Data.** A `modkit extract full` table **simulated** from the same levels (seed 1): 40 reads across
the region, each CpG called modified with the probability given by its bulk level; written to
`provenance/example_modkit_extract.tsv`. It exercises the per-read path; it is not nanopore data.

**Input checks.** 0 calls without a reference position; call type 5mC present; coverage counted as
calls at the site because modkit lists every tested base (NOTE recording this).

**Analysis.** `summarizeModCalls()` gives a bmData with 50 sites, 40 reads each
(`tables/sites.tsv`); `plotModReads()` draws the reads (`figures/reads_region.pdf`);
`modCooccurrence()` of the two halves of the region: 40 reads span both, but the test is NA
(`tables/cooccurrence.tsv`) - in the first half every read counts as modified (0.61-0.75 of its
CpGs modified, above the 0.5 cutoff), so a margin of the 2x2 table is empty. The script reports
this as "not testable", not as "no association". The test adjusted for each read's level
outside the two halves (`adjust = "read"`, columns `adj_odds_ratio`, `adj_p_value`) is NA for the
same reason.

**What can be said.** "The per-read path runs from the table to a bmData, a read plot and a
co-occurrence test; whether the two halves are modified on the same molecules cannot be tested here."
