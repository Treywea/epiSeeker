---
name: range-spatial
description: Turn spatial domains of spatial ATAC or CUT&Tag spots (a fragments file, a barcode-to-domain table and the spot positions) into one region set per domain with the epiSeeker adapters, after checking that the domains are not an artefact - that sequencing depth, chip rows or chip columns do not predict them (the check of the paper's Fig 2) - then annotate and compare the domains, and finish with interpret() and provenance. Stops when the domains are explained by depth. Use when a user has spatial chromatin data with domains or clusters of spots and asks how the domains differ or what they are.
---

# range-spatial: spatial domains with epiSeeker

Domains of spatial chromatin data are made upstream (LSI / clustering of spots), and they can follow
**sequencing depth** or the **rows and columns of the chip** instead of tissue structure. In the
paper's Fig 2 check (spatial ATAC of P21 mouse brain, GSM6043256), depth predicted earlier domains
with balanced accuracy 0.73-0.89, and even after removing depth-correlated components depth still
explained R^2 = 0.52 of the domains. So this skill checks first, then annotates. The script
`scripts/run_range_spatial.R` runs the steps below; prefer it. Shared steps:
`../shared/group-steps.R`.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`) and equal to `GenomeInfoDb::genome(TxDb)`.
2. **Fragments, groups and positions** all given; seqlevels of the fragments match the TxDb.
3. **Barcodes**: >= 90% of the group table found in the fragments (the `-1` suffix is handled), and
   >= 90% of in-tissue spots (`readSpotPositions()`) have a domain.
4. **Domain size**: fewer than 50 spots or 200 regions is flagged (noisy pseudo-bulk, unreliable
   build check).

## 2. Analysis (the artefact check first, then epiSeeker functions only)
- Artefact check, before any annotation. Per spot: fragments (depth), chip row, chip column.
  Balanced accuracy of LDA predicting the domain (5-fold x 5 repeats, folds stratified by domain,
  as Fig 2) from log10 depth, from the row (one-hot), from the column (one-hot); kappa =
  (accuracy - chance) / (1 - chance); R^2 of log10 depth ~ domain (`tables/artefact_check.tsv`,
  `figures/artefact_check.pdf`, `spatial_depth_log10.pdf` next to `spatial_groups.pdf`). With
  `--spot_scores` (per-spot biology, e.g. marker-module or gene-activity scores; csv, barcode
  first) the same for them.
- `readFragments()`, `fragmentsToGroups(fragments, groups, tile =, min_fragments =, fold =)`,
  `annotateGroups(groups, TxDb =, annoDb =, resolution = "spatial", coords =)` -> `csAnnoList`;
  `groupInfo()`, `plotAnnoBar()`, `plotDistToTSS()`, `plotSpatialGroups()`.

## 3. Correctness checks
- **Depth**: kappa >= 0.4 or R^2 >= 0.4 -> **STOP** ("domains explained by sequencing depth");
  re-cluster with depth removed (drop depth-correlated components, downsample spots) or pass
  `--allow_depth_confound` to go on, and then no anatomical reading. 0.15-0.4 -> WARN.
- **Rows / columns**: kappa >= 0.4 -> WARN: a channel artefact (DBiT-type chips) or banded
  anatomy; compare with the histology. A domain with >= 50% of its spots in one row or column
  (>= 10 spots) -> WARN (stripe).
- **Spot scores** (if given): WARN unless they predict the domain >= 0.05 better than depth, rows
  and columns each (Fig 2 rule SP3).
- `annoChecks()` of every domain (build suspect -> stop); every warning of an epiSeeker call as a
  WARN row.

## 4. interpret() last
`clusterProfiler::interpret(x, context =, assay =, question =, model =)` on the `csAnnoList`; the
script adds the artefact-check result to `context`. Experimental method; skipped without a model.

## 5. Provenance and outputs
`figures/` (spatial_groups, spatial_depth_log10, artefact_check, annotation_bar, distance_to_tss),
`tables/` (spots, artefact_check, groups, annotation_stat), `report/report.md`, `provenance/`
(run.txt, checks.tsv, sessionInfo.txt, interpret_provenance.txt, objects.rds; the simulated inputs).

## Do not say
- That a domain is an anatomical structure because it is spatially coherent: stripes and depth
  gradients are coherent too. Only histology or independent markers name a structure.
- Anything anatomical when the depth check warned or stopped.
- "every spot / each cell ...": a domain's regions are a pseudo-bulk.
- "domain A regulates / drives ..." from accessibility.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/range-spatial/scripts/run_range_spatial.R", package = "epiSeeker"))')
Rscript "$S" --example --out spatial_example     # simulated spots, depth independent of domain
Rscript "$S" --fragments fragments.tsv.gz --groups spot_domains.csv --positions tissue_positions_list.csv \
  --genome mm10 --txdb TxDb.Mmusculus.UCSC.mm10.knownGene --assay "spatial ATAC" --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
