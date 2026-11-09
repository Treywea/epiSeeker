---
name: range-singlecell
description: Turn clusters of single cells (scATAC, single-cell CUT&Tag) or spatial domains (spatial ATAC / CUT&Tag spots) into one region set per group with the epiSeeker adapters (fragmentsToGroups from a fragments file and a barcode-to-group table, or peaksToGroups from ArchR / Signac / Cell Ranger marker tables), then annotate and compare the groups like bulk samples, with checks for barcode matching, small groups, build and seqlevels mismatches, and finish with interpret() and provenance. Use when a user has a fragments file or per-cluster peaks and asks how clusters or spatial domains differ in where their open chromatin or marks lie.
---

# range-singlecell: groups of cells or spots with epiSeeker

The adapters turn each group into a pseudo-bulk region set; after that, annotation, comparison
and `interpret()` are the bulk functions. The script `scripts/run_range_singlecell.R` runs the
steps below; prefer it. Clustering cells or spots is upstream (Signac, ArchR, SnapATAC2, a spatial
pipeline); bring the barcode-to-group table.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`) and equal to `GenomeInfoDb::genome(TxDb)`.
2. **Seqlevels** of the fragments share a style with the TxDb.
3. **Barcodes.** `fragmentsToGroups()` matches barcodes exactly, or without a `-1`-like suffix when
   fewer than half match, and reports how many matched (`metadata(groups)$match`). Fewer than 90%
   of the group table found (the adapters' own evaluation rule) means the wrong table, sample or
   suffix: find out before going on. For spatial data, the in-tissue spots of
   `readSpotPositions()` should have a group.
4. **Group size.** A group of few cells is a noisy pseudo-bulk; fewer than 50 members (the
   adapters' evaluation rule) or 200 regions is flagged, and below 200 regions the build check is
   unreliable.
5. **Tile width** fits the assay: about 500 bp for ATAC, about 5000 bp for a broad histone mark.

## 2. Analysis (epiSeeker functions only)
- `readFragments(file)`; `fragmentsToGroups(fragments, groups, tile =, min_fragments = 5,
  specific = TRUE, fold = 2)`: per group, tiles with enough fragments and enriched over the other
  groups. Or `peaksToGroups(marker_table)` for per-cluster peaks (`ArchR::getMarkers()`,
  `Signac::FindMarkers()`, Cell Ranger ATAC differential tables).
- `annotateGroups(groups, TxDb =, annoDb =, resolution = "single-cell" | "spatial", coords =)` ->
  a `csAnnoList`; `groupInfo()` gives resolution, members and regions per group.
- `plotAnnoBar(x)`, `plotDistToTSS(x)`; spatial: `plotSpatialGroups(coords)`.

## 3. Correctness checks
- `annoChecks()` of every group (build suspect -> **stop**; dropped regions -> report).
- The barcode, group-size and spatial-coverage checks above, as rows of `provenance/checks.tsv`.
- Every warning of an epiSeeker call is recorded as a WARN row.

## 4. interpret() last
`clusterProfiler::interpret(x, context =, assay =, question =, model =)` on the `csAnnoList`.
The digest states the resolution and members of each group, then compares distributions, widths,
distances to TSS, shared and group-specific promoter genes and top regions. The `csAnnoList`
method is **experimental** (not evaluated against a baseline). Without a model the script skips it.

## 5. Provenance and outputs
`figures/` (annotation_bar, distance_to_tss, spatial_groups), `tables/` (groups,
annotation_stat), `report/report.md`, `provenance/` (run.txt, checks.tsv, sessionInfo.txt,
interpret_provenance.txt, objects.rds; the simulated inputs of the examples).

## Do not say
- Anything about "every cell" or "each spot": a group's regions are a pseudo-bulk.
- A cell type named from the regions alone, as fact; it is a hypothesis to check with markers.
- "cluster A regulates / drives ..." from accessibility.
- Anything biological about the examples: their fragments are simulated from bundled peaks.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/range-singlecell/scripts/run_range_singlecell.R", package = "epiSeeker"))')
Rscript "$S" --example --out sc_example                 # simulated scATAC, 3 groups
Rscript "$S" --example_spatial --out spatial_example    # the same cells as spots on a grid
Rscript "$S" --fragments fragments.tsv.gz --groups clusters.csv --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db --assay ATAC --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
