# epiSeeker agent skills

Instructions and scripts that let an AI agent run epiSeeker analyses the way the package is
meant to be used: the same functions a person calls, in the right order, with the package's
correctness checks read before any result is used, and `interpret()` as the last step.

| skill | data | adapters / classes |
|---|---|---|
| `range-tf` | narrow peaks: TF ChIP-seq, CUT&RUN | `csAnno`, `csAnnoList`, `peakOverlapResult` |
| `range-histone` | broad histone marks, merged into domains | `csAnno`, `seq2gene()`, `getTagMatrix(type = "body")` |
| `range-atac` | ATAC / DNase peaks, fragments | `csAnno`, `readFragments()`, `getTagMatrix()` |
| `range-singlecell` | scATAC clusters | `fragmentsToGroups()`, `peaksToGroups()`, `csAnnoList` |
| `range-spatial` | spatial domains, checked against depth / chip rows / columns | `fragmentsToGroups()`, `readSpotPositions()`, `csAnnoList` |
| `base-methylation` | single-base levels; per-read 5mC / m6A tables | `bmData`, `readModCalls()` |
| `base-singlemolecule` | per-read Fiber-seq / modkit tables; co-occurrence adjusted for per-read openness | `readModCalls()`, `modReadMatrix()`, `modCooccurrence()` |
| `interaction-loops` | chromatin loops (BEDPE) and peaks | `loopAnno`, `linkPeaksToGenes()` |
| `report` | the folders of other skill runs | (collects; no analysis, no model) |

Every skill is a folder in the Agent Skills format: `SKILL.md` (YAML front matter with `name` and
`description`, then the instructions), `scripts/` (one runnable R script) and
`references/worked-example.md` (what the script gives on the package's own example data). All
`SKILL.md` follow one structure: input checks (genome build, seqlevels, assay) -> analysis with
epiSeeker functions -> built-in correctness checks -> `interpret()` last -> provenance, then the
statements the skill does not allow. `shared/skill-helpers.R` holds the bookkeeping the scripts
share (arguments, check table, warnings, provenance); `shared/range-steps.R` (region sets) and
`shared/group-steps.R` (groups of cells or spots) hold the steps several skills run the same way.
None of them has analysis of its own. The examples of `range-tf`, `range-histone`, `range-atac`
and `base-singlemolecule` use two small ENCODE region sets in `extdata/encode/` (source and
licence in `SOURCE.txt`).

Use them from the installed package:
```sh
Rscript -e 'cat(system.file("skills", package = "epiSeeker"))'   # where they are
cp -r "$(Rscript -e 'cat(system.file("skills", package = "epiSeeker"))')"/range-tf ~/.claude/skills/
```
The scripts load their helpers with `system.file()`, so a copied skill folder still runs against
the installed package. Each script takes `--example` (package data, no model needed) and writes
`figures/`, `tables/`, `report/` and `provenance/` under `--out`. `interpret()` runs only when a
model is given (`--model provider:model` or `EPISEEKER_MODEL`), through aisdk.
