# epiSeeker agent skills

Instructions and scripts that let an AI agent run epiSeeker analyses the way the package is
meant to be used: the same functions a person calls, in the right order, with the package's
correctness checks read before any result is used, and `interpret()` as the last step.

| skill | data | adapters / classes |
|---|---|---|
| `range-tf` | narrow peaks: TF ChIP-seq, CUT&RUN, ATAC summits | `csAnno`, `csAnnoList`, `peakOverlapResult` |
| `base-methylation` | single-base levels; per-read 5mC / m6A tables | `bmData`, `readModCalls()` |
| `interaction-loops` | chromatin loops (BEDPE) and peaks | `loopAnno`, `linkPeaksToGenes()` |
| `range-singlecell` | scATAC clusters, spatial domains | `fragmentsToGroups()`, `peaksToGroups()`, `csAnnoList` |

Every skill is a folder in the Agent Skills format: `SKILL.md` (YAML front matter with `name` and
`description`, then the instructions), `scripts/` (one runnable R script) and
`references/worked-example.md` (what the script gives on the package's own example data). All
`SKILL.md` follow one structure: input checks (genome build, seqlevels, assay) -> analysis with
epiSeeker functions -> built-in correctness checks -> `interpret()` last -> provenance, then the
statements the skill does not allow. `shared/skill-helpers.R` holds the bookkeeping the scripts
share (arguments, check table, warnings, provenance); it has no analysis in it.

Use them from the installed package:
```sh
Rscript -e 'cat(system.file("skills", package = "epiSeeker"))'   # where they are
cp -r "$(Rscript -e 'cat(system.file("skills", package = "epiSeeker"))')"/range-tf ~/.claude/skills/
```
The scripts load their helpers with `system.file()`, so a copied skill folder still runs against
the installed package. Each script takes `--example` (package data, no model needed) and writes
`figures/`, `tables/`, `report/` and `provenance/` under `--out`. `interpret()` runs only when a
model is given (`--model provider:model` or `EPISEEKER_MODEL`), through aisdk.
