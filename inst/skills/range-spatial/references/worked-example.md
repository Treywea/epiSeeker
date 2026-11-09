# range-spatial worked example (simulated spots)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-spatial/scripts/run_range_spatial.R", package = "epiSeeker"))')
Rscript "$S" --example --out spatial_example [--model provider:model]
```

**Data.** The simulated cells of `range-singlecell` (fragments from the bundled peaks; groups
AR_like 60, PcG_like 60, AR_0M_like 20) as spots on a 14 x 10 grid, with 100-300 fragments kept
per spot at random (independent of the domain). Domains are concentric: AR_0M_like the 20 spots
nearest the centre, AR_like the next 60, PcG_like the outer 60
(`provenance/example_fragments.tsv.gz`, `example_groups.csv`, `example_tissue_positions_list.csv`).
Not real data.

**Input checks.** Genome = TxDb; 140 of 140 barcodes matched (without `-1`); 100% of in-tissue
spots have a domain. AR_0M_like (20 spots, 170 regions) and PcG_like (86 regions) are flagged as
small (NOTE): with 100-300 fragments per spot, fewer tiles pass `min_fragments`.

**Artefact check** (`tables/artefact_check.tsv`, `figures/artefact_check.pdf`):

| predictor | balanced accuracy | chance | kappa |
|---|---|---|---|
| depth | 0.363 | 0.333 | 0.05 |
| row | 0.381 | 0.333 | 0.07 |
| column | 0.616 | 0.333 | 0.42 |

R^2 of log10 depth ~ domain 0.005. Depth: PASS (domains not explained by depth). Columns: **WARN**
(kappa 0.42 >= 0.4): on a 14-column grid the outer domain owns the edge columns, so the column
alone predicts it. That is what the check is for - it cannot tell banded anatomy from a channel
artefact, the histology has to. No stripe (largest share of a domain in one row or column 30%).

**Domains** (`tables/groups.tsv`, `annotation_stat.tsv`): AR_0M_like 170 regions, 18.8% within
3 kb of a TSS; AR_like 402, 15.7%; PcG_like 86, 59.3%. Figures: `spatial_groups.pdf`,
`spatial_depth_log10.pdf`, `annotation_bar.pdf`, `distance_to_tss.pdf`.

**interpret()** (last; only with `--model`): the context carries the depth-check result.

**What can be said.** "Depth does not predict the domains (kappa 0.05); the outer domain's tiles
are promoter-rich (59.3%), the inner two mostly distal." **Not:** "the outer domain is a tissue
layer", or anything biological - the spots are simulated. On the paper's real spatial ATAC data
(Fig 2, GSM6043256) the same check stops the run: depth explains the domains.
