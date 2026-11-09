# range-singlecell worked example (package data, simulated cells)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-singlecell/scripts/run_range_singlecell.R", package = "epiSeeker"))')
Rscript "$S" --example --out sc_example                 [--model provider:model]
```

**Data.** A fragments file **simulated** from the bundled peaks (seed 1): 140 cells, 300 fragments
each, 70% inside the peaks of their group and 30% anywhere on chr1-22. Groups: `AR_like` (60 cells,
ARmo_1nM peaks), `PcG_like` (60, CBX7 peaks), `AR_0M_like` (20, ARmo_0M peaks). Fragment barcodes
end in `-1`, the group table's do not (`provenance/example_fragments.tsv.gz`,
`example_groups.csv`). Not real data. (The same cells as spatial spots are the example of
`range-spatial`.)

**Input checks.** Genome hg38 = TxDb; fragments in UCSC style; 140 of 140 barcodes matched after
dropping the `-1` suffix (PASS, >= 90%).
`AR_0M_like` is flagged as small (NOTE: 20 cells, below the 50 of the adapters' evaluation rule).

**Analysis.** `fragmentsToGroups(tile = 500, min_fragments = 5, fold = 2)` then
`annotateGroups()` (`tables/groups.tsv`, `annotation_stat.tsv`):

| group | cells | regions | % within 3 kb of a TSS |
|---|---|---|---|
| AR_0M_like | 20 | 438 | 18.5 |
| AR_like | 60 | 1131 | 15.9 |
| PcG_like | 60 | 443 | 58.7 |

The groups recover what they were simulated from: the CBX7-like group is promoter-rich, the AR-like
groups are not. Figures: `annotation_bar.pdf`, `distance_to_tss.pdf`.

**Correctness.** `annoChecks()` of every group: no region beyond chromosome ends, none dropped;
every group has >= 200 regions, so the build check is reliable.

**interpret()** (last) on the `csAnnoList`: the digest begins "Resolution: single-cell;
pseudo-bulk of AR_0M_like 20 cells, AR_like 60 cells, PcG_like 60 cells", then compares the groups.

**What can be said.** "PcG_like's enriched tiles are promoter-rich (58.7% within 3 kb of a TSS), the
two AR-like groups' are mostly distal (15.9-18.5%); AR_0M_like is a small group." **Not:** "every
cell of PcG_like ...", a cell type as fact, or any biology - the cells are simulated.
