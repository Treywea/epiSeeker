---
name: report
description: Collect the outputs and provenance of several epiSeeker skill runs (range-tf, range-histone, range-atac, range-singlecell, range-spatial, base-methylation, base-singlemolecule, interaction-loops) into one report - a table of runs, every STOP and WARN of their checks, the interpret() status of each, cross-run consistency (one genome, one epiSeeker version) and an md5 check that each run's outputs are unchanged since the run - without calling a model or computing new numbers. Use when a user wants one document for an analysis made of several skill runs, or to see at a glance what in those runs needs attention.
---

# report: one report from several skill runs

Every other skill writes `report/report.md` and `provenance/` into its own folder. This skill reads
those folders and puts them in one report, checks first. It computes nothing new and calls no
model. The script `scripts/run_report.R` does it; prefer it.

## 1. Input checks (stop on failure)
1. Each folder given (`--runs a,b,c` or `--runs_dir`) **is a skill run**: it has
   `provenance/run.txt` and `provenance/checks.tsv`; otherwise stop.
2. **Outputs unchanged since the run**: the md5 of every output listed in `run.txt` is recomputed;
   a file changed or missing -> WARN (results edited by hand, or a later run wrote into the folder).

## 2. Analysis (collection only)
- `tables/runs.tsv`: per run the skill, status (completed / STOPPED), counts of PASS / NOTE / WARN /
  STOP, interpret() status, epiSeeker version, genome, outputs changed, seconds.
- `tables/checks_all.tsv` (every check row of every run) and `tables/attention.tsv` (STOP and WARN).
- `figures/check_status.pdf`: checks per run by status.
- `report/report.md`: the run table, every STOP and WARN, then each run's own report (its check
  table replaced by the attention list) with the paths of its figures, tables and interpretation.

## 3. Correctness checks
- **Runs stopped by a check** -> WARN naming them: their results are absent, not negative.
- **Runs declare one genome** and **used one epiSeeker version** -> WARN otherwise (TxDb
  differences -> NOTE).
- **interpret() answers complete** (none INCOMPLETE) and **interpretations without flagged
  statements** (the phrase flags of the runs) -> WARN otherwise.

## 4. interpret() last
Not called here. Each run's interpret() answer, its evidence digest and its provenance
(`clusterProfiler::provenance()`) stay in that run's folder; the report links them and states their
status. Summarising several runs with a model would be a new, unchecked analysis.

## 5. Provenance and outputs
`provenance/run.txt` lists every collected `run.txt` with its md5, so the report can be traced back
to the exact runs; `checks.tsv`, `sessionInfo.txt`, `objects.rds` as in the other skills.

## Do not say
- Any number not in the runs' own outputs; no pooled p-value or combined ranking across runs.
- That a stopped run found "nothing"; it was not analysed.
- A run's result without its WARN rows next to it.
- That outputs are a run's results when the md5 check says they changed.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/report/scripts/run_report.R", package = "epiSeeker"))')
Rscript "$S" --example --out report_example      # runs 3 small skill examples (one stops), then collects them
Rscript "$S" --runs out/range_tf,out/loops,out/methylation --out out/report
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
