# report worked example (three small runs, one of them stopped)

```sh
S=$(Rscript -e 'cat(system.file("skills/report/scripts/run_report.R", package = "epiSeeker"))')
Rscript "$S" --example --out report_example
```

**What it does.** Runs, without a model, `base-methylation --example`, `base-singlemolecule
--example`, and `base-methylation` on the bundled table with its levels multiplied by 100 (a
planted error), into `report_example/runs/`; then collects the three folders. About 6 minutes,
most of it the two analysis runs.

**Input checks.** All three folders are skill runs (PASS); every output's md5 is as recorded in its
`run.txt` (PASS).

**Runs** (`tables/runs.tsv`):

| run | skill | status | PASS | NOTE | WARN | STOP | interpret() |
|---|---|---|---|---|---|---|---|
| base-methylation | base-methylation | completed | 6 | 4 | 0 | 0 | skipped (no model) |
| base-singlemolecule | base-singlemolecule | completed | 8 | 4 | 1 | 0 | skipped (no model) |
| base-methylation_percent | base-methylation | STOPPED by a check | 1 | 1 | 0 | 1 | not run |

**Cross-run checks.** One epiSeeker version, one genome (hg38) (PASS); "runs stopped by a check":
**WARN** naming base-methylation_percent - its analysis stopped, so nothing of it is in the report;
no incomplete or flagged interpretation (PASS).

**Needs attention** (`tables/attention.tsv`): the per-read openness WARN of base-singlemolecule
(unadjusted co-occurrence confounded) and the STOP of the percent run ("Methylation 0-100; looks
like percent: divide by 100"). Then each run's own report follows, its headings demoted, its check
table replaced by this list. `figures/check_status.pdf` shows the checks per run.

**What can be said.** "Two of the three runs completed; the third stopped because its levels were
in percent." **Not:** "the third sample showed no methylation", or any number that is not in a
run's own outputs.
