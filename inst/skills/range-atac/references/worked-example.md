# range-atac worked examples (ENCODE data, chr1-2; simulated fragments)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-atac/scripts/run_range_atac.R", package = "epiSeeker"))')
Rscript "$S" --example --out atac_example                       [--model provider:model]
Rscript "$S" --example_fragments --out atac_fragments_example   [--model provider:model]
```

**Data.** DNase-seq peaks of LNCaP after 1 nM R1881 for 12 h (ENCODE ENCSR000EKU, file
ENCFF570MWX, GRCh38), chr1 and chr2 only: 19,077 peaks (`extdata/encode/`, `SOURCE.txt`). The
second run adds 60,000 fragments **simulated** from those peaks (60% start in a peak; lengths from
a 40-150 bp and a 180-250 bp mode; `provenance/example_fragments.tsv.gz`).

**Input checks.** Genome = TxDb (PASS); seqlevels of peaks and fragments UCSC (PASS); median peak
width 180 bp (PASS).

**Accessibility** (`tables/annotation_stat.tsv`): 47.8% of peaks within 3 kb of a TSS (37.3%
within 1 kb), 32.6% intronic, 13.3% distal intergenic; median |distance to TSS| 4,478 bp.
`annoChecks()`: 0 peaks beyond chromosome ends, 0 dropped. Figures: `annotation_bar.pdf`,
`distance_to_tss.pdf`, `tss_profile.pdf`.

**Fragments** (second run, `tables/fragments_qc.tsv`): TSS enrichment 8.66 (PASS next to ENCODE's
> 7 "ideal" for ATAC on GRCh38), FRiP 0.605 (PASS next to 0.3), median fragment 132 bp, 57.7%
below 147 bp; `figures/tss_insertions.pdf`. These numbers describe the simulation (fragments were
placed in peaks, half of which are at promoters), not a library.

**interpret()** (last; only with `--model`) on the annotated peaks, assay "DNase-seq".

**What can be said.** "Half of the LNCaP DNase peaks on chr1-2 lie within 3 kb of a TSS; the rest
are intronic or distal." **Not:** "AR binds these sites", "these are active enhancers", or a QC
verdict on the simulated fragments.
