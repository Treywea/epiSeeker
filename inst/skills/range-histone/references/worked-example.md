# range-histone worked example (ENCODE data, chr1-2)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-histone/scripts/run_range_histone.R", package = "epiSeeker"))')
Rscript "$S" --example --out histone_example [--model provider:model]
```

**Data.** H3K27me3 "replicated peaks" of adult prostate gland tissue (ENCODE ENCSR690CSD, file
ENCFF198QRW, GRCh38), chr1 and chr2 only: 2,385 peaks (`extdata/encode/`, source and licence in
`SOURCE.txt`). Declared: `--genome hg38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`,
assay "H3K27me3 (broad histone mark)", merge gap 5,000 bp.

**Input checks.** Genome = TxDb (PASS); seqlevels UCSC (PASS); no genome recorded in the bed
(NOTE). Median peak width 354 bp - inside the range of TF peaks (NOTE: the ENCODE histone pipeline
calls peaks in narrow mode, so the domains are built by merging).

**Domains** (`tables/domains_summary.tsv`): 2,385 peaks -> 1,275 domains within 5 kb; width
median 363 bp (IQR 234-1,802, max 41,562), 50 domains of >= 10 kb, 2.47 Mb covered. Many domains
are single peaks: the summary depends on the gap, which is why the report states it.
810 genes have a domain over the promoter (+/-3 kb) or the body (`seq2gene()`,
`tables/genes_under_domains.tsv`); 33.2% of domains are annotated to a promoter.
`annoChecks()`: 0 domains beyond chromosome ends, 0 dropped (>= 200 domains, so reliable).
Figures: `annotation_bar.pdf`, `distance_to_tss.pdf`, `tss_profile_wide.pdf` (+/-10 kb),
`gene_body_profile.pdf` (scaled bodies, 20% flanks).

**interpret()** (last; only with `--model`) on the annotated domains; the context says they are
domains merged within 5,000 bp, not binding sites. The digest starts "Regions: 1275".

**What can be said.** "On chr1-2, prostate H3K27me3 covers 2.5 Mb in 1,275 domains (merged within
5 kb), a third of them at promoters; 810 genes lie under a domain." **Not:** "H3K27me3 represses
these 810 genes", any motif or binding site, or a domain count without the gap.
