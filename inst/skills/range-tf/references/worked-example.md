# range-tf worked example (package data)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-tf/scripts/run_range_tf.R", package = "epiSeeker"))')
Rscript "$S" --example --out range_tf_example [--model provider:model]
```

**Data.** The bundled GEO peak sets (`getSampleFiles()`), lifted from hg19 to hg38 for the package:
AR ChIP-seq in LNCaP-ARmo prostate cancer cells at 0, 1 and 100 nM DHT (ARmo_0M, ARmo_1nM,
ARmo_100nM; GSM1174480-2) and CBX6 / CBX7 (GSM1295076-7). Query: ARmo_1nM; compared with ARmo_0M
and ARmo_100nM; overlap tested against the other four sets in a background of DNase-seq peaks of
LNCaP after 1 nM R1881 (ENCODE ENCFF570MWX, chr1-chr2 only; `extdata/encode/SOURCE.txt`).
Declared: `--genome hg38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`, assay TF, question
"Which transcription factor was profiled?" (the context names the cell line and treatment, not the
factor).

**Input checks.** Declared genome = TxDb genome (PASS). No genome recorded in the bed files
(NOTE: the out-of-bound check is what protects against another build). Seqlevels UCSC (PASS).
Median peak width 333.5 / 86 / 147 bp, inside or below the 184-556 bp of 20 ENCODE TF sets (PASS).

**Annotation** (`tables/annotation_stat.tsv`): ARmo_1nM 2288 peaks, 17.8% within 3 kb of a TSS
(8.9% within 1 kb), 34.2% other introns, 23.3% distal intergenic; ARmo_0M 800 peaks, 16.2%;
ARmo_100nM 1351 peaks, 17.8%. `annoChecks()`: 0 peaks beyond chromosome ends and 0 dropped for all
three (each >= 200 peaks, so the build check is reliable).

**Overlap** (`tables/overlap.tsv`, `method = "fisher"`, ranked by `jaccard`). The background covers
chr1-chr2, so the test runs there (NOTE: 365 of 2,288 query peaks, 826 of 5,118 target peaks);
3% of the 19,077 background regions overlap the query or a target (PASS: the universe is broader
than the tested sets).

| target | qLen | tLen | N_OL | jaccard | odds_ratio | p.adjust |
|---|---|---|---|---|---|---|
| GSM1174482_ARmo_100nM_peaks | 365 | 217 | 123 | 0.156 | 113 | < 1e-300 |
| GSM1174480_ARmo_0M_peaks | 365 | 105 | 3 | 0.0023 | 34.0 | 1.5e-07 |
| GSM1295077_CBX7_BF_ChipSeq_mergedReps_peaks | 365 | 278 | 2 | 0.00054 | 0.35 | 0.99 |
| GSM1295076_CBX6_BF_ChipSeq_mergedReps_peaks | 365 | 226 | 1 | 0.00035 | 0.18 | 0.99 |

With open chromatin of the same cells as the universe the test has something to compare with: the
AR peaks of the same cells at 100 nM overlap the query far beyond chance (odds ratio 113), the
3 shared peaks with ARmo_0M are also more than expected (0.2), and the CBX6/CBX7 sets overlap less
than expected. (With the union of the five bundled sets as universe, batch 1 of the skills got
p.adjust = 1 for all four: every unit was hit by some set.) Rank by `jaccard`, report `odds_ratio`
and `p.adjust` next to it; a p-value of 0 is a floating-point underflow, not a measurement.

**interpret()** (last; only with `--model`). The digest it sends starts with the overview (2288
regions, hg38, "Checks: genome build consistent with TxDb", declared assay TF, the context and the
question), then the top 200 regions (ranked by distance to TSS: the bed files carry no signal
column), promoter genes, the overlap table above ranked by jaccard, distribution, distances and
widths (`report/evidence_digest.txt`).

**What can be said.** "The peaks are mostly intronic and distal (17.8% within 3 kb of a TSS); on
chr1-2, within LNCaP open chromatin, they overlap the AR peaks of the same cells at 100 nM far more
than expected (odds ratio 113, jaccard 0.156) and CBX6/CBX7 peaks less than expected, consistent
with an AR ChIP." **Not:** "AR regulates <nearest genes>", "the overlap proves the factor is AR",
a ranking by p-value, or a genome-wide statement from the chr1-2 test.
