# range-tf worked example (package data)

```sh
S=$(Rscript -e 'cat(system.file("skills/range-tf/scripts/run_range_tf.R", package = "epiSeeker"))')
Rscript "$S" --example --out range_tf_example [--model provider:model]
```

**Data.** The bundled GEO peak sets (`getSampleFiles()`), lifted from hg19 to hg38 for the package:
AR ChIP-seq in prostate cancer cells at 0, 1 and 100 nM androgen (ARmo_0M, ARmo_1nM, ARmo_100nM;
GSM1174480-2) and CBX6 / CBX7 (GSM1295076-7). Query: ARmo_1nM; compared with ARmo_0M and
ARmo_100nM; overlap tested against the other four sets. Declared: `--genome hg38`,
`TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`, assay TF, question "Which transcription
factor was profiled?" (the context names the cell line and treatment, not the factor).

**Input checks.** Declared genome = TxDb genome (PASS). No genome recorded in the bed files
(NOTE: the out-of-bound check is what protects against another build). Seqlevels UCSC (PASS).
Median peak width 333.5 / 86 / 147 bp, inside or below the 184-556 bp of 20 ENCODE TF sets (PASS).

**Annotation** (`tables/annotation_stat.tsv`): ARmo_1nM 2288 peaks, 17.8% within 3 kb of a TSS
(8.9% within 1 kb), 34.2% other introns, 23.3% distal intergenic; ARmo_0M 800 peaks, 16.2%;
ARmo_100nM 1351 peaks, 17.8%. `annoChecks()`: 0 peaks beyond chromosome ends and 0 dropped for all
three (each >= 200 peaks, so the build check is reliable).

**Overlap** (`tables/overlap.tsv`, `method = "fisher"`, ranked by `jaccard`):

| target | N_OL | jaccard | odds_ratio | p.adjust |
|---|---|---|---|---|
| ARmo_100nM | 655 | 0.136 | 0.422 | 1 |
| ARmo_0M | 25 | 0.00324 | 0.515 | 1 |
| CBX7_BF | 8 | 0.000569 | 0.0259 | 1 |
| CBX6_BF | 4 | 0.000331 | 0.149 | 1 |

The only background the package ships is the union of these five sets, and the script flags it
(WARN: 100% of the 5,758 background units overlap the query or a target). In such a universe every
unit is hit by some set, so the test has nothing to compare with: all p.adjust are 1 and the odds
ratios fall below 1. This is the conservative failure, the mirror image of the whole-genome null.
The ranking by effect size still works - the AR set of the same cells (ARmo_100nM) is first by a
factor of 40 in jaccard - and that is all this example supports. With real data, pass
open chromatin of the cell type or the union of a large reference collection as `--background`.

**interpret()** (last; only with `--model`). The digest it sends starts with the overview (2288
regions, hg38, "Checks: genome build consistent with TxDb", declared assay TF, the context and the
question), then the top 200 regions (ranked by distance to TSS: the bed files carry no signal
column), promoter genes, the overlap table above ranked by jaccard, distribution, distances and
widths (`report/evidence_digest.txt`).

**What can be said.** "The peaks are mostly intronic and distal (17.8% within 3 kb of a TSS); they
overlap the AR peaks of the same cells at 100 nM far more than CBX6/CBX7 peaks (jaccard 0.136 vs
< 0.001), consistent with an AR ChIP; the bundled background is too small for p-values."
**Not:** "AR regulates <nearest genes>", "the overlap proves the factor is AR", any p-value.
