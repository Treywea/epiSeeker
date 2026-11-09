# base-singlemolecule worked examples (simulated Fiber-seq)

```sh
S=$(Rscript -e 'cat(system.file("skills/base-singlemolecule/scripts/run_base_singlemolecule.R", package = "epiSeeker"))')
Rscript "$S" --example --out sm_example            [--model provider:model]
Rscript "$S" --example_control --out sm_control    [--model provider:model]
```

**Data.** A fibertools `ft extract --all` table **simulated** over chr1:1,410,000-1,440,000: 200
reads of 12-25 kb; the 18 elements are the LNCaP DNase peaks of the window (`extdata/encode/`,
merged within 100 bp). Each read has a labelling efficiency z ~ N(0, 1): an element is open (an
MSP over it) with probability plogis(logit(base) + 1.5 z), linker MSPs elsewhere occur at a rate
exp(0.8 z). Given the read, elements open **independently**: any co-occurrence is the confound.
The control run sets z = 0 for every read. Not real data.

**Input checks.** 200 reads, 0 unplaced calls; types m6A, nuc, msp; 18 of 18 elements spanned by
reads; elements' seqlevels match the reads. m6A summarised to 3,749 sites (coverage = spanning
reads, fibertools rule).

**Openness and the confound** (`tables/element_openness.tsv`, `reads.tsv`): the median element is
open on 47.3% of spanning reads. Per-read covariate = fraction of the read outside the elements
covered by MSPs (median 0.095, IQR 0.062-0.153). Median Spearman(open at an element, covariate)
= **0.48** -> WARN "per-read openness tracks element openness". Control: 0.01 -> PASS.

**Co-occurrence** (`tables/cooccurrence.tsv`, 102 pairs 1-10 kb apart, all testable):

| run | OR > 1 and p < 0.05, unadjusted (`modCooccurrence()`) | adjusted for the covariate (`modCooccurrence(adjust = "read")`) |
|---|---|---|
| example (reads differ) | 73 of 102 | 7 of 102 |
| control (equal reads) | 0 of 102 | 0 of 102 |

Every unadjusted "co-open" pair here is an artefact of reads that are open everywhere. The
adjustment removes most of it; 7 remain because the covariate measures labelling with noise
(linker MSPs are a sample of it), so adjusted results are better, not perfect: compare with
distance-matched control pairs too. `figures/cooccurrence_adjusted_vs_unadjusted.pdf`,
`reads_region.pdf`.

**interpret()** (last; only with `--model`) on the m6A levels; co-occurrence is not in its digest.

**What can be said.** "Unadjusted, 73 of 102 element pairs look co-open; adjusted for per-read
openness, 7." **Not:** "the elements open in a coordinated way", or any biology - the reads are
simulated. (On real COLO829BL Fiber-seq, MYC window: 7 of 44 FIRE pairs unadjusted, 1 of 44
adjusted; `checks/20261107_adapter-single-molecule` in the paper records.)
