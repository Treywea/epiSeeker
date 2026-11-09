# interaction-loops worked example (package data, toy loops)

```sh
S=$(Rscript -e 'cat(system.file("skills/interaction-loops/scripts/run_interaction_loops.R", package = "epiSeeker"))')
Rscript "$S" --example --out loops_example [--model provider:model]
```

**Data.** The package ships no loops, so the example builds **toy** loops from the bundled peaks:
each AR peak (ARmo_1nM, 2288 peaks) is joined to the nearest CBX7 peak 20 kb - 1 Mb away on the same
chromosome, anchors 5 kb (`provenance/example_loops.bedpe`). They exercise the code; their
genes mean nothing. Peaks: the same ARmo_1nM set. Settings as in the paper's Fig 4: protein-coding
genes only (20,599 of org.Hs.eg.db), `--strict` (`anchor_flank = 1000`, `distal_distance = 3000`,
`exclude_own_anchor = TRUE`).

**Input checks.** Genome hg38 = TxDb (PASS); no genome recorded on the BEDPE or the bed (NOTE);
seqlevels UCSC for anchors and peaks (PASS); 1157 loops read, 0 rows dropped.

**Analysis.** `annotateLoops()`: 1157 loops between 1823 anchors, median span 390 kb;
promoter-promoter 149 (12.9%), promoter-distal 664 (57.4%), distal-distal 344 (29.7%)
(`tables/loop_classes.tsv`). `linkPeaksToGenes()` for the 2288 peaks: promoter 274, no loop 1016,
loop without a gene 410, loop includes the nearest gene 110, loop to other genes 478
(`tables/peak_links.tsv`).

**Correctness.** `annoChecks()` of anchors and peaks: no peak beyond chromosome ends, none dropped;
peaks and anchors annotated with the same TxDb (PASS); 0 loops between chromosomes.

**interpret()** (last). The digest lists genes at promoter anchors by number of loops (651 genes),
the link classes, the genes reached only through loops, classes, anchor distribution and spans,
and ends with its fixed caveat ("Genes at the other end of a loop are candidate genes, not
established targets ... precision 0.220 vs 0.219 ...").

**What can be said.** "478 of the 2014 peaks outside promoters reach, through a loop, genes other
than their nearest gene; these loop-linked genes are candidates." **Not:** "TNFSF11 is the target
of the peak at chr13:42924211", "the loop shows that the peak regulates ...", or anything about
these genes at all - the loops are toys.

The phrase flags are tested with a mock answer that contains "target gene", "regulates" and
"proves": the run completes and `provenance/checks.tsv` carries a WARN "forbidden statements in
the interpretation -- loop partner called the target: ...".
