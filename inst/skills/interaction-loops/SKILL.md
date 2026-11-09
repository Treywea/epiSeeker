---
name: interaction-loops
description: Annotate chromatin loops (Hi-C, HiChIP, ChIA-PET BEDPE) and link peaks to genes through them with the epiSeeker R package - loop classes (promoter-promoter, promoter-distal, distal-distal), spans, genes at promoter anchors, and for each peak whether loops point to its nearest gene or to others - with checks for build and seqlevels mismatches and peaks annotated with another TxDb, then interpret() and provenance. Use when a user asks what their loops connect, which genes distal peaks reach through loops, or how loop links compare with the nearest gene. Loop partners are reported as candidate genes, never as targets.
---

# interaction-loops: loops and peak-to-gene links with epiSeeker

The script `scripts/run_interaction_loops.R` runs the steps below; prefer it over re-typing them.
Normalisation, loop calling and contact matrices belong to upstream tools (HiCExperiment, Juicer,
cooler); this skill starts from called loops.

## The rule of this skill: a loop partner is a candidate, not the target
In the epiSeeker paper's test (K562 Pol II ChIA-PET loops against 9,830 CRISPRi-tested
enhancer-gene pairs, 450 of them regulated), genes reached through real loops were regulated
**no more often** than genes reached through loops rewired at random with the same anchors and
the same distances: precision 0.220 vs 0.219 (enrichment 1.00, p = 0.47), in every distance bin
and in every loop subset. Loops did reach more tested pairs: recall of nearest gene or loop
0.740 vs 0.567 for the nearest gene alone. So loops widen the list of candidate genes; they do not
tell which candidate is regulated. Write "loop-linked candidate gene", never "target gene", and
never rank candidates as if a loop made one more likely than another gene at the same distance.

## 1. Input checks (stop on failure)
1. **Genome build** declared (`--genome`), equal to `GenomeInfoDb::genome(TxDb)`; loops and peaks
   must be on that build.
2. **Seqlevels** of loop anchors and peaks share a style with the TxDb (`--fix_seqlevels` converts,
   recorded).
3. **BEDPE**: `readLoops()` reads 0-based BEDPE (also gzipped, Juicer `#` headers) and counts rows
   without coordinates (`metadata(gi)$loopSource`).
4. **Gene set**: count only the genes you mean. The paper used protein-coding genes
   (`gene_ids` from `org.*.eg.db` `GENETYPE == "protein-coding"`, `--protein_coding`); with all
   genes, LOC/LINC/snoRNA promoters on wide anchors inflate the links.

## 2. Analysis (epiSeeker functions only)
- `annotateLoops(loops, TxDb =, annoDb =, gene_ids =)` -> a `loopAnno`: anchors annotated with
  `annotateSeq()`, promoter anchors and their genes, loop class and span; `plotLoopAnno(type =
  "class" / "span")`.
- `annotateSeq(peaks, TxDb = <the same TxDb>)`, then `linkPeaksToGenes(peakAnno, loopAnno,
  gene_ids =)`: one row per peak with its nearest gene, the genes at the other ends of its loops and
  a class ("promoter", "no loop", "loop, no gene", "loop includes nearest", "loop, other genes");
  `plotLoopLinks()`. `--strict` uses the paper's definitions (`anchor_flank = 1000`,
  `distal_distance = 3000`, `exclude_own_anchor = TRUE`).

## 3. Correctness checks
- `annoChecks(loopAnno@anchors)` and `annoChecks(peakAnno)`: build suspect -> **stop**; dropped
  anchors or peaks -> report.
- `annoSource()` of the peaks and of the anchors must name the same TxDb and genome
  (`linkPeaksToGenes()` assumes it).
- Loops between chromosomes and loops whose anchors are closer than two anchor widths are counted
  and reported.

## 4. interpret() last
`clusterProfiler::interpret(loopAnno, context =, question =, links = <linkPeaksToGenes result>,
model =)`. The digest lists genes at promoter anchors by number of loops, the peak link classes,
the genes reached only through loops, loop classes, anchor distribution and spans, and ends with
a fixed caveat: loop partners are candidate genes (the Fig 4 rewiring result above, 0.220 vs
0.219), so the rule reaches the model also when interpret() is called directly. The method is
**experimental** (no comparison with a no-model baseline). The script flags (does not edit) any
"target", "regulates" or "proves" in the answer.

## 5. Provenance and outputs
`figures/` (loop_classes, loop_span, peak_links), `tables/` (loops_annotated, loop_classes,
peak_links), `report/report.md` (with the rule above), `provenance/` (run.txt, checks.tsv,
sessionInfo.txt, interpret_provenance.txt, objects.rds; the toy BEDPE of the example).

## Do not say
- "gene X is the target of this enhancer / peak" because a loop joins them.
- "the loop shows that the peak regulates X", or that loop links are more precise than distance.
- "the nearest gene is wrong" because a loop points elsewhere; both are candidates.
- Anything about the toy example's genes: its loops are built from peak positions and carry no biology.

## Run it
```sh
S=$(Rscript -e 'cat(system.file("skills/interaction-loops/scripts/run_interaction_loops.R", package = "epiSeeker"))')
Rscript "$S" --example --out loops_example
Rscript "$S" --loops loops.bedpe --peaks atac.bed --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene --annodb org.Hs.eg.db --protein_coding --strict --out out
```
Exit status 1 with "STOP at ..." means a check fired. Worked example: `references/worked-example.md`.
