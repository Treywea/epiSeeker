## Lift the hg19 example peak files to hg38 (2026-11-03).
##
## The GEO example files (GSE40740 CBX6/CBX7, GSE48308 ARmo) and inst/extdata/sample_peaks.txt
## are hg19 peaks, while the vignette and the tests annotate with TxDb.Hsapiens.UCSC.hg38.knownGene:
## CBX6_BF had 58% of its peaks in promoters (<=1 kb) against hg19 and 15% against hg38, and
## its Reactome enrichment disappeared. The files are lifted with the UCSC chain so that every
## example is hg38; the hg38 examples already shipped (demo_peak, demo_bmdata) are unchanged.
##
## A peak is kept only when it maps to one range of the same width on the same chromosome, so
## that column offsets such as the MACS summit stay valid; the others are dropped and counted.
## The order of the peaks and every other column are kept.
##
## Run from the package root: Rscript data-raw/lift_sample_data_hg38.R
## Then rebuild tagMatrix and peakAnnoList with data-raw/example_data.R.

library(rtracklayer)
library(GenomicRanges)

chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
## set EPISEEKER_CHAIN_DIR to keep the chain between runs; a temporary directory otherwise
chain_gz <- file.path(Sys.getenv("EPISEEKER_CHAIN_DIR", tempdir()), basename(chain_url))
if (!file.exists(chain_gz)) utils::download.file(chain_url, chain_gz, mode = "wb", quiet = TRUE)
chain_file <- sub("\\.gz$", "", chain_gz)
if (!file.exists(chain_file)) system2("gunzip", c("-k", chain_gz))
chain <- import.chain(chain_file)

## readPeakFile() adds 1 to the second column of every file, so the range it reads is
## [col2 + 1, col3]; that range is lifted and col2 is written back as start - 1.
lift_table <- function(df) {
    gr <- GRanges(df[[1]], IRanges(df[[2]] + 1, df[[3]]))
    hits <- liftOver(gr, chain)
    one <- lengths(hits) == 1
    to <- rep(NA_integer_, length(gr)); s <- e <- to; chr <- rep(NA_character_, length(gr))
    u <- unlist(hits[one])
    s[one] <- start(u); e[one] <- end(u); chr[one] <- as.character(seqnames(u))
    keep <- one & !is.na(s) & (e - s) == (df[[3]] - df[[2]] - 1) & chr == as.character(df[[1]])
    out <- df[keep, , drop = FALSE]
    out[[1]] <- chr[keep]
    out[[2]] <- s[keep] - 1L
    out[[3]] <- e[keep]
    list(table = out, n = nrow(df), unmapped = sum(lengths(hits) == 0),
         split = sum(lengths(hits) > 1), changed = sum(one & !keep))
}

report <- function(name, r) {
    cat(sprintf("%-55s %5d peaks: kept %5d, unmapped %d, split %d, width or chromosome changed %d\n",
                name, r$n, nrow(r$table), r$unmapped, r$split, r$changed))
}

## the five GEO bed files (no header)
for (f in list.files("inst/extdata/GEO_sample_data", full.names = TRUE)) {
    df <- utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE)
    r <- lift_table(df)
    report(basename(f), r)
    con <- gzfile(f, "w")
    utils::write.table(r$table, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
}

## sample_peaks.txt (MACS table with a header; `length` = end - start + 1 as before)
f <- "inst/extdata/sample_peaks.txt"
df <- utils::read.delim(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
r <- lift_table(df)
report(basename(f), r)
utils::write.table(r$table, f, sep = "\t", quote = FALSE, row.names = FALSE)
