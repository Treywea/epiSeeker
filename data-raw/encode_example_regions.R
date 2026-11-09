## Make inst/extdata/encode/*.bed.gz: chr1-chr2 subsets of two ENCODE GRCh38 files, used by the
## examples of the agent skills (see inst/extdata/encode/SOURCE.txt for source and licence).
## Run from the package root with network access.
files <- c(LNCaP_R1881_DNase_ENCFF570MWX = "ENCFF570MWX",      # ENCSR000EKU, md5 11fc2eb70efc494a00217350f0e83035
           prostate_H3K27me3_ENCFF198QRW = "ENCFF198QRW")      # ENCSR690CSD, md5 db7ed4b10bdf4f8eb4a7e37d14e50ea0
dir.create("inst/extdata/encode", showWarnings = FALSE, recursive = TRUE)
for (nm in names(files)) {
    acc <- files[[nm]]
    tmp <- tempfile(fileext = ".bed.gz")
    utils::download.file(sprintf("https://www.encodeproject.org/files/%s/@@download/%s.bed.gz", acc, acc), tmp, mode = "wb")
    message(acc, " md5 ", tools::md5sum(tmp))
    d <- utils::read.delim(tmp, header = FALSE)[, 1:3]
    d <- d[d$V1 %in% c("chr1", "chr2"), ]
    d <- unique(d[order(match(d$V1, c("chr1", "chr2")), d$V2, d$V3), ])
    con <- gzfile(file.path("inst/extdata/encode", paste0(nm, "_chr1-2.bed.gz")), "w", compression = 9)
    utils::write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
}
