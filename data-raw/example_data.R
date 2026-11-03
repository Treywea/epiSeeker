#------------------------------------------ChIP-seq example------------------------------------------#
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
set.seed(929)
# data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6418464
data <- readPeakFile("./GSM6418464_H441_FOXA2_Rep1.narrowPeak.bed.gz")
demo_peak_list <- lapply(paste0("chr", 1:22), function(chr) {
    chr_data <- data[seqnames(data) == chr]
    n <- length(chr_data)
    if (n > 10) {
        chr_data[sample(seq_len(n), 10)]
    } else {
        chr_data
    }
})
demo_peak <- do.call(c, demo_peak_list)
peakAnno <- annotateSeq(demo_peak, TxDb = txdb)
seq2gene_result <- seq2gene(demo_peak, tssRegion = c(-1000, 1000), flankDistance = 3000, txdb)

write.table(demo_peak[1:50], file = "./inst/extdata/demo_peak.txt", quote = FALSE, sep = "\t")

usethis::use_data(peakAnno, overwrite = TRUE, compress = "xz")
usethis::use_data(demo_peak, overwrite = TRUE, compress = "xz")
usethis::use_data(seq2gene_result, overwrite = TRUE, compress = "xz")

#------------------------------------------Bisulfite-Seq example------------------------------------------#
# Bisulfite-Seq from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6940395, build by hg38
bisseq <- data.table::fread("./GSM6940395_fresh_acinar_methyl_1.bismark.cov.gz", data.table = FALSE)


# select a windows
demo_bisseq <- bisseq[bisseq$V1 == "22" & bisseq$V2 >= 10525991 & bisseq$V3 <= 10526342, ]
demo_bisseq$V1 <- paste0("chr", demo_bisseq$V1)
demo_bisseq$Cov <- as.numeric(demo_bisseq$V5) + as.numeric(demo_bisseq$V6)
demo_bisseq <- demo_bisseq[, c(1, 2, 7, 4)]
colnames(demo_bisseq) <- c("chr", "pos", "Cov", "Methylation")
demo_bisseq$Methylation <- demo_bisseq$Methylation * 0.01

demo_bmdata <- makeBmDataFromData(data = list(acinar_methyl = demo_bisseq), sampleNames = "acinar_methyl")
usethis::use_data(demo_bmdata, overwrite = TRUE, compress = "xz")
write.table(demo_bisseq, file = "./inst/extdata/demo_bisseq.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#------------------------------------------Human ref motif------------------------------------------#
library(RSQLite)
library(TFBSTools)
opts_base <- list()
opts_base[["collection"]] <- "CORE"
opts_base[["all_versions"]] <- FALSE
opts_human <- opts_base
opts_human[["species"]] <- "Homo sapiens"
opts_human[["tax_group"]] <- "vertebrates"
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024::JASPAR2024()))
pwm_obj <- TFBSTools::getMatrixSet(sq24, opts_human)

usethis::use_data(pwm_obj, overwrite = TRUE, compress = "xz")

#------------------------------------------tagMatrix result------------------------------------------#
## the GEO example files are hg38 since 2026-11-03 (lifted from hg19 with
## data-raw/lift_sample_data_hg38.R); run that first
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile(getSampleFiles()[[4]])
tagMatrix <- getTagMatrix(peak,
    type = "start_site", by = "gene",
    upstream = 3000, downstream = 3000,
    TxDb = txdb, weightCol = "V5", nbin = 500
)
usethis::use_data(tagMatrix, overwrite = TRUE, compress = "xz")

#------------------------------------------peakAnnoList------------------------------------------#
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
files <- getSampleFiles()
peakAnnoList <- lapply(files, annotateSeq,
    TxDb = txdb, annoDb = "org.Hs.eg.db",
    tssRegion = c(-3000, 3000), verbose = FALSE
)
names(peakAnnoList) <- names(files)
usethis::use_data(peakAnnoList, overwrite = TRUE, compress = "xz")
