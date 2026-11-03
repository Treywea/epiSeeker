#' Annotate peaks
#'
#' @title annotateSeq
#' @param peak peak file or GRanges object
#' @param tssRegion Region Range of TSS
#' @param TxDb TxDb or EnsDb annotation object
#' @param level one of transcript and gene
#' @param assignGenomicAnnotation logical, assign peak genomic annotation or not
#' @param genomicAnnotationPriority genomic annotation priority
#' @param annoDb annotation package
#' @param addFlankGeneInfo logical, add flanking gene information from the peaks.
#'   The resulting `flank_gene_distances` column reports 0 when the peak
#'   overlaps the feature range (at `level = "transcript"` this is the whole
#'   transcript, so peaks inside a transcript body get 0) and the signed
#'   distance to the feature TSS otherwise.
#' @param flankDistance distance of flanking sequence
#' @param sameStrand logical, whether find nearest/overlap gene in the same strand
#' @param ignoreOverlap logical, whether ignore overlap of TSS with peak
#' @param ignoreUpstream logical, if True only annotate gene at the 3' of the peak.
#' @param ignoreDownstream logical, if True only annotate gene at the 5' of the peak.
#' @param overlap one of 'TSS' or 'all', if overlap="all", then gene overlap with peak will be reported as nearest gene, no matter the overlap is at TSS region or not.
#' @param verbose print message or not
#' @param columns names of columns to be obtained from database
#' @return data.frame or GRanges object with columns of:
#'
#' all columns provided by input.
#'
#' annotation: genomic feature of the peak, for instance if the peak is
#' located in 5'UTR, it will annotated by 5'UTR. Possible annotation is
#' Promoter-TSS, Exon, 5' UTR, 3' UTR, Intron, and Intergenic.
#'
#' geneChr: Chromosome of the nearest gene
#'
#' geneStart: gene start
#'
#' geneEnd: gene end
#'
#' geneLength: gene length
#'
#' geneStrand: gene strand
#'
#' geneId: entrezgene ID
#'
#' distanceToTSS: distance from peak to gene TSS
#'
#' if addFlankGeneInfo is TRUE, extra columns will be included:
#'
#' flank_geneIds: semicolon-separated gene IDs within `flankDistance` of the peak
#'
#' flank_gene_distances: semicolon-separated distances of these genes. A value
#' of 0 means that the peak overlaps the feature range: at level = "transcript"
#' the feature is the whole transcript, so peaks inside a transcript body
#' always get 0 regardless of their distance to the TSS. For non-overlapping
#' features the signed distance to the feature TSS is reported. See
#' https://github.com/YuLab-SMU/ChIPseeker/issues/235
#'
#' if annoDb is provided, extra column will be included:
#'
#' ENSEMBL: ensembl ID of the nearest gene
#'
#' SYMBOL: gene symbol
#'
#' GENENAME: full gene name
#' @import GenomeInfoDb
#' @importFrom methods new
#' @examples
#' data(peakAnno)
#' peakAnno
#' @seealso [plotAnnoBar()] [plotAnnoPie()] [plotDistToTSS()]
#' @export
#' @author G Yu
annotateSeq <- function(peak,
                        tssRegion = c(-3000, 3000),
                        TxDb = NULL,
                        level = "transcript",
                        assignGenomicAnnotation = TRUE,
                        genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                        annoDb = NULL,
                        addFlankGeneInfo = FALSE,
                        flankDistance = 5000,
                        sameStrand = FALSE,
                        ignoreOverlap = FALSE,
                        ignoreUpstream = FALSE,
                        ignoreDownstream = FALSE,
                        overlap = "TSS",
                        verbose = TRUE,
                        columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")) {
    is_GRanges_of_TxDb <- FALSE
    if (is(TxDb, "GRanges")) {
        is_GRanges_of_TxDb <- TRUE
        assignGenomicAnnotation <- FALSE
        annoDb <- NULL
        addFlankGeneInfo <- FALSE
        message("#\n#.. 'TxDb' is a self-defined 'GRanges' object...\n#")
        message("#.. Some parameters of 'annotateSeq' will be disable,")
        message("#.. including:")
        message("#..\tlevel, assignGenomicAnnotation, genomicAnnotationPriority,")
        message("#..\tannoDb, addFlankGeneInfo and flankDistance.")
        message("#\n#.. Some plotting functions are designed for visualizing genomic annotation")
        message("#.. and will not be available for the output object.\n#")
    }

    if (is_GRanges_of_TxDb) {
        level <- "USER_DEFINED"
    } else {
        level <- match.arg(level, c("transcript", "gene"))
    }

    if (assignGenomicAnnotation && all(genomicAnnotationPriority %in% c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")) == FALSE) {
        stop('genomicAnnotationPriority should be any order of c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")')
    }

    if (is(peak, "GRanges")) {
        ## this test will be TRUE
        ## when peak is an instance of class/subclass of "GRanges"
        input <- "gr"
        peak.gr <- peak
    } else {
        input <- "file"
        peak.gr <- loadPeak(peak, verbose)
    }

    peakNum <- length(peak.gr)

    if (verbose) {
        message(
            ">> preparing features information...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
    }

    if (is_GRanges_of_TxDb) {
        features <- TxDb
    } else {
        TxDb <- loadTxDb(TxDb)

        if (level == "transcript") {
            features <- getGene(TxDb, by = "transcript")
        } else {
            features <- getGene(TxDb, by = "gene")
        }
    }
## An empty input returns an empty annotation. The nearest feature lookup
    ## below and the flank/annotation steps after it are not defined for zero
    ## peaks and used to fail with unrelated errors ("Error: invalid
    ## subscript", "replacement has 1 row, data has 0").
    if (peakNum == 0) {
        anno <- peak.gr
        anno@seqinfo <- seqinfo(TxDb)[seqlevels(anno)]

        detail <- data.frame(
            genic = logical(0), Intergenic = logical(0),
            Promoter = logical(0), fiveUTR = logical(0),
            threeUTR = logical(0), Exon = logical(0), Intron = logical(0),
            downstream = logical(0), distal_intergenic = logical(0)
        )
        annoStat <- data.frame(
            Feature = factor(character(0)),
            Frequency = numeric(0)
        )

        return(new("csAnno",
            anno = anno,
            tssRegion = tssRegion,
            level = level,
            hasGenomicAnnotation = assignGenomicAnnotation,
            detailGenomicAnnotation = detail,
            annoStat = annoStat,
            peakNum = peakNum,
            annoSource = .make_anno_source(TxDb, annoDb)
        ))
    }


    if (verbose) {
        message(
            ">> identifying nearest features...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
    }

    ## nearest features
    idx.dist <- getNearestFeatureIndicesAndDistances(peak.gr, features,
        sameStrand, ignoreOverlap,
        ignoreUpstream, ignoreDownstream,
        overlap = overlap
    )

    if (verbose) {
        message(
            ">> calculating distance from peak to TSS...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
    }
    ## distance
    distance <- idx.dist$distance

    ## update peak, remove un-map peak if exists.
    peak.gr <- idx.dist$peak

    n.dropped <- peakNum - length(peak.gr)
    if (n.dropped > 0) {
        warning(n.dropped, " of ", peakNum, " peaks were dropped, ",
                "since no feature of 'TxDb' can be found for them. ",
                "This usually happens for peaks located on contigs/scaffolds ",
                "that carry no gene, or when the seqlevels style of the peaks ",
                "does not match the one of 'TxDb'.",
                call. = FALSE)
    }

    ## nothing left to annotate: report the most common cause instead of
    ## failing later on with an unrelated error (issue #238 of ChIPseeker)
    if (peakNum > 0 && length(peak.gr) == 0) {
        stop("all ", peakNum, " peaks were dropped, so there is nothing to ",
             "annotate. None of the seqlevels of 'peak' matches the seqlevels ",
             "of 'TxDb' (e.g. 'chr1' versus 'NC_000001.11'). Check the style ",
             "of both with GenomeInfoDb::seqlevelsStyle(seqlevels(x)) and align ",
             "them with GenomeInfoDb::seqlevelsStyle(peak) <- 'UCSC' when one ",
             "side uses '1' and the other 'chr1'; names such as ",
             "'NC_000001.11' are not mapped automatically and have to be ",
             "renamed explicitly, for example seqlevels(peak) <- ",
             "sub('^NC_0*(\\\\d+)\\..*$', 'chr\\\\1', seqlevels(peak)).",
             call. = FALSE)
    }

    ## annotation
    if (assignGenomicAnnotation == TRUE) {
        if (verbose) {
            message(
                ">> assigning genomic annotation...\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n"
            )
        }

        anno <- getGenomicAnnotation(peak.gr, distance, tssRegion, TxDb, level, genomicAnnotationPriority, sameStrand = sameStrand)
        annotation <- anno[["annotation"]]
        detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]

        ## Keep the feature level metadata tied to the feature that supplied
        ## the exon/intron annotation when genes or isoforms overlap
        ## (issue #252).
        if (level == "transcript") {
            aligned <- .align_transcript_annotation(
                peak.gr, features, idx.dist$index, distance,
                anno[["annotationFeatureId"]]
            )
        } else if (level == "gene") {
            aligned <- .align_annotation_feature(
                peak.gr, features, idx.dist$index, distance,
                anno[["annotationFeatureGene"]], idColumn = "gene_id"
            )
        } else {
            aligned <- NULL
        }

        if (!is.null(aligned)) {
            idx.dist$index <- aligned$index
            distance <- aligned$distance
        }
    } else {
        annotation <- NULL
        detailGenomicAnnotation <- NULL
    }

    ## append annotation to peak.gr
    if (!is.null(annotation)) {
        mcols(peak.gr)[["annotation"]] <- annotation
    }


    has_nearest_idx <- which(idx.dist$index <= length(features))
    nearestFeatures <- features[idx.dist$index[has_nearest_idx]]

    ## duplicated names since more than 1 peak may annotated by only 1 gene
    names(nearestFeatures) <- NULL
    nearestFeatures.df <- as.data.frame(nearestFeatures)
    if (is_GRanges_of_TxDb) {
        colnames(nearestFeatures.df)[seq_len(5)] <- c(
            "geneChr", "geneStart", "geneEnd",
            "geneLength", "geneStrand"
        )
    } else if (level == "transcript") {
        if (is(TxDb, "EnsDb")) {
            nearestFeatures.df <- nearestFeatures.df[, c(
                "seqnames", "start",
                "end", "width",
                "strand", "gene_id",
                "tx_id", "tx_biotype"
            ),
            drop = FALSE
            ]
            colnames(nearestFeatures.df) <- c(
                "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand",
                "geneId", "transcriptId", "transcriptBiotype"
            )
        } else {
            colnames(nearestFeatures.df) <- c(
                "geneChr", "geneStart", "geneEnd",
                "geneLength", "geneStrand",
                "geneId", "transcriptId"
            )
            nearestFeatures.df$geneId <- TXID2EG(
                as.character(nearestFeatures.df$geneId),
                geneIdOnly = TRUE
            )
        }
    } else {
        if (is(TxDb, "EnsDb")) {
            nearestFeatures.df <- nearestFeatures.df[, c(
                "seqnames", "start",
                "end", "width",
                "strand", "gene_id",
                "gene_biotype"
            ),
            drop = FALSE
            ]
            colnames(nearestFeatures.df) <- c(
                "geneChr", "geneStart", "geneEnd",
                "geneLength", "geneStrand",
                "geneId", "geneBiotype"
            )
        } else {
            colnames(nearestFeatures.df) <- c(
                "geneChr", "geneStart", "geneEnd",
                "geneLength", "geneStrand",
                "geneId"
            )
        }
    }

    for (cn in colnames(nearestFeatures.df)) {
        v <- nearestFeatures.df[[cn]]
        ## as.data.frame() returns 'seqnames' and 'strand' as factors; assigning
        ## a factor into mcols() drops the class and keeps only the integer
        ## codes, so that geneChr/geneStrand were reported as integers
        ## (e.g. 1/2 instead of +/-). Turn factors into characters first.
        if (is.factor(v)) v <- as.character(v)
        mcols(peak.gr)[[cn]][has_nearest_idx] <- unlist(v)
    }

    mcols(peak.gr)[["distanceToTSS"]] <- distance

    if (!is.null(annoDb)) {
        if (verbose) {
            message(
                ">> adding gene annotation...\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n"
            )
        }
        .idtype <- IDType(TxDb)
        if (length(.idtype) == 0 || is.na(.idtype) || is.null(.idtype)) {
            n <- length(peak.gr)
            if (n > 100) {
                n <- 100
            }
            sampleID <- peak.gr$geneId[seq_len(n)]

            if (all(grepl("^ENS", sampleID))) {
                .idtype <- "Ensembl Gene ID"
            } else if (all(grepl("^\\d+$", sampleID))) {
                .idtype <- "Entrez Gene ID"
            } else {
                warning("Unknown ID type, gene annotation will not be added...")
                .idtype <- NA
            }
        }

        if (!is.na(.idtype)) {
            peak.gr %<>% addGeneAnno(annoDb, .idtype, columns)
        }
    }

    if (addFlankGeneInfo == TRUE) {
        if (verbose) {
            message(
                ">> adding flank feature information from peaks...\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n"
            )
        }

        flankInfo <- getAllFlankingGene(peak.gr, features, level, flankDistance)

        if (level == "transcript") {
            mcols(peak.gr)[["flank_txIds"]] <- NA
            mcols(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
        }

        mcols(peak.gr)[["flank_geneIds"]] <- NA
        mcols(peak.gr)[["flank_gene_distances"]] <- NA

        mcols(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
        mcols(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances
    }

    if (!is_GRanges_of_TxDb) {
        if (verbose) {
            message(
                ">> assigning chromosome lengths\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n"
            )
        }

        peak.gr@seqinfo <- seqinfo(TxDb)[names(seqlengths(peak.gr))]
    }

    if (verbose) {
        message(
            ">> done...\t\t\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
    }

    if (assignGenomicAnnotation) {
        res <- new("csAnno",
            anno = peak.gr,
            tssRegion = tssRegion,
            level = level,
            hasGenomicAnnotation = TRUE,
            detailGenomicAnnotation = detailGenomicAnnotation,
            annoStat = getGenomicAnnoStat(peak.gr),
            peakNum = peakNum,
            annoSource = .make_anno_source(TxDb, annoDb)
        )
    } else {
        res <- new("csAnno",
            anno = peak.gr,
            tssRegion = tssRegion,
            level = level,
            hasGenomicAnnotation = FALSE,
            peakNum = peakNum,
            annoSource = .make_anno_source(TxDb, annoDb)
        )
    }

    return(res)
}


#' dropAnno
#'
#' drop annotation exceeding distanceToTSS_cutoff
#' @title dropAnno
#' @param csAnno output of annotateSeq
#' @param distanceToTSS_cutoff distance to TSS cutoff
#' @return csAnno object
#' @export
#' @examples
#' data(peakAnno)
#' dropAnno(peakAnno)
#' @author Guangchuang Yu
dropAnno <- function(csAnno, distanceToTSS_cutoff = 10000) {
    idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) < distanceToTSS_cutoff)
    csAnno@anno <- csAnno@anno[idx]
    csAnno@peakNum <- length(idx)
    if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- getGenomicAnnoStat(csAnno@anno)
        csAnno@detailGenomicAnnotation <- csAnno@detailGenomicAnnotation[idx, ]
    }
    csAnno
}
