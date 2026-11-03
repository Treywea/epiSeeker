#' @keywords internal
"_PACKAGE"

#' @title Information Datasets
#'
#' @description ucsc genome version, precalculated data and gsm information
#' @section Provenance:
#' The `gsminfo` dataset was constructed programmatically from public
#' resources in the NCBI GEO and UCSC Genome Browser databases.
#' The data generation pipeline is implemented in
#' `data-raw/` (see `prepareGSMInfo()` in the package source).
#'
#' Briefly, GEO metadata were retrieved using the `GEOmetadb` SQLite
#' database and `GEOquery`. The latest GEOmetadb SQLite file was downloaded
#' via `getSQLiteFile()` or, if unavailable, directly from
#' <http://starbuck1.s3.amazonaws.com/sradb/GEOmetadb.sqlite.gz>.
#' Platform (GPL) records were queried to identify platforms associated with
#' high-throughput sequencing experiments. For each sequencing platform, the
#' corresponding GSM records were obtained using `Meta(getGEO())`.
#' Supplementary BED-like files for each GSM were collected using
#' `getGSMsuppFile()` and `batchGetGSMsuppFile()`.
#'
#' Additional metadata fields (title, organism, extract protocol, characteristics,
#' data processing description, submission date, and supplementary file URLs)
#' were extracted from GSM SOFT files downloaded using `GEOquery`.
#' Genome assembly versions for each GSM were inferred using the function
#' `getGenomicVersion()`, which matches UCSC genome labels to either
#' the data processing description or the supplementary file names, using the
#' reference table provided in the internal dataset `ucsc_release`.
#'
#' PubMed IDs associated with each GEO series (GSE) were obtained from the
#' `gse` table in GEOmetadb. All GSM-level metadata were merged, cleaned,
#' and converted to ASCII using `iconv()` to remove non-ASCII characters.
#'
#' Finally, newly processed GSM entries were appended to any preexisting
#' `gsminfo` object stored in the package, deduplicated, and saved as
#' `gsminfo.rda` with `compress="xz"`.
#'
#' Thus, `gsminfo` represents a curated, reproducibly constructed metadata
#' table summarizing GEO high-throughput sequencing samples, including organism,
#' platform, experimental descriptions, processing information, genome versions,
#' supplementary BED file locations, and associated PubMed IDs.
#'
#' @name gsminfo
#' @aliases ucsc_release
#' @section Data structure:
#' A data frame with one row per GSM sample and the following columns:
#' \describe{
#'   \item{`series_id`}{GEO series accession (GSE).}
#'   \item{`gsm`}{GEO sample accession (GSM).}
#'   \item{`gpl`}{GEO platform accession (GPL).}
#'   \item{`organism`}{Organism name (e.g., *Mus musculus*).}
#'   \item{`title`}{Sample title as provided in GEO.}
#'   \item{`characteristics`}{Experiment-specific metadata such as cell type, treatment, or antibody.}
#'   \item{`source_name`}{Source material for sequencing, typically cell or tissue type.}
#'   \item{`extract_protocol`}{Detailed wet-lab protocol for chromatin extraction, immunoprecipitation, and library preparation as reported in GEO.}
#'   \item{`description`}{Antibody information or additional sample description.}
#'   \item{`data_processing`}{Bioinformatics processing description including aligner, genome build, peak calling method, and filtering steps.}
#'   \item{`submission_date`}{Date when the sample was submitted to GEO.}
#'   \item{`supplementary_file`}{URL to supplementary processed files (e.g., BED).}
#'   \item{`genomeVersion`}{Genome assembly used in the processed data (e.g., mm8, hg19).}
#'   \item{`pubmed_id`}{PMID of the reference publication associated with the dataset.}
#' }
#'
#' @format A data frame with `n` rows (GSM samples) and 14 columns.
#' @docType data
#' @keywords datasets
#' @return data frame
NULL


#' @title Example data of peak annotation
#'
#' @description A `csAnno` object representing the annotation result of the example peak set `demo_peak`.
#' Peaks were annotated using the function `annotateSeq()` in `epiSeeker`.
#' @section Provenance:
#' Input peaks were taken from the example dataset `demo_peak`.
#' Annotation was generated using `epiSeeker::annotateSeq()`.
#' @section Data structure:
#' A `csAnno` S4 object with the following slots:
#' \describe{
#'   \item{`anno`}{A `GRanges` object containing the annotated peaks,
#'          including peak coordinates, basic peak metrics, and gene-based annotation fields.}
#'   \item{`tssRegion`}{Numeric vector of length two defining the upstream
#'          and downstream window used for TSS annotation.}
#'   \item{`level`}{Character string indicating whether annotation was
#'          performed at the `"transcript"` or `"gene"` level.}
#'   \item{`hasGenomicAnnotation`}{Logical value indicating whether
#'          detailed genomic annotation (promoter, exon, intron, etc.) was computed.}
#'   \item{`detailGenomicAnnotation`}{A data frame providing per-peak
#'          binary indicators for genomic categories.}
#'   \item{`annoStat`}{A data frame summarizing annotation category
#'          frequencies across the annotated peak set.}
#'   \item{`peakNum`}{Total number of annotated peaks.}
#' }
#' @name peakAnno
#' @docType data
#' @return csAnno object
#' @format A `csAnno` object containing 220 annotated peaks.
NULL


#' @title demo peak file
#'
#' @description  Peak in Grange object. See data-raw/example_data.R
#' @section Provenance:
#' The demo peaks were extracted from GSM6418464 in the GEO database
#' (\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6418464}).
#' @section Data structure:
#' A \code{GRanges} object with 220 genomic ranges and the following metadata columns:
#' \describe{
#'   \item{\code{seqnames}}{chr name}
#'   \item{\code{ranges}}{Peak ranges}
#'   \item{\code{strand}}{strand information}
#'   \item{\code{mcol}}{output from MACS2}
#' }
#' @name demo_peak
#' @docType data
#' @format A \code{GRanges} object with 200 rows and 5 metadata columns.
#' @return Grange object
NULL


#' @title Example data of a list of peak annotation
#'
#' @description A list of \code{csAnno} objects obtained by annotating multiple peak
#' files using \code{epiSeeker::annotateSeq()}.
#' See data-raw/example_data.R
#' @section Provenance:
#' The example peak annotation list was generated using several example peak
#' files returned by \code{getSampleFiles()}. Each peak file was annotated
#' using \code{epiSeeker::annotateSeq()} with \code{TxDb.Hsapiens.UCSC.hg38.knownGene},
#' \code{annoDb = "org.Hs.eg.db"} and \code{tssRegion = c(-3000, 3000)}.
#' @section Data structure:
#' A named list where each element is a \code{csAnno} S4 object produced by \code{annotateSeq()}.
#' @name peakAnnoList
#' @docType data
#' @format A a list of \code{csAnno} objects.
#' @return list of csAnno object
NULL


#' @title Example data of tagMatrix
#'
#' @description tagMatrix result used to demonstrate TSS enrichment visualization and
#' tag distribution plotting functions in \pkg{epiSeeker}. See data-raw/example_data.R
#' @section Provenance:
#' The tag matrix was generated using a sample peak file obtained from
#' \code{getSampleFiles()[[4]]}.
#' Peaks were imported via \code{readPeakFile()} and processed using
#' \code{epiSeeker::getTagMatrix()} with the following settings:
#' \itemize{
#'   \item Transcript database: \code{TxDb.Hsapiens.UCSC.hg38.knownGene}
#'   \item Annotation mode: \code{type = "start_site"}, \code{by = "gene"}
#'   \item TSS window: upstream 3000 bp, downstream 3000 bp
#'   \item Peak weight: column \code{"V5"} of the peak file
#'   \item Number of bins: \code{nbin = 500}
#' }
#' @section Data structure:
#' A numeric matrix in which:
#' \describe{
#'   \item{Rows}{Represent individual genes contributing tags around their TSS.}
#'   \item{Columns}{Represent evenly spaced bins across the TSS window from
#'         -3000 bp to +3000 bp (500 bins total).}
#' }
#'
#'
#' @format A numeric matrix with \code{n} genes × 500 bins.
#' @name tagMatrix
#' @docType data
#' @return matrix
NULL


#' @title motif reference for Homo sapiens
#'
#' @description A collection of transcription factor position weight matrices (PWMs)
#' retrieved from the JASPAR 2024 database.
#' This dataset is used to demonstrate motif enrichment, motif scanning,
#' and peak–motif association analyses in \pkg{epiSeeker}.
#' See data-raw/example_data.R
#' @section Provenance:
#' The PWM set was obtained using the JASPAR 2024 SQLite database bundled in
#' the \pkg{JASPAR2024} package.
#' Matrices were retrieved using \pkg{TFBSTools} with the following parameters:
#' \itemize{
#'   \item \code{collection = "CORE"}
#'   \item \code{all_versions = FALSE}
#'   \item \code{species = "Homo sapiens"}
#'   \item \code{tax_group = "vertebrates"}
#' }
#' @section Data structure:
#' A \code{TFBSTools::PWMatrixList} (or \code{PFMatrixList}) object containing
#' one PWM per transcription factor.
#' Each matrix stores nucleotide position weights across the TF binding motif,
#' with rows representing \code{A, C, G, T} and columns representing motif positions.
#' @format A \code{PFMatrixList} object containing PWMs for multiple human
#' transcription factors from the JASPAR 2024 CORE collection.
#' @name pwm_obj
#' @docType data
#' @return pwm_obj
NULL

#' @title demo base modification data
#'
#' @description A small example \code{bmData} object representing cytosine methylation
#' measurements from Bisulfite-Seq data.
#' This dataset is intended for demonstrating base-modification visualization,
#' regional methylation profiling, and \code{epiSeeker} workflows operating on
#' \code{bmData} objects.
#' See data-raw/example_data.R
#'
#' @section Provenance:
#' The example dataset was constructed from publicly available Bisulfite-Seq
#' data (GEO accession: \code{GSM6940395}, genome build: hg38).
#' The raw methylation coverage file (\code{*.bismark.cov.gz}) was imported
#' using \pkg{data.table::fread()}.
#'
#' A small genomic window on chromosome 22
#' (\code{[10525991, 10526342]}) was selected to create a lightweight example
#' dataset. The data were processed as follows:
#' \enumerate{
#'   \item Filter records where \code{chrom == 22} and positions fall within
#'         the chosen window.
#'   \item Convert chromosome name to UCSC style (\code{"chr22"}).
#'   \item Compute total coverage as: \code{Cov = methylated + unmethylated}.
#'   \item Extract columns: chromosome, position, coverage, and methylation
#'         percentage.
#'   \item Convert methylation percentage to a fraction.
#' }
#' @section Data structure:
#' A \code{bmData} S4 object containing one sample (\code{"acinar_methyl"}).
#' Each entry stores:
#' \describe{
#'   \item{\code{chr}}{Chromosome in UCSC format (e.g. \code{"chr22"}).}
#'   \item{\code{pos}}{Genomic coordinate of the cytosine.}
#'   \item{\code{Cov}}{Total read coverage at the site.}
#'   \item{\code{Methylation}}{Methylation level as a fraction (0–1).}
#' }
#' @name demo_bmdata
#' @return bmData object
#' @format A \code{bmData} object containing one sample.
#' @docType data
NULL

#' @title Result of seq2gene
#'
#' @description A character vector of gene IDs returned by \code{seq2gene()}, representing
#' genes associated with a subset of peaks.
#' This dataset is used to illustrate peak-to-gene mapping and regulatory
#' region annotation workflows in \pkg{epiSeeker}.
#' See data-raw/example_data.R
#' @name seq2gene_result
#' @section Data structure:
#' A character vector of gene identifiers (ENTREZ IDs) representing genes
#' linked to the example peak set via TSS proximity or flanking-gene search.
#'
#' @section Provenance:
#' The example peak set \code{demo_peak} was constructed by sampling up to
#' 10 peaks per autosome (chr1–chr22) from the ChIP-seq dataset
#' \code{GSM6418464}.
#' Peaks were imported using \code{readPeakFile()}, subset by chromosome,
#' and combined into a single \code{GRanges} object.
#'
#' The gene-level associations were then computed directly using:
#' \preformatted{
#' seq2gene_result <- seq2gene(
#'     demo_peak,
#'     tssRegion = c(-1000, 1000),
#'     flankDistance = 3000,
#'     txdb
#' )
#' }
#'
#' The resulting character vector of gene IDs was saved via
#' \code{data-raw/example_data.R}.
#' @format A character vector of gene IDs generated by
#' \code{seq2gene()} from the subset of peaks derived from demo_peak.
#' @docType data
#' @return vector of gene names
NULL

#' Name of the epiSeeker cache environment (internal static variable)
#' @format character vector
epiSeekerCache <- "epiSeekerEnv"
