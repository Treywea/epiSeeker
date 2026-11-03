#' Interpret annotated regions with a language model
#'
#' The `csAnno` method of [clusterProfiler::interpret()]. It reads an annotation from
#' [annotateSeq()] (and, optionally, an [enrichPeakOverlap()] result), builds an evidence
#' digest without calling a model, and asks the model one structured question: what this set
#' of regions looks like, whether it fits the declared assay, which genes and programs it
#' points to, and, when `question` is given, a direct answer with ranked candidates.
#'
#' The digest is what the model sees, in this order: top regions with their nearest genes,
#' promoter genes, overlap with reference samples, genomic distribution, distance to TSS and
#' region width. Every section says how many items it sends and how many there are, and is
#' hashed into the execution record, which [clusterProfiler::provenance()] returns.
#'
#' A call counts as successful only when the model finished (`finish_reason == "stop"`) and
#' returned a non-empty object; otherwise the fields are `NA`, a warning is raised and
#' `provenance(res)$ok` is `FALSE`.
#'
#' @param x A `csAnno` object from [annotateSeq()].
#' @param context Experimental background, e.g. "ChIP-seq in mouse liver, HFD vs chow".
#' @param assay Declared assay or target, e.g. "TF", "H3K4me3", "H3K27me3", "ATAC".
#' @param question A question to answer directly, e.g. "Which transcription factor was
#'   profiled?". Without one, `answer` is "None".
#' @param overlap Optional result of [enrichPeakOverlap()] for the same regions. Samples with
#'   `p.adjust < 0.05` are sent, ranked by the first of `jaccard`, `frac_query`,
#'   `fold_enrichment` or `N_OL` the table has; a `label` column, if present, names them,
#'   otherwise GEO samples are named from `gsminfo`.
#' @param n_regions Number of top regions sent. Default 200.
#' @param signal Column of the regions used to rank them. `NULL` takes the first of
#'   `signalValue`, `signal`, `score`, `fold_enrichment`, `V7`, `V5`; without any, regions are
#'   ranked by distance to the TSS.
#' @param model Model, as `"provider:model"` or an aisdk model object. `NULL` uses the default
#'   set with `aisdk::set_model()`.
#' @param max_tokens Maximum tokens of the response. Default 4096.
#' @param temperature Sampling temperature. `NULL` (default) does not send one.
#' @param reasoning Optional reasoning configuration passed to the provider; not sent when `NULL`.
#' @param ... Not used.
#' @return A `seq_interpretation` list with `answer`, `candidates`, `binding_profile`,
#'   `assay_consistency`, `target_programs`, `cofactors` and `narrative`. The evidence digest
#'   is in `attr(, "evidence")` and the execution record in `attr(, "provenance")`.
#' @name interpret-csAnno
#' @aliases interpret,csAnno-method
#' @importFrom clusterProfiler interpret
#' @exportMethod interpret
#' @examples
#' \dontrun{
#' peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
#' anno <- annotateSeq(peak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
#'                     annoDb = "org.Hs.eg.db")
#' res <- interpret(anno, context = "ChIP-seq in a human cell line", assay = "TF",
#'                  question = "Which transcription factor was profiled?",
#'                  model = "openai:gpt-4o-mini")
#' res$answer
#' clusterProfiler::provenance(res)
#' }
setMethod("interpret", "csAnno", function(x, context = NULL, assay = NULL, question = NULL,
                                          overlap = NULL, n_regions = 200, signal = NULL,
                                          model = NULL, max_tokens = 4096, temperature = NULL,
                                          reasoning = NULL, ...) {
    if (!requireNamespace("aisdk", quietly = TRUE)) {
        stop("interpret() needs the 'aisdk' package: install.packages('aisdk')")
    }
    ev <- .seq_digest(x, context = context, assay = assay, question = question,
                      overlap = overlap, n_regions = n_regions, signal = signal)
    prov <- .prov_new()
    gen <- tryCatch(
        .seq_generate(prov, model = model, prompt = ev$text, schema = .seq_schema(),
                      system = .seq_system_prompt(), temperature = temperature,
                      max_tokens = max_tokens, reasoning = reasoning),
        error = function(e) NULL
    )
    obj <- if (.seq_call_ok(gen)) gen$object else NULL
    if (is.null(obj)) {
        warning("The model returned no complete answer; see clusterProfiler::provenance() of the result.",
                call. = FALSE)
        obj <- stats::setNames(as.list(rep(NA_character_, length(.seq_fields))), .seq_fields)
        obj$candidates <- list()
        obj$assay_consistency <- list(verdict = NA_character_, reason = NA_character_)
    }
    res <- structure(obj[.seq_fields], class = c("seq_interpretation", "list"))
    attr(res, "evidence") <- ev
    attr(res, "provenance") <- .prov_finalize(
        prov, entry_point = "interpret,csAnno",
        evidence_sha256 = .sha256(ev$sha256),
        extra = list(evidence_sections_sha256 = ev$sha256, totals = ev$totals)
    )
    res
})

.seq_fields <- c("answer", "candidates", "binding_profile", "assay_consistency",
                 "target_programs", "cofactors", "narrative")

## Seven fields, all required: a relay that drops long outputs loses optional fields first.
## The direct answer and the ranked candidates come first, the narrative last.
.seq_schema <- function() {
    aisdk::z_object(
        answer = aisdk::z_string(paste(
            "Direct answer to the question in a few words, or 'None' when no question was asked.")),
        candidates = aisdk::z_array(aisdk::z_object(
            label = aisdk::z_string("A candidate answer, e.g. a factor, mark or cell type."),
            support = aisdk::z_string(paste(
                "The evidence for it, cited as it appears in the input: reference sample ids,",
                "gene symbols or distribution categories."))
        ), description = "Ranked candidates, best first; empty when no question was asked."),
        binding_profile = aisdk::z_string(paste(
            "What the regions look like (promoter-proximal or distal, narrow or broad),",
            "citing the percentages and widths given.")),
        assay_consistency = aisdk::z_object(
            verdict = aisdk::z_enum(c("Consistent", "Inconsistent", "Not assessable")),
            reason = aisdk::z_string("Why, citing the numbers given.")
        ),
        target_programs = aisdk::z_string(
            "Biological programs the genes point to, naming the genes from the input."),
        cofactors = aisdk::z_string(paste(
            "Co-binding factors or cell types suggested by the overlap with reference samples,",
            "citing the sample ids; 'None' when no overlap was given.")),
        narrative = aisdk::z_string("A short paragraph that ties the evidence together.")
    )
}

.seq_system_prompt <- function() {
    paste(
        "You interpret a set of genomic regions (peaks from ChIP-seq, ATAC-seq, CUT&Tag or similar)",
        "from a deterministic summary of them. Use only the evidence given:",
        "- Name only genes, regions and reference samples that appear in the input.",
        "- Quote percentages, widths and counts as given; do not compute new ones.",
        "- Rank reference samples by the order given; a lower-ranked sample is weaker evidence.",
        "- Promoter-wide reference sets overlap most promoter-rich queries; do not take such an",
        "  overlap alone as the identity of the factor.",
        "- Say 'Not assessable' when the evidence cannot decide, rather than guessing.",
        sep = "\n"
    )
}

## One structured call, recorded whether it succeeds or not.
.seq_generate <- function(collector, model, prompt, schema, system, temperature, max_tokens,
                          reasoning = NULL, stage = "interpret") {
    t0 <- Sys.time()
    args <- list(model = model, prompt = prompt, schema = schema, schema_name = "seq_interpretation",
                 system = system, max_tokens = max_tokens)
    ## not sent when NULL: gateways have been seen to ignore temperature, and an unsent
    ## parameter cannot be mistaken for an applied one
    if (!is.null(temperature)) args$temperature <- temperature
    if (!is.null(reasoning)) args$reasoning <- reasoning
    cap <- tryCatch(list(gen = do.call(aisdk::generate_object, args), err = NULL),
                    error = function(e) list(gen = NULL, err = conditionMessage(e)))
    .prov_record(collector, stage, model, system, prompt, "seq_interpretation", temperature,
                 max_tokens, gen = cap$gen, t0 = t0, error_message = cap$err, reasoning = reasoning)
    cap$gen
}

## success = the model finished and returned something; aisdk's "blocked" marks a transport
## failure, not an answer
.seq_call_ok <- function(gen) {
    !is.null(gen) && identical(gen$finish_reason, "stop") && length(gen$object) > 0
}

#' @method print seq_interpretation
#' @export
print.seq_interpretation <- function(x, ...) {
    cat("<epiSeeker interpretation>\n")
    cat("answer            :", format(x$answer), "\n")
    if (length(x$candidates) && !all(is.na(x$candidates))) {
        for (i in seq_along(x$candidates)) {
            cat(sprintf("  %d. %s -- %s\n", i, x$candidates[[i]]$label, x$candidates[[i]]$support))
        }
    }
    cat("assay consistency :", format(x$assay_consistency$verdict),
        "-", format(x$assay_consistency$reason), "\n")
    for (f in c("binding_profile", "target_programs", "cofactors", "narrative")) {
        cat(sprintf("%-18s: %s\n", f, format(x[[f]])))
    }
    p <- attr(x, "provenance")
    if (!is.null(p)) cat("provenance        :", if (isTRUE(p$ok)) "ok" else "INCOMPLETE", "\n")
    invisible(x)
}
