#' Interpret annotated regions with a language model
#'
#' Methods of [clusterProfiler::interpret()] for epiSeeker objects. Each builds an evidence
#' digest without calling a model, then makes one structured call. The digest is what the
#' model sees: every section says how many items it sends and how many there are, and is
#' hashed into the execution record, which [clusterProfiler::provenance()] returns.
#'
#' * `csAnno` (from [annotateSeq()]): what this set of regions looks like, whether it fits
#'   the declared assay, which genes and programs it points to, and, when `question` is given,
#'   a direct answer with ranked candidates. The digest lists, in this order, the top regions
#'   with their nearest genes, promoter genes, overlap with reference samples, genomic
#'   distribution, distance to TSS and region width.
#' * `GRanges`: annotated with [annotateSeq()] first, then interpreted as a `csAnno`. `TxDb`
#'   must be given, so that regions are never annotated against a default genome.
#' * `csAnnoList` (experimental, not yet evaluated against a baseline): compares the groups --
#'   their distributions, widths, distances to TSS, shared and group-specific promoter genes,
#'   and top regions.
#' * `peakOverlapResult`: not interpreted alone; pass it as `overlap` of the `csAnno` method.
#'
#' A call counts as successful only when the model finished (`finish_reason == "stop"`) and
#' returned a non-empty object; otherwise the fields are `NA`, a warning is raised and
#' `provenance(res)$ok` is `FALSE`.
#'
#' @param x A `csAnno`, `GRanges` or `csAnnoList` object.
#' @param context Experimental background, e.g. "ChIP-seq in mouse liver, HFD vs chow".
#' @param assay Declared assay or target, e.g. "TF", "H3K4me3", "H3K27me3", "ATAC".
#' @param question A question to answer directly, e.g. "Which transcription factor was
#'   profiled?". Without one, `answer` is "None".
#' @param overlap Optional result of [enrichPeakOverlap()] for the same regions. Samples with
#'   `p.adjust < 0.05` are sent, ranked by the first of `jaccard`, `frac_query`,
#'   `fold_enrichment` or `N_OL` the table has; a `label` column, if present, names them,
#'   otherwise GEO samples are named from `gsminfo`.
#' @param n_regions Number of top regions sent (per group for a `csAnnoList`). Default 200
#'   for a `csAnno`, 30 per group for a `csAnnoList`.
#' @param signal Column of the regions used to rank them. `NULL` takes the first of
#'   `signalValue`, `signal`, `score`, `fold_enrichment`, `V7`, `V5`; without any, regions are
#'   ranked by distance to the TSS.
#' @param model Model, as `"provider:model"` or an aisdk model object. `NULL` uses the default
#'   set with `aisdk::set_model()`.
#' @param max_tokens Maximum tokens of the response. Default 4096.
#' @param temperature Sampling temperature. `NULL` (default) does not send one.
#' @param reasoning Optional reasoning configuration passed to the provider; not sent when `NULL`.
#' @param TxDb,annoDb,tssRegion,level For a `GRanges`: passed to [annotateSeq()]. `TxDb` is
#'   required.
#' @param ... For a `GRanges`, further arguments of the `csAnno` method; otherwise not used.
#' @return A `seq_interpretation` list. For a `csAnno`: `answer`, `candidates`,
#'   `binding_profile`, `assay_consistency`, `target_programs`, `cofactors`, `narrative`. For
#'   a `csAnnoList`: `answer`, `candidates`, `group_profiles`, `differences`, `shared`,
#'   `target_programs`, `narrative`. The evidence digest is in `attr(, "evidence")` and the
#'   execution record in `attr(, "provenance")`.
#' @name interpret-epiSeeker
#' @aliases interpret,csAnno-method interpret,GRanges-method interpret,csAnnoList-method interpret,peakOverlapResult-method
#' @importFrom clusterProfiler interpret
#' @exportMethod interpret
#' @examples
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' peak <- readPeakFile(system.file("extdata", "sample_peaks.txt", package = "epiSeeker"))
#' anno <- annotateSeq(peak, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
#' ## the call below needs a model; set one up with aisdk (an API key) to run it
#' if (requireNamespace("aisdk", quietly = TRUE) && nzchar(Sys.getenv("OPENAI_API_KEY"))) {
#'     res <- interpret(anno, context = "ChIP-seq in a human cell line", assay = "TF",
#'                      question = "Which transcription factor was profiled?",
#'                      model = "openai:gpt-4o-mini")
#'     res$answer
#'     clusterProfiler::provenance(res)
#'     ## straight from the peaks
#'     res <- interpret(peak, TxDb = txdb, annoDb = "org.Hs.eg.db", assay = "TF",
#'                      model = "openai:gpt-4o-mini")
#' }
NULL

#' @rdname interpret-epiSeeker
setMethod("interpret", "csAnno", function(x, context = NULL, assay = NULL, question = NULL,
                                          overlap = NULL, n_regions = 200, signal = NULL,
                                          model = NULL, max_tokens = 4096, temperature = NULL,
                                          reasoning = NULL, ...) {
    ev <- .seq_digest(x, context = context, assay = assay, question = question,
                      overlap = overlap, n_regions = n_regions, signal = signal)
    .seq_interpret_run(ev, schema = .seq_schema(), entry_point = "interpret,csAnno",
                       model = model, max_tokens = max_tokens, temperature = temperature,
                       reasoning = reasoning)
})

#' @rdname interpret-epiSeeker
setMethod("interpret", "GRanges", function(x, TxDb = NULL, annoDb = NULL,
                                           tssRegion = c(-3000, 3000), level = "transcript", ...) {
    if (is.null(TxDb)) {
        stop("Give `TxDb` for the genome of these regions; interpret() does not annotate ",
             "against a default genome. Or annotate first with annotateSeq() and pass the csAnno.")
    }
    db <- if (is.null(annoDb)) "NULL" else if (is.character(annoDb)) annoDb else class(annoDb)[1]
    message("Annotating ", length(x), " regions with annotateSeq(tssRegion = c(", tssRegion[1],
            ", ", tssRegion[2], "), level = \"", level, "\", annoDb = ", db, ") before interpreting.")
    anno <- annotateSeq(x, TxDb = TxDb, annoDb = annoDb, tssRegion = tssRegion, level = level,
                        verbose = FALSE)
    interpret(anno, ...)
})

#' @rdname interpret-epiSeeker
setMethod("interpret", "csAnnoList", function(x, context = NULL, assay = NULL, question = NULL,
                                              n_regions = 30, signal = NULL, model = NULL,
                                              max_tokens = 4096, temperature = NULL,
                                              reasoning = NULL, ...) {
    ev <- .seq_digest_list(x, context = context, assay = assay, question = question,
                           n_regions = n_regions, signal = signal)
    .seq_interpret_run(ev, schema = .seq_schema_list(), entry_point = "interpret,csAnnoList",
                       model = model, max_tokens = max_tokens, temperature = temperature,
                       reasoning = reasoning)
})

#' @rdname interpret-epiSeeker
setMethod("interpret", "peakOverlapResult", function(x, ...) {
    stop("An overlap table is evidence about a set of regions, not something to interpret ",
         "alone. Pass it with the regions: interpret(<csAnno>, overlap = <this table>).")
})

## ---- shared skeleton: one model call, success check and execution record ---------------

.seq_interpret_run <- function(ev, schema, entry_point, model, max_tokens, temperature, reasoning) {
    if (!requireNamespace("aisdk", quietly = TRUE)) {
        stop("interpret() needs the 'aisdk' package: install.packages('aisdk')")
    }
    prov <- .prov_new()
    gen <- tryCatch(
        .seq_generate(prov, model = model, prompt = ev$text, schema = schema,
                      system = .seq_system_prompt(), temperature = temperature,
                      max_tokens = max_tokens, reasoning = reasoning),
        error = function(e) NULL
    )
    fields <- names(schema$properties)
    obj <- if (.seq_call_ok(gen)) gen$object else NULL
    if (is.null(obj)) {
        warning("The model returned no complete answer; see clusterProfiler::provenance() of the result.",
                call. = FALSE)
        obj <- lapply(schema$properties, .na_of)
    }
    res <- structure(obj[fields], names = fields, class = c("seq_interpretation", "list"))
    attr(res, "evidence") <- ev
    attr(res, "provenance") <- .prov_finalize(
        prov, entry_point = entry_point,
        evidence_sha256 = .sha256(ev$sha256),
        extra = list(evidence_sections_sha256 = ev$sha256, totals = ev$totals)
    )
    res
}

## the empty value of a schema field, used when the call failed
.na_of <- function(s) {
    if (identical(s$type, "array")) return(list())
    if (identical(s$type, "object")) return(lapply(s$properties, .na_of))
    NA_character_
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

.seq_system_prompt <- function() {
    paste(
        "You interpret genomic regions (peaks from ChIP-seq, ATAC-seq, CUT&Tag or similar)",
        "from a deterministic summary of them. Use only the evidence given:",
        "- Name only genes, regions, groups and reference samples that appear in the input.",
        "- Quote percentages, widths and counts as given; do not compute new ones.",
        "- Rank reference samples by the order given; a lower-ranked sample is weaker evidence.",
        "- Promoter-wide reference sets overlap most promoter-rich queries; do not take such an",
        "  overlap alone as the identity of the factor.",
        "- If the input reports a GENOME BUILD SUSPECT or dropped peaks, say so first and treat",
        "  statements that rest on the annotation (genes, distribution) as unreliable.",
        "- Say 'Not assessable' when the evidence cannot decide, rather than guessing.",
        sep = "\n"
    )
}

## ---- schemas: at most 8 fields, all required; the direct answer first, the narrative last --

.candidates_field <- function() {
    aisdk::z_array(aisdk::z_object(
        label = aisdk::z_string("A candidate answer, e.g. a factor, mark, cell type or group."),
        support = aisdk::z_string(paste(
            "The evidence for it, cited as it appears in the input: reference sample ids,",
            "gene symbols, group names or distribution categories."))
    ), description = "Ranked candidates, best first; empty when no question was asked.")
}

.answer_field <- function() {
    aisdk::z_string("Direct answer to the question in a few words, or 'None' when no question was asked.")
}

.seq_schema <- function() {
    aisdk::z_object(
        answer = .answer_field(),
        candidates = .candidates_field(),
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

.seq_schema_list <- function() {
    aisdk::z_object(
        answer = .answer_field(),
        candidates = .candidates_field(),
        group_profiles = aisdk::z_array(aisdk::z_object(
            group = aisdk::z_string("Group name as given."),
            profile = aisdk::z_string("What its regions look like, citing the numbers given.")
        ), description = "One entry per group, in the order given."),
        differences = aisdk::z_string(
            "The main differences between groups, citing numbers and group-specific genes."),
        shared = aisdk::z_string("What the groups have in common, citing the shared genes."),
        target_programs = aisdk::z_string(
            "Biological programs the group-specific genes point to, naming genes and groups."),
        narrative = aisdk::z_string("A short paragraph that ties the comparison together.")
    )
}

#' @method print seq_interpretation
#' @export
print.seq_interpretation <- function(x, ...) {
    cat("<epiSeeker interpretation>\n")
    cat("answer            :", format(x$answer), "\n")
    if (length(x$candidates)) {
        for (i in seq_along(x$candidates)) {
            cat(sprintf("  %d. %s -- %s\n", i, x$candidates[[i]]$label, x$candidates[[i]]$support))
        }
    }
    if (!is.null(x$assay_consistency)) {
        cat("assay consistency :", format(x$assay_consistency$verdict),
            "-", format(x$assay_consistency$reason), "\n")
    }
    if (length(x$group_profiles)) {
        for (g in x$group_profiles) cat(sprintf("  [%s] %s\n", g$group, g$profile))
    }
    for (f in intersect(c("binding_profile", "differences", "shared", "target_programs",
                          "cofactors", "narrative"), names(x))) {
        cat(sprintf("%-18s: %s\n", f, format(x[[f]])))
    }
    p <- attr(x, "provenance")
    if (!is.null(p)) cat("provenance        :", if (isTRUE(p$ok)) "ok" else "INCOMPLETE", "\n")
    invisible(x)
}
