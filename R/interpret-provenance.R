## Execution record of interpret(<csAnno>), in the format of clusterProfiler's interpret layer
## (copied from clusterProfiler R/interpret.R @ 2792236, Treywea/clusterProfiler
## feat/interpret-s4-generic), so that clusterProfiler::provenance() and its print method read
## an epiSeeker result unchanged. Keep the field names in step with clusterProfiler.

.prov_new <- function() {
    e <- new.env(parent = emptyenv())
    e$calls <- list()
    e$t_start <- Sys.time()
    e
}

## model id, provider and endpoint host of a "provider:model" string or an aisdk model object
.model_identity <- function(model) {
    out <- list(model_id = NA_character_, provider = NA_character_, endpoint = NA_character_)
    if (is.null(model)) return(out)
    if (is.character(model) && length(model) == 1L) {
        out$model_id <- model
        parts <- strsplit(model, ":", fixed = TRUE)[[1]]
        if (length(parts) > 1L) {
            out$provider <- parts[1]
            out$model_id <- paste(parts[-1], collapse = ":")
        }
        return(out)
    }
    grab <- function(...) {
        for (f in list(...)) {
            v <- tryCatch(f(), error = function(e) NULL)
            if (!is.null(v) && is.character(v) && length(v) == 1L && nzchar(v)) return(v)
        }
        NA_character_
    }
    out$model_id <- grab(function() model$modelId, function() model$model_id, function() model$model)
    out$provider <- grab(function() model$provider, function() model$providerName)
    url <- grab(function() model$config$baseURL, function() model$config$base_url)
    ## host only: a full URL can carry a key in its query string
    if (!is.na(url)) out$endpoint <- sub("^([a-z]+://[^/]+).*$", "\\1", url)
    out
}

## Record one model call. Never throws: the record must not be able to break a run.
.prov_record <- function(collector, stage, model, system, prompt, schema_name, temperature,
                         max_tokens, gen = NULL, t0 = NULL, error_message = NULL, reasoning = NULL) {
    if (is.null(collector)) return(invisible(NULL))
    tryCatch({
        id <- .model_identity(model)
        raw <- tryCatch(gen$raw_text, error = function(e) NULL)
        obj <- tryCatch(gen$object, error = function(e) NULL)
        fin <- tryCatch(gen$finish_reason, error = function(e) NULL)
        usage <- tryCatch(gen$usage, error = function(e) NULL)
        served <- tryCatch(gen$raw_response$model, error = function(e) NULL)
        if (!is.character(served) || length(served) != 1L) served <- NA_character_
        ## stricter than clusterProfiler's: an object without finish_reason "stop" is incomplete
        status <- if (!is.null(error_message)) "call_failed" else
            if (!is.null(obj) && identical(fin, "stop")) "ok" else
            if (!is.null(obj)) "incomplete" else "unparsed"
        num <- function(v) if (is.null(v) || !is.numeric(v) || length(v) != 1L) NA_integer_ else as.integer(v)
        collector$calls[[length(collector$calls) + 1L]] <- list(
            stage = stage,
            status = status,
            error_message = if (is.null(error_message)) NA_character_ else error_message,
            model_requested = id$model_id,
            model_served = served,
            model_matches_request = if (is.na(served)) NA else identical(served, id$model_id),
            provider = id$provider,
            endpoint = id$endpoint,
            schema_name = schema_name,
            temperature_requested = if (is.null(temperature)) NA_real_ else temperature,
            max_tokens_requested = max_tokens,
            temperature_applied = NA_real_,
            max_tokens_applied = NA_real_,
            reasoning_requested = if (is.null(reasoning)) NA_character_ else
                paste(utils::capture.output(utils::str(reasoning)), collapse = " "),
            finish_reason = if (is.null(fin)) NA_character_ else fin,
            response_chars = if (is.null(raw)) NA_integer_ else nchar(raw),
            prompt_tokens = num(usage$prompt_tokens),
            completion_tokens = num(usage$completion_tokens),
            total_tokens = num(usage$total_tokens),
            system_sha256 = .sha256(system),
            prompt_sha256 = .sha256(prompt),
            response_sha256 = .sha256(raw),
            runtime_sec = if (is.null(t0)) NA_real_ else
                round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
            timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        )
    }, error = function(e) invisible(NULL))
    invisible(NULL)
}

.pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)

.prov_finalize <- function(collector, entry_point, evidence_sha256 = NA_character_, extra = list()) {
    calls <- collector$calls
    st <- vapply(calls, function(cl) cl$status, character(1))
    prov <- c(list(
        entry_point = entry_point,
        ok = length(st) > 0 && all(st == "ok"),
        n_calls = length(calls),
        n_calls_failed = sum(st != "ok"),
        evidence_sha256 = evidence_sha256,
        started_utc = format(collector$t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        runtime_sec = round(as.numeric(difftime(Sys.time(), collector$t_start, units = "secs")), 2),
        software = list(
            r = paste0(R.version$major, ".", R.version$minor),
            epiSeeker = .pkg_ver("epiSeeker"),
            clusterProfiler = .pkg_ver("clusterProfiler"),
            aisdk = .pkg_ver("aisdk")
        ),
        calls = calls
    ), extra)
    class(prov) <- c("interpret_provenance", "list")
    prov
}
