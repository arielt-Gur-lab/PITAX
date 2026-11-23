# ============================================================
# Saved-project re-trimming without re-uploading original AB1 files.
# Reuses the EXACT trim decision as trim_one_ab1(). No dye/base recalling.
# ============================================================

pitax_saved_retrim_requirements <- function(result) {
  if (!is.list(result)) return("Missing processed read.")
  raw <- if (is.character(result$raw_seq) && length(result$raw_seq) == 1L) result$raw_seq else ""
  if (is.na(raw) || !nzchar(raw) || grepl("[^ACGTN]", raw)) {
    return("Stored original base calls (raw_seq) are missing or invalid.")
  }
  n <- nchar(raw)
  trace <- result$trace
  if (is.null(trace) || length(dim(trace)) != 2L || nrow(trace) < 3L || ncol(trace) < 4L) {
    return("Stored four-channel chromatogram is not available.")
  }
  pos <- suppressWarnings(as.integer(result$peak_pos))
  if (length(pos) != n || any(!is.finite(pos)) || any(pos < 1L | pos > nrow(trace))) {
    return("Stored peak positions do not match original base calls.")
  }
  cmap <- suppressWarnings(as.integer(result$channel_map[c("A", "C", "G", "T")]))
  if (length(cmap) != 4L || any(!is.finite(cmap)) ||
      any(cmap < 1L | cmap > ncol(trace)) || anyDuplicated(cmap)) {
    return("Original trace channel mapping is missing or invalid.")
  }
  if (!is.data.frame(result$summary) || nrow(result$summary) != 1L) {
    return("Stored trimming summary is missing.")
  }
  settings <- result$processing_settings
  required <- c("window", "min_peak_ratio", "min_len_before_collapse",
                "min_relative_signal", "bad_run_windows", "absolute_max_base_index",
                "min_usable_len")
  if (!is.list(settings) || !all(required %in% names(settings))) {
    return("Saved per-read trimming settings are missing; original AB1 re-upload is required.")
  }
  vals <- vapply(required, function(key) {
    value <- suppressWarnings(as.numeric(settings[[key]])[1])
    if (!length(value)) NA_real_ else value
  }, numeric(1))
  if (any(!is.finite(vals)) || any(vals <= 0)) {
    return("Stored trimming settings contain invalid parameters.")
  }
  if (vals[["window"]] > n) {
    return("Original read is shorter than the saved trimming window.")
  }
  NULL
}

pitax_saved_retrim_before <- function(result) {
  result <- ensure_curation_state(result)
  list(
    raw_seq = result$raw_seq,
    seq = result$seq,
    start = result$curation$trim_start,
    end = result$curation$trim_end,
    auto_start = result$curation$auto_trim_start,
    auto_end = result$curation$auto_trim_end,
    base_edits = result$curation$base_edits,
    revision = result$curation$revision,
    processing_settings = result$processing_settings
  )
}

pitax_retrim_saved_read <- function(result) {
  problem <- pitax_saved_retrim_requirements(result)
  if (!is.null(problem)) stop(problem, call. = FALSE)

  result <- ensure_curation_state(result)
  old <- pitax_saved_retrim_before(result)
  old_summary <- result$summary
  settings <- result$processing_settings
  raw <- as.character(result$raw_seq)
  metrics <- calculate_metrics(
    raw, as.matrix(result$trace),
    as.integer(result$peak_pos), result$channel_map
  )
  trim <- pitax_trim_decision(raw, metrics, settings)
  result$metrics <- metrics

  # Keep the original machine calls, AB1 evidence and read/assay assignment.
  # The newly computed bounds become a fresh AUTOMATIC trim baseline; manual
  # base edits are kept at immutable raw coordinates, but older manual trim
  # bounds must not silently override a user-requested re-trim.
  curation <- result$curation
  curation$auto_trim_start <- trim$start_index
  curation$auto_trim_end <- trim$end_index
  curation$auto_seq <- trim$trimmed_seq
  curation$auto_status <- trim$status
  curation$trim_start <- trim$start_index
  curation$trim_end <- trim$end_index
  result$curation <- curation

  # These are the same columns used when loading a fresh AB1 file. Other
  # metadata columns from a loaded project are kept intact.
  updates <- list(
    raw_length = nchar(raw),
    trimmed_length = trim$trimmed_length,
    trim_start = trim$start_index,
    trim_end = trim$end_index,
    collapse_index = trim$collapse_index,
    reason = trim$reason,
    mean_called_signal_trimmed = trim$mean_signal,
    median_peak_ratio_trimmed = trim$median_ratio,
    status = trim$status,
    error_message = ""
  )
  for (key in names(updates)) result$summary[[key]] <- updates[[key]]
  result <- curation_rebuild(result, settings)

  bounds_changed <- !identical(as.integer(old$start), as.integer(result$curation$trim_start)) ||
    !identical(as.integer(old$end), as.integer(result$curation$trim_end)) ||
    !identical(as.integer(old$auto_start), as.integer(result$curation$auto_trim_start)) ||
    !identical(as.integer(old$auto_end), as.integer(result$curation$auto_trim_end))
  sequence_changed <- !identical(as.character(old$seq), as.character(result$seq))
  if (bounds_changed || sequence_changed) {
    # Previous manual edits and audit trail survive, but old undo/redo
    # snapshots used the old auto-trim baseline and are no longer safe.
    history <- if (is.list(result$retrim_history)) result$retrim_history else list()
    history[[length(history) + 1L]] <- list(
      at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      previous_summary = old_summary,
      previous_sequence = old$seq,
      previous_trim = c(start = old$start, end = old$end),
      previous_auto_trim = c(start = old$auto_start, end = old$auto_end),
      previous_revision = old$revision
    )
    result$retrim_history <- history
    result$curation$undo_stack <- list()
    result$curation$redo_stack <- list()
    previous_revision <- suppressWarnings(as.integer(result$curation$revision)[1])
    if (!length(previous_revision) || !is.finite(previous_revision)) previous_revision <- 0L
    result$curation$revision <- previous_revision + 1L
    row <- data.frame(
      Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      Transaction_ID = curation_new_transaction_id(),
      Revision = as.integer(result$curation$revision),
      Action = "Re-trim saved chromatogram",
      Position = NA_integer_,
      Before = paste0(old$start, "-", old$end),
      After = paste0(result$curation$trim_start, "-", result$curation$trim_end),
      Method = "Stored AB1 trace; original base calls unchanged",
      Evidence = paste0("Algorithm: shared pitax_trim_decision; reason: ", trim$reason),
      Details = "Existing manual base edits kept; manual trim bounds reset; previous snapshot in retrim_history",
      stringsAsFactors = FALSE
    )
    result$curation$audit_log <- rbind(result$curation$audit_log, row)
    result <- curation_rebuild(result, settings)
  }
  result
}

# Atomic in-memory preview. Any invalid legacy read stops the batch BEFORE
# the server changes results, names, assignments or downstream evidence.
pitax_retrim_saved_batch <- function(results, selected_ids = names(results)) {
  if (!is.list(results) || !length(results)) stop("No saved trimmed reads are available.", call. = FALSE)
  selected_ids <- unique(as.character(selected_ids))
  if (!length(selected_ids) || any(!selected_ids %in% names(results))) {
    stop("Select existing saved reads to re-trim.", call. = FALSE)
  }
  out <- list()
  rows <- list()
  before <- list()
  errors <- character()
  for (id in selected_ids) {
    old <- results[[id]]
    err <- pitax_saved_retrim_requirements(old)
    if (!is.null(err)) {
      errors <- c(errors, paste0(id, ": ", err))
      next
    }
    new <- tryCatch(pitax_retrim_saved_read(old), error = function(e) e)
    if (inherits(new, "error")) {
      errors <- c(errors, paste0(id, ": ", conditionMessage(new)))
      next
    }
    before[[id]] <- pitax_saved_retrim_before(old)
    out[[id]] <- new
    summary_start <- function(x) {
      pos <- suppressWarnings(as.integer(x$summary$trim_start[1]))
      if (length(pos) && is.finite(pos)) as.character(pos) else "-"
    }
    summary_end <- function(x) {
      pos <- suppressWarnings(as.integer(x$summary$trim_end[1]))
      if (length(pos) && is.finite(pos)) as.character(pos) else "-"
    }
    rows[[id]] <- data.frame(
      Read = id,
      Before = paste(summary_start(old), summary_end(old), sep = "-"),
      After = paste(summary_start(new), summary_end(new), sep = "-"),
      Previous_length = nchar(as.character(old$seq)),
      New_length = nchar(as.character(new$seq)),
      Sequence_changed = !identical(as.character(old$seq), as.character(new$seq)),
      stringsAsFactors = FALSE
    )
  }
  if (length(errors)) {
    stop(paste(c("Re-trim aborted; no results changed:", errors), collapse = "\n"), call. = FALSE)
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  list(
    results = out,
    original = before,
    summary = summary,
    changed_sources = as.character(summary$Read[summary$Sequence_changed])
  )
}
