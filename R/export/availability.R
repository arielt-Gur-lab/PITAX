# ============================================================
# Export availability and extractors
# ============================================================

pitax_export_status <- function(status, reason = "", ready_ids = character(), blocked_ids = character(),
                                record_count = 0L) {
  list(
    status = as.character(status)[1],
    reason = as.character(reason)[1],
    ready_ids = unique(as.character(ready_ids)),
    blocked_ids = unique(as.character(blocked_ids)),
    record_count = as.integer(record_count)[1]
  )
}

pitax_export_assignment_ready <- function(snap) {
  asg <- snap$read_assignments
  if (!is.data.frame(asg) || !nrow(asg)) return(FALSE)
  err <- tryCatch(stage2_identity_error(asg), error = function(e) conditionMessage(e))
  is.null(err)
}

pitax_export_active_consensus_revision <- function(snap, analysis_id) {
  cs <- snap$consensus_set
  if (!is.list(cs) || !is.list(cs$records) || !analysis_id %in% names(cs$records)) return(NA_integer_)
  rec <- cs$records[[analysis_id]]
  if (is.list(rec$curation) && !is.null(rec$curation$revision)) {
    return(suppressWarnings(as.integer(rec$curation$revision)[1]))
  }
  NA_integer_
}

pitax_export_blast_job_is_active <- function(snap, job_row) {
  analysis_id <- as.character(job_row$original_name)[1]
  job_rev <- suppressWarnings(as.integer(job_row$consensus_revision)[1])
  active_rev <- pitax_export_active_consensus_revision(snap, analysis_id)
  if (!is.na(job_rev) && !is.na(active_rev) && !identical(job_rev, active_rev)) return(FALSE)
  status <- toupper(as.character(job_row$status)[1])
  !identical(status, "STALE")
}

pitax_export_filter_scope_ids <- function(snap, selection) {
  asg <- if (is.data.frame(snap$read_assignments)) snap$read_assignments else data.frame()
  source_ids <- if (nrow(asg) && "Source_ID" %in% names(asg)) as.character(asg$Source_ID) else names(snap$results)
  if (!length(source_ids)) source_ids <- character()

  sample_mode <- if (is.null(selection$sample_mode)) "all" else as.character(selection$sample_mode)[1]
  selected_samples <- unique(as.character(selection$selected_samples %||% character()))
  if (identical(sample_mode, "current") && nzchar(as.character(snap$current_sample %||% ""))) {
    selected_samples <- as.character(snap$current_sample)[1]
  }
  if (identical(sample_mode, "all")) {
    keep <- source_ids
  } else {
    keep <- intersect(source_ids, selected_samples)
    if (!length(keep) && length(selected_samples)) keep <- selected_samples
  }

  isolates <- unique(as.character(selection$isolates %||% character()))
  loci <- unique(as.character(selection$loci %||% character()))
  assays <- unique(as.character(selection$assays %||% character()))
  if (nrow(asg) && length(keep)) {
    sub <- asg[asg$Source_ID %in% keep, , drop = FALSE]
    if (length(isolates)) sub <- sub[as.character(sub$Isolate) %in% isolates, , drop = FALSE]
    if (length(loci)) sub <- sub[as.character(sub$Locus) %in% loci, , drop = FALSE]
    if (length(assays) && "Assay_ID" %in% names(sub)) {
      sub <- sub[as.character(sub$Assay_ID) %in% assays, , drop = FALSE]
    }
    keep <- as.character(sub$Source_ID)
  }
  unique(keep)
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

pitax_export_final_name <- function(snap, source_id) {
  if (is.data.frame(snap$rename) && nrow(snap$rename) && "Original_name" %in% names(snap$rename)) {
    idx <- match(source_id, snap$rename$Original_name)
    if (!is.na(idx)) return(as.character(snap$rename$New_name[idx]))
  }
  if (is.data.frame(snap$read_assignments) && nrow(snap$read_assignments)) {
    idx <- match(source_id, snap$read_assignments$Source_ID)
    if (!is.na(idx) && "Final_Name" %in% names(snap$read_assignments)) {
      nm <- as.character(snap$read_assignments$Final_Name[idx])
      if (nzchar(nm)) return(nm)
    }
  }
  as.character(source_id)
}

pitax_export_settings_for_result <- function(result, fallback = list()) {
  if (is.list(result) && is.list(result$processing_settings)) result$processing_settings else fallback
}

# ---- Extractors -------------------------------------------------------------

pitax_export_extract_upload_file_inventory <- function(snap, scope_ids = NULL) {
  uf <- snap$upload_files
  if (!is.data.frame(uf) || !nrow(uf)) return(data.frame())
  out <- data.frame(
    Source_ID = as.character(uf$Source_ID),
    File = as.character(uf$name),
    Size_KB = round(as.numeric(uf$size) / 1024, 1),
    Bytes_available = as.logical(uf$Bytes_available),
    stringsAsFactors = FALSE
  )
  if (!is.null(scope_ids) && length(scope_ids)) out <- out[out$Source_ID %in% scope_ids, , drop = FALSE]
  out
}

pitax_export_extract_upload_raw_ab1 <- function(snap, scope_ids = NULL) {
  uf <- snap$upload_files
  if (!is.data.frame(uf) || !nrow(uf)) return(data.frame())
  out <- uf[as.logical(uf$Bytes_available) & file.exists(as.character(uf$datapath)), , drop = FALSE]
  if (!is.null(scope_ids) && length(scope_ids)) out <- out[out$Source_ID %in% scope_ids, , drop = FALSE]
  out
}

pitax_export_extract_assay_profiles_draft <- function(snap, scope_ids = NULL) {
  df <- snap$assay_profiles
  if (!is.data.frame(df)) return(data.frame())
  df
}

pitax_export_extract_project_defaults_draft <- function(snap, scope_ids = NULL) {
  pitax_export_named_list_to_df(snap$project_defaults)
}

pitax_export_extract_processing_settings_applied <- function(snap, scope_ids = NULL) {
  rows <- list()
  for (sid in names(snap$results)) {
    if (!is.null(scope_ids) && length(scope_ids) && !sid %in% scope_ids) next
    ps <- snap$results[[sid]]$processing_settings
    if (!is.list(ps) || !length(ps)) next
    for (nm in names(ps)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Source_ID = sid, Setting = nm,
        Value = paste(as.character(unlist(ps[[nm]])), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pitax_export_extract_assign_identity_table <- function(snap, scope_ids = NULL) {
  df <- snap$read_assignments
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids)) df <- df[df$Source_ID %in% scope_ids, , drop = FALSE]
  df
}

pitax_export_extract_assign_rename_map <- function(snap, scope_ids = NULL) {
  df <- snap$rename
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "Original_name" %in% names(df)) {
    df <- df[df$Original_name %in% scope_ids, , drop = FALSE]
  }
  df
}

pitax_export_extract_assign_architecture <- function(snap, scope_ids = NULL) {
  arch <- snap$architecture
  if (!is.list(arch)) return(data.frame())
  rows <- list()
  for (tbl in c("assays", "isolates", "loci", "reads")) {
    df <- arch[[tbl]]
    if (!is.data.frame(df) || !nrow(df)) next
    key_col <- names(df)[1]
    for (i in seq_len(nrow(df))) {
      rows[[length(rows) + 1L]] <- data.frame(
        Table = tbl,
        Key = as.character(df[[key_col]][i]),
        Payload = paste(vapply(names(df), function(nm) paste0(nm, "=", as.character(df[[nm]][i])), character(1)), collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pitax_export_extract_qc_summary <- function(snap, scope_ids = NULL) {
  df <- snap$summary
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  df$final_name <- vapply(as.character(df$sample_id), function(sid) pitax_export_final_name(snap, sid), character(1))
  if (!is.null(scope_ids) && length(scope_ids)) df <- df[df$sample_id %in% scope_ids, , drop = FALSE]
  df
}

pitax_export_extract_qc_ab1_run_evidence <- function(snap, scope_ids = NULL) {
  if (!length(snap$results)) return(data.frame())
  df <- tryCatch(ab1_evidence_run_summary(snap$results), error = function(e) data.frame())
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if ("Sample" %in% names(df)) {
    df$Final_Name <- vapply(as.character(df$Sample), function(sid) pitax_export_final_name(snap, sid), character(1))
    if (!is.null(scope_ids) && length(scope_ids)) df <- df[df$Sample %in% scope_ids, , drop = FALSE]
  }
  df
}

pitax_export_extract_qc_ab1_base_evidence <- function(snap, scope_ids = NULL) {
  rows <- list()
  for (sid in names(snap$results)) {
    if (!is.null(scope_ids) && length(scope_ids) && !sid %in% scope_ids) next
    result <- snap$results[[sid]]
    if (!is.list(result) || !is.list(result$ab1_evidence)) next
    detail <- result$ab1_evidence$detail
    if (!is.data.frame(detail) || !nrow(detail)) next
    detail$Sample_ID <- sid
    detail$Final_Name <- pitax_export_final_name(snap, sid)
    rows[[length(rows) + 1L]] <- detail
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pitax_export_extract_qc_peak_flags <- function(snap, scope_ids = NULL) {
  if (!length(snap$results)) return(data.frame())
  df <- tryCatch(
    collect_ambiguous_peak_flags(snap$results, scope = "trimmed", settings = snap$settings),
    error = function(e) data.frame()
  )
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "Sample" %in% names(df)) {
    df <- df[df$Sample %in% scope_ids, , drop = FALSE]
  }
  df
}

pitax_export_extract_qc_active_curation <- function(snap, scope_ids = NULL) {
  df <- tryCatch(collect_curation_log(snap$results), error = function(e) data.frame())
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  # Active = latest revision actions that are not undo/redo bookkeeping alone.
  if ("Action" %in% names(df)) {
    df <- df[!grepl("^Undo|^Redo", as.character(df$Action), ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(scope_ids) && length(scope_ids) && "Sample" %in% names(df)) {
    df <- df[df$Sample %in% scope_ids, , drop = FALSE]
  }
  if (nrow(df)) df$Export_class <- "active_curation"
  df
}

pitax_export_extract_qc_historical_audit <- function(snap, scope_ids = NULL) {
  df <- tryCatch(collect_curation_log(snap$results), error = function(e) data.frame())
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "Sample" %in% names(df)) {
    df <- df[df$Sample %in% scope_ids, , drop = FALSE]
  }
  if (nrow(df)) df$Export_class <- "historical_audit"
  df
}

pitax_export_extract_qc_trimmed_sequences <- function(snap, scope_ids = NULL) {
  out <- list()
  for (sid in names(snap$results)) {
    if (!is.null(scope_ids) && length(scope_ids) && !sid %in% scope_ids) next
    result <- snap$results[[sid]]
    seq <- if (is.null(result$seq)) "" else as.character(result$seq)[1]
    asg <- if (is.data.frame(snap$read_assignments) && nrow(snap$read_assignments)) {
      snap$read_assignments[match(sid, snap$read_assignments$Source_ID), , drop = FALSE]
    } else NULL
    out[[sid]] <- list(
      id = sid,
      final_name = pitax_export_final_name(snap, sid),
      sequence = seq,
      sequence_role = "trimmed_curated_read",
      direction = if (!is.null(asg) && nrow(asg) && "Direction" %in% names(asg)) as.character(asg$Direction[1]) else "",
      isolate = if (!is.null(asg) && nrow(asg) && "Isolate" %in% names(asg)) as.character(asg$Isolate[1]) else "",
      locus = if (!is.null(asg) && nrow(asg) && "Locus" %in% names(asg)) as.character(asg$Locus[1]) else "",
      missing = !nzchar(seq)
    )
  }
  out
}

pitax_export_extract_qc_plot_targets <- function(snap, scope_ids = NULL) {
  ids <- names(snap$results)
  if (!is.null(scope_ids) && length(scope_ids)) ids <- intersect(ids, scope_ids)
  ids
}

pitax_export_extract_consensus_summary <- function(snap, scope_ids = NULL) {
  cs <- snap$consensus_set
  if (!is.list(cs) || !is.data.frame(cs$summary) || !nrow(cs$summary)) return(data.frame())
  df <- cs$summary
  if (!is.null(scope_ids) && length(scope_ids) && is.list(cs$records)) {
    keep <- character()
    for (id in names(cs$records)) {
      src <- as.character(cs$records[[id]]$source_read_ids)
      if (any(src %in% scope_ids) || id %in% scope_ids) keep <- c(keep, id)
    }
    if ("Consensus_ID" %in% names(df)) df <- df[df$Consensus_ID %in% keep, , drop = FALSE]
  }
  df
}

pitax_export_extract_consensus_analysis_sequences <- function(snap, scope_ids = NULL) {
  cs <- snap$consensus_set
  out <- list()
  if (!is.list(cs) || !is.list(cs$records)) return(out)
  for (id in names(cs$records)) {
    rec <- cs$records[[id]]
    src <- as.character(rec$source_read_ids)
    if (!is.null(scope_ids) && length(scope_ids) && !any(src %in% scope_ids) && !id %in% scope_ids) next
    seq <- if (is.null(rec$sequence)) "" else as.character(rec$sequence)[1]
    status <- if (is.null(rec$status)) "" else as.character(rec$status)[1]
    role <- if (identical(status, "INDEPENDENT_READ") || identical(status, "SINGLE_READ")) {
      "analysis_oriented"
    } else {
      "consensus"
    }
    out[[id]] <- list(
      id = id,
      final_name = if (nzchar(as.character(rec$final_name %||% ""))) as.character(rec$final_name) else id,
      sequence = seq,
      sequence_role = role,
      direction = "oriented",
      isolate = as.character(rec$isolate %||% ""),
      locus = as.character(rec$locus %||% ""),
      status = status,
      missing = !nzchar(seq),
      approved_final = identical(status, "READY") || identical(status, "INDEPENDENT_READ") || identical(status, "SINGLE_READ")
    )
  }
  out
}

pitax_export_extract_consensus_column_evidence <- function(snap, scope_ids = NULL) {
  cs <- snap$consensus_set
  rows <- list()
  if (!is.list(cs) || !is.list(cs$records)) return(data.frame())
  for (id in names(cs$records)) {
    rec <- cs$records[[id]]
    src <- as.character(rec$source_read_ids)
    if (!is.null(scope_ids) && length(scope_ids) && !any(src %in% scope_ids) && !id %in% scope_ids) next
    ev <- rec$evidence
    if (!is.data.frame(ev) || !nrow(ev)) next
    ev$Consensus_ID <- id
    ev$Final_Name <- as.character(rec$final_name %||% id)
    rows[[length(rows) + 1L]] <- ev
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pitax_export_extract_consensus_conflicts <- function(snap, scope_ids = NULL) {
  ev <- pitax_export_extract_consensus_column_evidence(snap, scope_ids)
  if (!nrow(ev)) return(data.frame())
  if ("Needs_Review" %in% names(ev)) {
    keep <- as.logical(ev$Needs_Review)
    keep[is.na(keep)] <- FALSE
    return(ev[keep, , drop = FALSE])
  }
  if ("Decision" %in% names(ev)) {
    return(ev[grepl("Conflict|IUPAC|review", as.character(ev$Decision), ignore.case = TRUE), , drop = FALSE])
  }
  ev[0, , drop = FALSE]
}

pitax_export_extract_consensus_pairwise_alignments <- function(snap, scope_ids = NULL) {
  cs <- snap$consensus_set
  out <- list()
  if (!is.list(cs) || !is.list(cs$records)) return(out)
  for (id in names(cs$records)) {
    rec <- cs$records[[id]]
    if (!is.list(rec$alignment)) next
    src <- as.character(rec$source_read_ids)
    if (!is.null(scope_ids) && length(scope_ids) && !any(src %in% scope_ids) && !id %in% scope_ids) next
    out[[id]] <- list(
      id = id,
      final_name = as.character(rec$final_name %||% id),
      text = paste(
        paste0("consensus_id: ", id),
        paste0("final_name: ", as.character(rec$final_name %||% "")),
        paste0("status: ", as.character(rec$status %||% "")),
        paste0("Forward: ", as.character(rec$alignment$forward_aligned %||% "")),
        paste0("Reverse: ", as.character(rec$alignment$reverse_aligned %||% "")),
        paste0("Consensus: ", as.character(rec$sequence %||% "")),
        sep = "\n"
      )
    )
  }
  out
}

pitax_export_extract_consensus_active_curation <- function(snap, scope_ids = NULL) {
  cs <- snap$consensus_set
  rows <- list()
  if (!is.list(cs) || !is.list(cs$records)) return(data.frame())
  for (id in names(cs$records)) {
    rec <- cs$records[[id]]
    src <- as.character(rec$source_read_ids)
    if (!is.null(scope_ids) && length(scope_ids) && !any(src %in% scope_ids) && !id %in% scope_ids) next
    log <- if (is.list(rec$curation)) rec$curation$audit_log else NULL
    if (!is.data.frame(log) || !nrow(log)) next
    log$Consensus_ID <- id
    log$Export_class <- "active_curation"
    rows[[length(rows) + 1L]] <- log
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pitax_export_extract_consensus_historical_audit <- function(snap, scope_ids = NULL) {
  df <- pitax_export_extract_consensus_active_curation(snap, scope_ids)
  if (!nrow(df)) return(data.frame())
  df$Export_class <- "historical_audit"
  df
}

pitax_export_extract_blast_jobs <- function(snap, scope_ids = NULL) {
  df <- snap$blast_jobs
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    df <- df[df$original_name %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  df
}

pitax_export_analysis_ids_for_sources <- function(snap, scope_ids) {
  cs <- snap$consensus_set
  ids <- character()
  if (is.list(cs) && is.list(cs$records)) {
    for (id in names(cs$records)) {
      src <- as.character(cs$records[[id]]$source_read_ids)
      if (any(src %in% scope_ids) || id %in% scope_ids) ids <- c(ids, id)
    }
  }
  unique(c(ids, scope_ids))
}

pitax_export_extract_blast_hits <- function(snap, scope_ids = NULL) {
  df <- snap$blast_hits
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  jobs <- pitax_export_extract_blast_jobs(snap, scope_ids)
  if (nrow(jobs) && "rid" %in% names(df)) {
    active_rids <- as.character(jobs$rid[vapply(seq_len(nrow(jobs)), function(i) {
      pitax_export_blast_job_is_active(snap, jobs[i, , drop = FALSE])
    }, logical(1))])
    df <- df[df$rid %in% active_rids, , drop = FALSE]
  }
  if (!is.null(scope_ids) && length(scope_ids) && "original_name" %in% names(df)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    df <- df[df$original_name %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  df
}

pitax_export_extract_blast_raw <- function(snap, scope_ids = NULL) {
  jobs <- pitax_export_extract_blast_jobs(snap, scope_ids)
  out <- list()
  if (!nrow(jobs) || !length(snap$blast_raw)) return(out)
  for (i in seq_len(nrow(jobs))) {
    rid <- as.character(jobs$rid[i])
    if (!nzchar(rid) || is.null(snap$blast_raw[[rid]])) next
    out[[rid]] <- list(
      id = rid,
      final_name = as.character(jobs$final_name[i]),
      text = as.character(snap$blast_raw[[rid]])
    )
  }
  out
}

pitax_export_extract_blast_query_sequences <- function(snap, scope_ids = NULL) {
  pitax_export_extract_consensus_analysis_sequences(snap, scope_ids)
}

pitax_export_extract_taxonomy_summary <- function(snap, scope_ids = NULL) {
  df <- snap$taxonomy_summary
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "original_name" %in% names(df)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    df <- df[df$original_name %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  df
}

pitax_export_extract_taxonomy_hits <- function(snap, scope_ids = NULL) {
  df <- snap$taxonomy_hits
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "original_name" %in% names(df)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    df <- df[df$original_name %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  df
}

pitax_export_extract_taxonomy_species_evidence <- function(snap, scope_ids = NULL) {
  df <- snap$taxonomy_counts
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  if (!is.null(scope_ids) && length(scope_ids) && "original_name" %in% names(df)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    df <- df[df$original_name %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  df
}

pitax_export_extract_taxonomy_team_summary <- function(snap, scope_ids = NULL) {
  # Minimal team summary from already-computed taxonomy + analysis summary only.
  analysis <- tryCatch(stage3_analysis_summary(snap$consensus_set), error = function(e) data.frame())
  tax <- snap$taxonomy_summary
  if (!is.data.frame(analysis) || !nrow(analysis)) return(data.frame())
  out <- data.frame(
    Sample = as.character(analysis$final_name),
    Original_sample = as.character(analysis$sample_id),
    Target = if ("target" %in% names(analysis)) as.character(analysis$target) else "",
    Trimmed_length = if ("trimmed_length" %in% names(analysis)) analysis$trimmed_length else NA_integer_,
    Identification = NA_character_,
    Identification_level = NA_character_,
    Overall_confidence = NA_character_,
    Best_molecular_match = NA_character_,
    RID = NA_character_,
    Comment = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.data.frame(tax) && nrow(tax) && "original_name" %in% names(tax)) {
    for (i in seq_len(nrow(out))) {
      j <- which(tax$original_name == out$Original_sample[i])
      if (!length(j)) next
      j <- j[1]
      if ("recommended_identification" %in% names(tax)) out$Identification[i] <- as.character(tax$recommended_identification[j])
      if ("recommended_level" %in% names(tax)) out$Identification_level[i] <- as.character(tax$recommended_level[j])
      if ("confidence" %in% names(tax)) out$Overall_confidence[i] <- as.character(tax$confidence[j])
      if ("best_molecular_match" %in% names(tax)) out$Best_molecular_match[i] <- as.character(tax$best_molecular_match[j])
      if ("rid" %in% names(tax)) out$RID[i] <- as.character(tax$rid[j])
      if ("decision_reason" %in% names(tax)) out$Comment[i] <- as.character(tax$decision_reason[j])
    }
  }
  if (!is.null(scope_ids) && length(scope_ids)) {
    analysis_ids <- pitax_export_analysis_ids_for_sources(snap, scope_ids)
    out <- out[out$Original_sample %in% c(scope_ids, analysis_ids), , drop = FALSE]
  }
  out
}

pitax_export_extract_multilocus_profiles <- function(snap, scope_ids = NULL) {
  prof <- snap$multilocus_profile
  if (!is.list(prof) || !is.data.frame(prof$profiles)) return(data.frame())
  df <- prof$profiles
  isolates <- unique(as.character(selection_isolates_from_scope(snap, scope_ids)))
  if (length(isolates) && "Isolate" %in% names(df)) df <- df[df$Isolate %in% isolates, , drop = FALSE]
  df
}

selection_isolates_from_scope <- function(snap, scope_ids) {
  if (is.null(scope_ids) || !length(scope_ids)) return(character())
  asg <- snap$read_assignments
  if (!is.data.frame(asg) || !nrow(asg)) return(character())
  unique(as.character(asg$Isolate[asg$Source_ID %in% scope_ids]))
}

pitax_export_extract_multilocus_evidence <- function(snap, scope_ids = NULL) {
  prof <- snap$multilocus_profile
  if (!is.list(prof) || !is.data.frame(prof$evidence)) return(data.frame())
  df <- prof$evidence
  isolates <- unique(as.character(selection_isolates_from_scope(snap, scope_ids)))
  if (length(isolates) && "Isolate" %in% names(df)) df <- df[df$Isolate %in% isolates, , drop = FALSE]
  df
}

pitax_export_extract_multilocus_sources <- function(snap, scope_ids = NULL) {
  prof <- snap$multilocus_profile
  if (!is.list(prof) || !is.data.frame(prof$sources)) return(data.frame())
  prof$sources
}

pitax_export_extract_multilocus_sequences <- function(snap, scope_ids = NULL) {
  ev <- pitax_export_extract_multilocus_evidence(snap, scope_ids)
  out <- list()
  if (!nrow(ev) || !"Sequence" %in% names(ev)) return(out)
  for (i in seq_len(nrow(ev))) {
    id <- paste(ev$Isolate[i], ev$Locus[i], sep = "__")
    seq <- as.character(ev$Sequence[i])
    out[[id]] <- list(
      id = id,
      final_name = paste(ev$Isolate[i], ev$Locus[i], sep = "_"),
      sequence = seq,
      sequence_role = "multilocus_locus",
      direction = "oriented",
      isolate = as.character(ev$Isolate[i]),
      locus = as.character(ev$Locus[i]),
      missing = !nzchar(seq)
    )
  }
  out
}

pitax_export_extractor <- function(component_id) {
  switch(
    as.character(component_id)[1],
    "upload_file_inventory" = pitax_export_extract_upload_file_inventory,
    "upload_raw_ab1" = pitax_export_extract_upload_raw_ab1,
    "assay_profiles_draft" = pitax_export_extract_assay_profiles_draft,
    "project_defaults_draft" = pitax_export_extract_project_defaults_draft,
    "processing_settings_applied" = pitax_export_extract_processing_settings_applied,
    "assign_identity_table" = pitax_export_extract_assign_identity_table,
    "assign_rename_map" = pitax_export_extract_assign_rename_map,
    "assign_architecture" = pitax_export_extract_assign_architecture,
    "qc_summary" = pitax_export_extract_qc_summary,
    "qc_ab1_run_evidence" = pitax_export_extract_qc_ab1_run_evidence,
    "qc_ab1_base_evidence" = pitax_export_extract_qc_ab1_base_evidence,
    "qc_peak_flags" = pitax_export_extract_qc_peak_flags,
    "qc_active_curation" = pitax_export_extract_qc_active_curation,
    "qc_historical_audit" = pitax_export_extract_qc_historical_audit,
    "qc_trimmed_sequences" = pitax_export_extract_qc_trimmed_sequences,
    "qc_metrics_plots" = pitax_export_extract_qc_plot_targets,
    "qc_chromatograms" = pitax_export_extract_qc_plot_targets,
    "consensus_summary" = pitax_export_extract_consensus_summary,
    "consensus_analysis_sequences" = pitax_export_extract_consensus_analysis_sequences,
    "consensus_column_evidence" = pitax_export_extract_consensus_column_evidence,
    "consensus_pairwise_alignments" = pitax_export_extract_consensus_pairwise_alignments,
    "consensus_conflicts" = pitax_export_extract_consensus_conflicts,
    "consensus_active_curation" = pitax_export_extract_consensus_active_curation,
    "consensus_historical_audit" = pitax_export_extract_consensus_historical_audit,
    "blast_jobs" = pitax_export_extract_blast_jobs,
    "blast_hits" = pitax_export_extract_blast_hits,
    "blast_raw" = pitax_export_extract_blast_raw,
    "blast_query_sequences" = pitax_export_extract_blast_query_sequences,
    "taxonomy_team_summary" = pitax_export_extract_taxonomy_team_summary,
    "taxonomy_summary" = pitax_export_extract_taxonomy_summary,
    "taxonomy_hits" = pitax_export_extract_taxonomy_hits,
    "taxonomy_species_evidence" = pitax_export_extract_taxonomy_species_evidence,
    "taxonomy_score_plot" = function(snap, scope_ids = NULL) {
      df <- pitax_export_extract_taxonomy_hits(snap, scope_ids)
      if (!nrow(df) || !"original_name" %in% names(df)) return(character())
      unique(as.character(df$original_name))
    },
    "multilocus_profiles" = pitax_export_extract_multilocus_profiles,
    "multilocus_evidence" = pitax_export_extract_multilocus_evidence,
    "multilocus_sources" = pitax_export_extract_multilocus_sources,
    "multilocus_sequences" = pitax_export_extract_multilocus_sequences,
    "multilocus_plot" = function(snap, scope_ids = NULL) {
      df <- pitax_export_extract_multilocus_profiles(snap, scope_ids)
      if (!nrow(df) || !"Isolate" %in% names(df)) return(character())
      unique(as.character(df$Isolate))
    },
    NULL
  )
}

pitax_export_assess_component <- function(component_id, snap, selection = list()) {
  comp <- pitax_export_get_component(component_id)
  if (is.null(comp)) return(pitax_export_status("blocked", "Unknown export component."))
  scope_ids <- pitax_export_filter_scope_ids(snap, selection)
  extractor <- pitax_export_extractor(component_id)
  if (is.null(extractor)) return(pitax_export_status("blocked", "No extractor registered for component."))

  if (identical(component_id, "upload_file_inventory")) {
    df <- extractor(snap, scope_ids)
    if (!nrow(df)) return(pitax_export_status("not_created", "No uploaded files are present in this session."))
    return(pitax_export_status("available", "Upload inventory is available.", record_count = nrow(df)))
  }
  if (identical(component_id, "upload_raw_ab1")) {
    df <- extractor(snap, scope_ids)
    all_files <- pitax_export_extract_upload_file_inventory(snap, scope_ids)
    if (!nrow(all_files)) return(pitax_export_status("not_created", "No uploaded files are present."))
    if (!nrow(df)) {
      return(pitax_export_status("blocked", "Original AB1 bytes are no longer available. Raw files are not reconstructed from processed traces."))
    }
    if (nrow(df) < nrow(all_files)) {
      return(pitax_export_status("partial", "Original bytes are available for only some uploaded files.",
                                 ready_ids = as.character(df$Source_ID), record_count = nrow(df)))
    }
    return(pitax_export_status("available", "Original AB1 bytes are available.",
                               ready_ids = as.character(df$Source_ID), record_count = nrow(df)))
  }
  if (component_id %in% c("assay_profiles_draft", "project_defaults_draft")) {
    df <- extractor(snap, scope_ids)
    if (!nrow(df)) return(pitax_export_status("not_created", "Assay draft settings have not been defined."))
    return(pitax_export_status("available", "Draft assay settings are available.", record_count = nrow(df)))
  }
  if (identical(component_id, "processing_settings_applied")) {
    df <- extractor(snap, scope_ids)
    if (!nrow(df)) return(pitax_export_status("not_created", "No per-read processing settings have been stored yet. Run Trim & QC first."))
    return(pitax_export_status("available", "Applied processing settings are available.", record_count = nrow(df)))
  }
  if (component_id %in% c("assign_identity_table", "assign_rename_map", "assign_architecture")) {
    if (!pitax_export_assignment_ready(snap)) {
      return(pitax_export_status("blocked", "Read identity is incomplete or not unique."))
    }
    df <- extractor(snap, scope_ids)
    if (!nrow(df) && !identical(component_id, "assign_architecture")) {
      return(pitax_export_status("not_created", "No assignment rows are available."))
    }
    if (identical(component_id, "assign_architecture") && (!is.list(snap$architecture) || !length(snap$architecture))) {
      return(pitax_export_status("not_created", "Project architecture has not been built."))
    }
    n <- if (is.data.frame(df)) nrow(df) else 0L
    return(pitax_export_status("available", "Assignment data is available before trimming.", record_count = n))
  }
  if (startsWith(component_id, "qc_")) {
    if (!length(snap$results)) return(pitax_export_status("not_created", "Trim & QC has not been run."))
    data <- extractor(snap, scope_ids)
    n <- if (is.data.frame(data)) nrow(data) else if (is.list(data)) length(data) else length(data)
    if (!n) return(pitax_export_status("not_created", "No QC records match the current selection."))
    selected_missing <- setdiff(scope_ids, names(snap$results))
    if (length(scope_ids) && length(selected_missing) && length(intersect(scope_ids, names(snap$results)))) {
      return(pitax_export_status("partial", "Some selected samples are not yet processed.",
                                 ready_ids = intersect(scope_ids, names(snap$results)),
                                 blocked_ids = selected_missing, record_count = n))
    }
    return(pitax_export_status("available", "QC data is available for the selection.", record_count = n))
  }
  if (startsWith(component_id, "consensus_")) {
    cs <- snap$consensus_set
    if (!is.list(cs) || !length(cs$records)) {
      return(pitax_export_status("not_created", "Analysis sequences have not been built."))
    }
    if (!isTRUE(tryCatch(stage3_consensus_set_is_current(cs, snap$results), error = function(e) FALSE))) {
      return(pitax_export_status("stale", "Analysis sequences are not current for the curated source reads."))
    }
    if (identical(component_id, "consensus_pairwise_alignments") && identical(snap$project_mode, "simple")) {
      return(pitax_export_status("not_created", "Pairwise alignments exist only in paired Forward/Reverse mode."))
    }
    data <- extractor(snap, scope_ids)
    n <- if (is.data.frame(data)) nrow(data) else length(data)
    if (!n) return(pitax_export_status("not_created", "No consensus records match the current selection."))
    if (identical(comp$approval, "approved_result") && is.list(data) && !is.data.frame(data)) {
      blocked <- names(Filter(function(x) is.list(x) && isFALSE(x$approved_final), data))
      ready <- setdiff(names(data), blocked)
      if (length(blocked) && !length(ready)) {
        return(pitax_export_status("blocked", "Selected analysis sequences are not approved as final (review required or no reliable overlap).",
                                   blocked_ids = blocked, record_count = n))
      }
      if (length(blocked)) {
        return(pitax_export_status("partial", "Some analysis sequences remain under review and will not be labeled final.",
                                   ready_ids = ready, blocked_ids = blocked, record_count = length(ready)))
      }
    }
    return(pitax_export_status("available", "Consensus data is available.", record_count = n))
  }
  if (startsWith(component_id, "blast_")) {
    jobs <- snap$blast_jobs
    if (!is.data.frame(jobs) || !nrow(jobs)) {
      return(pitax_export_status("not_created", "No BLAST jobs have been created."))
    }
    data <- extractor(snap, scope_ids)
    n <- if (is.data.frame(data)) nrow(data) else length(data)
    if (!n) {
      if (identical(component_id, "blast_hits")) {
        scoped_jobs <- pitax_export_extract_blast_jobs(snap, scope_ids)
        if (nrow(scoped_jobs)) {
          statuses <- unique(toupper(as.character(scoped_jobs$status)))
          if (any(statuses %in% c("WAITING", "SUBMITTED"))) {
            return(pitax_export_status("partial", "BLAST jobs are still waiting; no hits are available yet.", record_count = 0L))
          }
          if (any(statuses %in% c("FAILED", "ERROR", "TIMEOUT"))) {
            return(pitax_export_status("partial", "BLAST jobs failed or timed out; no hits are available.", record_count = 0L))
          }
          return(pitax_export_status("partial", "BLAST completed without hits for the selection.", record_count = 0L))
        }
      }
      return(pitax_export_status("not_created", "No BLAST records match the current selection."))
    }
    if (is.data.frame(data) && "status" %in% names(data)) {
      stale <- toupper(as.character(data$status)) == "STALE"
      if (all(stale)) return(pitax_export_status("stale", "All matching BLAST jobs are stale for the active sequence revision."))
      if (any(stale)) {
        return(pitax_export_status("partial", "Some BLAST jobs are stale and are excluded from active results.",
                                   ready_ids = as.character(data$rid[!stale]),
                                   blocked_ids = as.character(data$rid[stale]),
                                   record_count = sum(!stale)))
      }
    }
    return(pitax_export_status("available", "BLAST data is available for the selection.", record_count = n))
  }
  if (startsWith(component_id, "taxonomy_")) {
    if (!is.data.frame(snap$taxonomy_summary) || !nrow(snap$taxonomy_summary)) {
      return(pitax_export_status("not_created", "Taxonomic analysis has not been computed."))
    }
    data <- extractor(snap, scope_ids)
    n <- if (is.data.frame(data)) nrow(data) else length(data)
    if (!n) return(pitax_export_status("not_created", "No taxonomic records match the current selection."))
    return(pitax_export_status("available", "Taxonomy data is available.", record_count = n))
  }
  if (startsWith(component_id, "multilocus_")) {
    prof <- snap$multilocus_profile
    if (!is.list(prof) || !is.data.frame(prof$profiles) || !nrow(prof$profiles)) {
      return(pitax_export_status("not_created", "Multi-locus profile has not been built."))
    }
    data <- extractor(snap, scope_ids)
    n <- if (is.data.frame(data)) nrow(data) else length(data)
    if (!n) return(pitax_export_status("not_created", "No multi-locus records match the current selection."))
    return(pitax_export_status("available", "Multi-locus data is available.", record_count = n))
  }
  pitax_export_status("blocked", "Unhandled component availability path.")
}

pitax_export_assess_all <- function(snap, selection = list(), component_ids = NULL) {
  ids <- if (is.null(component_ids)) pitax_export_component_ids() else as.character(component_ids)
  out <- lapply(ids, function(id) pitax_export_assess_component(id, snap, selection))
  names(out) <- ids
  out
}
