# ============================================================
# Export plan builder (no file writing)
# ============================================================

pitax_export_default_selection <- function(stage_id = "export", scope = c("stage", "project"),
                                           sample_mode = c("all", "current", "selected"),
                                           selected_samples = character(),
                                           isolates = character(), loci = character(), assays = character()) {
  scope <- match.arg(scope)
  sample_mode <- match.arg(sample_mode)
  list(
    stage_id = as.character(stage_id)[1],
    scope = scope,
    sample_mode = sample_mode,
    selected_samples = unique(as.character(selected_samples)),
    isolates = unique(as.character(isolates)),
    loci = unique(as.character(loci)),
    assays = unique(as.character(assays)),
    component_ids = pitax_export_default_component_ids(stage_id, scope),
    row_mode = "all",
    selected_rows = integer(),
    column_mode = "display",
    column_ids = character(),
    excel_layout = "sheet_per_component",
    sequence_layout = "combined",
    sequence_header = "short",
    include_plots = TRUE,
    package_name = "PITAX_export_package.zip"
  )
}

pitax_export_selectable_statuses <- function() c("available", "partial")

pitax_export_apply_row_column_filters <- function(df, comp, selection) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  row_mode <- as.character(selection$row_mode %||% "all")[1]
  if (identical(row_mode, "selected") && length(selection$selected_rows)) {
    idx <- as.integer(selection$selected_rows)
    idx <- idx[idx >= 1L & idx <= nrow(df)]
    df <- df[idx, , drop = FALSE]
  }
  cols <- names(df)
  column_mode <- as.character(selection$column_mode %||% "display")[1]
  if (identical(column_mode, "custom") && length(selection$column_ids)) {
    keep <- intersect(as.character(selection$column_ids), cols)
    if (length(keep)) df <- df[, keep, drop = FALSE]
  } else if (identical(column_mode, "display") && length(comp$columns)) {
    preferred <- vapply(comp$columns, function(c) c$id, character(1))
    preferred <- preferred[vapply(comp$columns, function(c) isTRUE(c$display_default), logical(1))]
    keep <- intersect(preferred, cols)
    if (length(keep)) df <- df[, keep, drop = FALSE]
  }
  df
}

pitax_export_planned_file <- function(path, component_id, record_ids = character(), kind = "table",
                                      optional = FALSE) {
  list(
    path = as.character(path)[1],
    component_id = as.character(component_id)[1],
    record_ids = unique(as.character(record_ids)),
    kind = as.character(kind)[1],
    optional = isTRUE(optional)
  )
}

pitax_export_build_plan <- function(snap, selection = pitax_export_default_selection()) {
  selection <- modifyList(pitax_export_default_selection(selection$stage_id %||% "export",
                                                         selection$scope %||% "stage"), selection)
  availability <- pitax_export_assess_all(snap, selection)
  requested <- unique(as.character(selection$component_ids))
  included <- character()
  omitted <- list()
  files <- list()
  warnings <- character()
  pending <- snap$pending_curation
  if (!is.null(pending)) {
    warnings <- c(warnings, "Unconfirmed curation edits are present and will NOT be applied during export.")
  }
  if (is.data.frame(snap$auto_correct_preview) && nrow(snap$auto_correct_preview)) {
    warnings <- c(warnings, "Auto-correction preview rows are present and will NOT be applied during export.")
  }

  excel_components <- list()
  sequence_components <- list()
  plot_components <- list()
  text_components <- list()
  raw_components <- list()

  for (id in requested) {
    av <- availability[[id]]
    if (is.null(av)) {
      omitted[[id]] <- list(reason = "Unknown component", status = "blocked")
      next
    }
    if (!av$status %in% pitax_export_selectable_statuses()) {
      omitted[[id]] <- list(reason = av$reason, status = av$status)
      next
    }
    if (identical(av$status, "partial") && identical(av$record_count, 0L) &&
        id %in% c("blast_hits")) {
      omitted[[id]] <- list(reason = av$reason, status = av$status)
      next
    }
    included <- c(included, id)
    comp <- pitax_export_get_component(id)
    if (identical(comp$kind, "table")) excel_components[[id]] <- comp
    else if (identical(comp$kind, "sequence")) sequence_components[[id]] <- comp
    else if (identical(comp$kind, "plot")) plot_components[[id]] <- comp
    else if (identical(comp$kind, "text")) text_components[[id]] <- comp
    else if (identical(comp$kind, "raw_file")) raw_components[[id]] <- comp
  }

  layout <- as.character(selection$excel_layout %||% "sheet_per_component")[1]
  if (length(excel_components)) {
    if (identical(layout, "workbook_per_isolate")) {
      isolates <- unique(as.character(selection$isolates))
      if (!length(isolates) && is.data.frame(snap$read_assignments) && nrow(snap$read_assignments)) {
        isolates <- unique(as.character(snap$read_assignments$Isolate))
        isolates <- isolates[nzchar(isolates)]
      }
      if (!length(isolates)) isolates <- "ALL"
      for (iso in isolates) {
        files[[length(files) + 1L]] <- pitax_export_planned_file(
          paste0("excel/", clean_fasta_name(iso), "_export.xlsx"),
          paste(names(excel_components), collapse = ","),
          record_ids = iso,
          kind = "excel"
        )
      }
    } else {
      files[[length(files) + 1L]] <- pitax_export_planned_file(
        "excel/PITAX_export_tables.xlsx",
        paste(names(excel_components), collapse = ","),
        kind = "excel"
      )
    }
  }

  seq_layout <- as.character(selection$sequence_layout %||% "combined")[1]
  for (id in names(sequence_components)) {
    extractor <- pitax_export_extractor(id)
    scope_ids <- pitax_export_filter_scope_ids(snap, selection)
    seqs <- extractor(snap, scope_ids)
    if (!length(seqs)) next
    if (identical(seq_layout, "combined")) {
      files[[length(files) + 1L]] <- pitax_export_planned_file(
        paste0("sequences/", id, ".fasta"), id, names(seqs), "fasta"
      )
    } else if (identical(seq_layout, "per_record")) {
      for (nm in names(seqs)) {
        stem <- clean_fasta_name(seqs[[nm]]$final_name)
        files[[length(files) + 1L]] <- pitax_export_planned_file(
          paste0("sequences/", id, "/", stem, ".fasta"), id, nm, "fasta"
        )
      }
    } else if (identical(seq_layout, "per_isolate")) {
      by_iso <- split(names(seqs), vapply(seqs, function(x) x$isolate %||% "unknown", character(1)))
      for (iso in names(by_iso)) {
        files[[length(files) + 1L]] <- pitax_export_planned_file(
          paste0("sequences/", id, "/", clean_fasta_name(iso), ".fasta"), id, by_iso[[iso]], "fasta"
        )
      }
    } else if (identical(seq_layout, "per_locus")) {
      by_loc <- split(names(seqs), vapply(seqs, function(x) x$locus %||% "unknown", character(1)))
      for (loc in names(by_loc)) {
        files[[length(files) + 1L]] <- pitax_export_planned_file(
          paste0("sequences/", id, "/", clean_fasta_name(loc), ".fasta"), id, by_loc[[loc]], "fasta"
        )
      }
    }
    missing_ids <- names(Filter(function(x) isTRUE(x$missing), seqs))
    if (length(missing_ids)) {
      files[[length(files) + 1L]] <- pitax_export_planned_file(
        paste0("sequences/", id, "_missing_sequences.txt"), id, missing_ids, "text"
      )
    }
  }

  if (isTRUE(selection$include_plots)) {
    for (id in names(plot_components)) {
      extractor <- pitax_export_extractor(id)
      scope_ids <- pitax_export_filter_scope_ids(snap, selection)
      targets <- extractor(snap, scope_ids)
      for (tg in targets) {
        files[[length(files) + 1L]] <- pitax_export_planned_file(
          paste0("plots/", id, "/", clean_fasta_name(tg), ".png"),
          id, tg, "plot", optional = TRUE
        )
      }
    }
  }

  for (id in names(text_components)) {
    extractor <- pitax_export_extractor(id)
    scope_ids <- pitax_export_filter_scope_ids(snap, selection)
    items <- extractor(snap, scope_ids)
    for (nm in names(items)) {
      stem <- clean_fasta_name(items[[nm]]$final_name %||% nm)
      files[[length(files) + 1L]] <- pitax_export_planned_file(
        paste0("text/", id, "/", stem, ".txt"), id, nm, "text"
      )
    }
  }

  for (id in names(raw_components)) {
    extractor <- pitax_export_extractor(id)
    scope_ids <- pitax_export_filter_scope_ids(snap, selection)
    df <- extractor(snap, scope_ids)
    if (!is.data.frame(df) || !nrow(df)) next
    for (i in seq_len(nrow(df))) {
      files[[length(files) + 1L]] <- pitax_export_planned_file(
        paste0("raw_ab1/", as.character(df$name[i])),
        id, as.character(df$Source_ID[i]), "raw_file"
      )
    }
  }

  files[[length(files) + 1L]] <- pitax_export_planned_file("RUN_INFO.txt", "run_info", kind = "text")

  preview_rows <- list()
  for (id in included) {
    av <- availability[[id]]
    preview_rows[[length(preview_rows) + 1L]] <- data.frame(
      Component = id,
      Label = pitax_export_get_component(id)$label,
      Status = av$status,
      Records = av$record_count,
      Reason = av$reason,
      stringsAsFactors = FALSE
    )
  }
  preview <- if (length(preview_rows)) do.call(rbind, preview_rows) else data.frame(
    Component = character(), Label = character(), Status = character(),
    Records = integer(), Reason = character(), stringsAsFactors = FALSE
  )

  list(
    snapshot_id = snap$snapshot_id,
    selection = selection,
    availability = availability,
    included_components = unique(included),
    omitted_components = omitted,
    warnings = warnings,
    files = files,
    preview = preview,
    excel_components = names(excel_components),
    sequence_components = names(sequence_components),
    plot_components = names(plot_components),
    text_components = names(text_components),
    raw_components = names(raw_components)
  )
}

pitax_export_validate_plan <- function(snap, plan) {
  if (!identical(as.character(snap$snapshot_id), as.character(plan$snapshot_id))) {
    return(list(ok = FALSE, error = "Export snapshot changed. Re-open Export and review the selection."))
  }
  fresh <- pitax_export_build_plan(snap, plan$selection)
  for (id in plan$included_components) {
    av <- fresh$availability[[id]]
    if (is.null(av) || !av$status %in% pitax_export_selectable_statuses()) {
      return(list(ok = FALSE, error = paste0("Component '", id, "' is no longer exportable: ", av$reason %||% "unknown")))
    }
  }
  list(ok = TRUE, plan = fresh)
}
