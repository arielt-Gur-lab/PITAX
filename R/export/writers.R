# ============================================================
# Export writers: Excel, FASTA, PNG, RUN_INFO, ZIP
# ============================================================

pitax_export_safe_sheet_name <- function(name, used = character()) {
  name <- as.character(name)[1]
  name <- gsub("[\\\\/:?*\\[\\]]", "_", name)
  name <- trimws(name)
  if (!nzchar(name)) name <- "Sheet"
  if (nchar(name) > 31) name <- substr(name, 1, 31)
  base <- name
  i <- 1L
  while (tolower(name) %in% tolower(used)) {
    suffix <- paste0("_", i)
    name <- paste0(substr(base, 1, max(1, 31 - nchar(suffix))), suffix)
    i <- i + 1L
  }
  name
}

pitax_export_coerce_column_types <- function(df, comp = NULL) {
  if (!is.data.frame(df) || !ncol(df)) return(df)
  type_map <- list()
  if (!is.null(comp) && length(comp$columns)) {
    for (col in comp$columns) type_map[[col$id]] <- col$type
  }
  for (nm in names(df)) {
    tp <- type_map[[nm]]
    if (identical(tp, "number")) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    } else if (identical(tp, "logical")) {
      df[[nm]] <- as.logical(df[[nm]])
    } else {
      # Keep identifiers as text even when they look numeric.
      if (identical(tp, "text") || nm %in% c(comp$link_keys, "rid", "accession", "Source_ID", "Sample", "Isolate")) {
        df[[nm]] <- as.character(df[[nm]])
      }
    }
  }
  df
}

pitax_export_write_excel_workbook <- function(path, snap, plan, isolate_filter = NULL) {
  comps <- plan$excel_components
  if (!length(comps)) return(invisible(NULL))
  wb <- openxlsx::createWorkbook(creator = "PITAX")
  header_style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold",
                                        halign = "center", valign = "center")
  wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  used_names <- character()
  layout <- as.character(plan$selection$excel_layout %||% "sheet_per_component")[1]
  scope_ids <- pitax_export_filter_scope_ids(snap, plan$selection)
  if (!is.null(isolate_filter) && nzchar(isolate_filter) && !identical(isolate_filter, "ALL")) {
    asg <- snap$read_assignments
    if (is.data.frame(asg) && nrow(asg)) {
      scope_ids <- intersect(scope_ids, as.character(asg$Source_ID[asg$Isolate == isolate_filter]))
    }
  }

  write_sheet <- function(sheet_name, df, comp) {
    sheet_name <- pitax_export_safe_sheet_name(sheet_name, used_names)
    used_names <<- c(used_names, sheet_name)
    if (is.null(df) || !is.data.frame(df) || !ncol(df)) {
      df <- data.frame(Message = "No data available", stringsAsFactors = FALSE)
    }
    df <- pitax_export_apply_row_column_filters(df, comp, plan$selection)
    df <- pitax_export_coerce_column_types(df, comp)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df, withFilter = nrow(df) > 0, headerStyle = header_style)
    openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
    if (nrow(df)) {
      openxlsx::addStyle(wb, sheet_name, wrap_style, rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)),
                         gridExpand = TRUE, stack = TRUE)
    }
    openxlsx::setColWidths(wb, sheet_name, cols = seq_len(max(1, ncol(df))), widths = "auto")
    invisible(sheet_name)
  }

  if (identical(layout, "single_sheet") && length(comps) == 1L) {
    id <- comps[1]
    comp <- pitax_export_get_component(id)
    df <- pitax_export_extractor(id)(snap, scope_ids)
    if (!is.data.frame(df)) stop("Component ", id, " did not return a table.", call. = FALSE)
    write_sheet(comp$label, df, comp)
  } else if (identical(layout, "sheet_per_sample")) {
    for (id in comps) {
      comp <- pitax_export_get_component(id)
      df <- pitax_export_extractor(id)(snap, scope_ids)
      if (!is.data.frame(df) || !nrow(df)) next
      sample_col <- intersect(c("sample_id", "Sample", "Sample_ID", "Source_ID", "original_name", "Original_sample"), names(df))
      if (!length(sample_col)) {
        write_sheet(comp$label, df, comp)
        next
      }
      sc <- sample_col[1]
      for (sid in unique(as.character(df[[sc]]))) {
        write_sheet(paste0(comp$id, "_", sid), df[df[[sc]] == sid, , drop = FALSE], comp)
      }
    }
  } else {
    # sheet_per_component and workbook_per_isolate content path
    for (id in comps) {
      comp <- pitax_export_get_component(id)
      df <- pitax_export_extractor(id)(snap, scope_ids)
      if (!is.data.frame(df)) stop("Component ", id, " did not return a table.", call. = FALSE)
      write_sheet(comp$label, df, comp)
    }
  }

  # Provenance sheet with link keys for active tables.
  link_rows <- list()
  for (id in comps) {
    av <- plan$availability[[id]]
    link_rows[[length(link_rows) + 1L]] <- data.frame(
      Component = id,
      Status = av$status,
      Records = av$record_count,
      Link_keys = paste(pitax_export_get_component(id)$link_keys, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  write_sheet("Provenance", do.call(rbind, link_rows), NULL)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

pitax_export_write_fasta_records <- function(path, records, header_mode = c("short", "rich")) {
  header_mode <- match.arg(header_mode)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  used <- character()
  lines <- character()
  missing_lines <- character()
  for (nm in names(records)) {
    rec <- records[[nm]]
    if (isTRUE(rec$missing) || !nzchar(as.character(rec$sequence %||% ""))) {
      missing_lines <- c(missing_lines, paste0(nm, "\t", rec$final_name %||% "", "\tmissing_sequence"))
      next
    }
    base <- clean_fasta_name(rec$final_name %||% nm)
    header_id <- base
    if (tolower(header_id) %in% tolower(used)) {
      header_id <- clean_fasta_name(paste0(base, "__", nm))
    }
    used <- c(used, header_id)
    header <- paste0(">", header_id)
    if (identical(header_mode, "rich")) {
      header <- paste0(
        header,
        " sequence_role=", clean_fasta_name(rec$sequence_role %||% "unknown"),
        " direction=", clean_fasta_name(rec$direction %||% "unknown"),
        " isolate=", clean_fasta_name(rec$isolate %||% ""),
        " locus=", clean_fasta_name(rec$locus %||% ""),
        " source_id=", clean_fasta_name(nm)
      )
      if (!is.null(rec$status)) header <- paste0(header, " status=", clean_fasta_name(rec$status))
    }
    lines <- c(lines, header, wrap_sequence(rec$sequence, 80))
  }
  writeLines(paste(lines, collapse = "\n"), path)
  list(path = path, missing_lines = missing_lines, headers = used)
}

pitax_export_write_chromatogram_png <- function(result, settings, file) {
  if (is.null(result$peak_pos) || !length(result$peak_pos) || is.null(result$trace)) {
    stop("Chromatogram data is missing for this read.")
  }
  n <- length(result$peak_pos)
  grDevices::png(file, width = 1800, height = 700, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_chromatogram(result, settings, start_base = 1L, visible_bases = max(50L, n))
  invisible(file)
}

pitax_export_write_taxonomy_score_png <- function(hits_df, file, sample_id) {
  grDevices::png(file, width = 1400, height = 900, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!is.data.frame(hits_df) || !nrow(hits_df)) {
    plot.new()
    title(main = paste0(sample_id, " - no taxonomy hits"))
    return(invisible(file))
  }
  df <- hits_df
  if ("analysis_rank" %in% names(df)) {
    x <- suppressWarnings(as.numeric(df$analysis_rank))
  } else if ("rank" %in% names(df)) {
    x <- suppressWarnings(as.numeric(df$rank))
  } else {
    x <- seq_len(nrow(df))
  }
  y <- if ("bit_score" %in% names(df)) suppressWarnings(as.numeric(df$bit_score)) else
    if ("bit_score_num" %in% names(df)) suppressWarnings(as.numeric(df$bit_score_num)) else rep(NA_real_, nrow(df))
  plot(x, y, type = "b", xlab = "Hit rank", ylab = "Bit score",
       main = paste0(sample_id, " - taxonomy score landscape"))
  invisible(file)
}

pitax_export_write_multilocus_png <- function(evidence_df, file, isolate) {
  grDevices::png(file, width = 1400, height = 900, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  df <- evidence_df
  if (!is.data.frame(df) || !nrow(df)) {
    plot.new()
    title(main = paste0(isolate, " - no multi-locus metrics"))
    return(invisible(file))
  }
  x <- if ("Best_Match_Identity" %in% names(df)) suppressWarnings(as.numeric(df$Best_Match_Identity)) else rep(NA_real_, nrow(df))
  y <- if ("Best_Match_Coverage" %in% names(df)) suppressWarnings(as.numeric(df$Best_Match_Coverage)) else rep(NA_real_, nrow(df))
  plot(x, y, xlab = "Best match identity (%)", ylab = "Best match query coverage (%)",
       main = paste0(isolate, " - multi-locus evidence"), pch = 19)
  if ("Locus" %in% names(df)) text(x, y, labels = as.character(df$Locus), pos = 3, cex = 0.8)
  invisible(file)
}

pitax_export_write_run_info <- function(path, snap, plan, written_files, failures = list()) {
  lines <- c(
    "PITAX export RUN_INFO",
    paste0("application_version: ", snap$app_version),
    paste0("project_schema_version: ", snap$schema_version),
    paste0("exported_at: ", format(snap$captured_at, "%Y-%m-%d %H:%M:%S")),
    paste0("timezone: ", snap$timezone),
    paste0("snapshot_id: ", snap$snapshot_id),
    paste0("source_stage: ", snap$source_stage),
    paste0("project_mode: ", snap$project_mode),
    paste0("selection_scope: ", plan$selection$scope),
    paste0("sample_mode: ", plan$selection$sample_mode),
    paste0("selected_samples: ", paste(plan$selection$selected_samples, collapse = ",")),
    paste0("isolates: ", paste(plan$selection$isolates, collapse = ",")),
    paste0("loci: ", paste(plan$selection$loci, collapse = ",")),
    paste0("assays: ", paste(plan$selection$assays, collapse = ",")),
    paste0("components_included: ", paste(plan$included_components, collapse = ",")),
    paste0("column_mode: ", plan$selection$column_mode),
    paste0("columns: ", paste(plan$selection$column_ids, collapse = ",")),
    paste0("row_mode: ", plan$selection$row_mode),
    paste0("excel_layout: ", plan$selection$excel_layout),
    paste0("sequence_layout: ", plan$selection$sequence_layout),
    paste0("sequence_header: ", plan$selection$sequence_header),
    ""
  )
  if (length(plan$omitted_components)) {
    lines <- c(lines, "omitted_or_blocked_components:")
    for (nm in names(plan$omitted_components)) {
      om <- plan$omitted_components[[nm]]
      lines <- c(lines, paste0("  - ", nm, " [", om$status, "] ", om$reason))
    }
    lines <- c(lines, "")
  }
  if (length(plan$warnings)) {
    lines <- c(lines, "warnings:")
    lines <- c(lines, paste0("  - ", plan$warnings))
    lines <- c(lines, "")
  }
  if (length(failures)) {
    lines <- c(lines, "component_or_file_failures:")
    for (nm in names(failures)) lines <- c(lines, paste0("  - ", nm, ": ", failures[[nm]]))
    lines <- c(lines, "")
  }
  # Active processing settings fingerprint
  lines <- c(lines, "processing_settings_by_source:")
  for (sid in names(snap$results)) {
    ps <- snap$results[[sid]]$processing_settings
    if (!is.list(ps)) next
    lines <- c(lines, paste0("  - ", sid, ": target=", as.character(ps$target %||% ""),
                             "; assay_id=", as.character(ps$assay_id %||% ""),
                             "; window=", as.character(ps$window %||% ""),
                             "; min_peak_ratio=", as.character(ps$min_peak_ratio %||% "")))
  }
  lines <- c(lines, "", "blast_rids:")
  if (is.data.frame(snap$blast_jobs) && nrow(snap$blast_jobs)) {
    for (i in seq_len(nrow(snap$blast_jobs))) {
      lines <- c(lines, paste0(
        "  - ", snap$blast_jobs$original_name[i],
        " rid=", snap$blast_jobs$rid[i],
        " status=", snap$blast_jobs$status[i],
        " consensus_revision=", snap$blast_jobs$consensus_revision[i]
      ))
    }
  } else {
    lines <- c(lines, "  - none")
  }
  lines <- c(lines, "", "written_files:")
  for (wf in written_files) {
    lines <- c(lines, paste0(
      "  - ", wf$path,
      " component=", wf$component_id,
      " records=", paste(wf$record_ids, collapse = ",")
    ))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  invisible(path)
}
