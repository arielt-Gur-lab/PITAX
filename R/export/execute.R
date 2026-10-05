# ============================================================
# Export package execution
# ============================================================

pitax_export_execute <- function(snap, plan, dest_zip) {
  validated <- pitax_export_validate_plan(snap, plan)
  if (!isTRUE(validated$ok)) {
    stop(validated$error, call. = FALSE)
  }
  plan <- validated$plan

  temp_dir <- tempfile("pitax_export_")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  written <- list()
  failures <- list()
  scope_ids <- pitax_export_filter_scope_ids(snap, plan$selection)

  # Excel workbooks
  excel_files <- Filter(function(f) identical(f$kind, "excel"), plan$files)
  for (ef in excel_files) {
    out_path <- file.path(temp_dir, ef$path)
    iso <- if (identical(plan$selection$excel_layout, "workbook_per_isolate") && length(ef$record_ids)) {
      ef$record_ids[1]
    } else NULL
    tryCatch({
      pitax_export_write_excel_workbook(out_path, snap, plan, isolate_filter = iso)
      written[[length(written) + 1L]] <- ef
    }, error = function(e) {
      stop("Required Excel export failed for ", ef$path, ": ", conditionMessage(e), call. = FALSE)
    })
  }

  # Sequences
  for (id in plan$sequence_components) {
    extractor <- pitax_export_extractor(id)
    seqs <- extractor(snap, scope_ids)
    related <- Filter(function(f) identical(f$component_id, id) && f$kind %in% c("fasta", "text"), plan$files)
    for (ff in related) {
      out_path <- file.path(temp_dir, ff$path)
      if (grepl("_missing_sequences\\.txt$", ff$path)) {
        lines <- vapply(ff$record_ids, function(rid) {
          rec <- seqs[[rid]]
          paste0(rid, "\t", rec$final_name %||% "", "\tmissing_sequence\trole=", rec$sequence_role %||% "")
        }, character(1))
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(lines, out_path)
        written[[length(written) + 1L]] <- ff
        next
      }
      subset_recs <- seqs[intersect(names(seqs), ff$record_ids)]
      if (!length(subset_recs) && length(ff$record_ids) == length(seqs)) subset_recs <- seqs
      # per-record / per-group plans already encode record ids
      if (!length(subset_recs) && length(ff$record_ids) == 1L && ff$record_ids %in% names(seqs)) {
        subset_recs <- seqs[ff$record_ids]
      }
      if (!length(subset_recs) && identical(plan$selection$sequence_layout, "combined")) subset_recs <- seqs
      ans <- tryCatch(
        pitax_export_write_fasta_records(out_path, subset_recs, plan$selection$sequence_header %||% "short"),
        error = function(e) {
          stop("Required FASTA export failed for ", ff$path, ": ", conditionMessage(e), call. = FALSE)
        }
      )
      written[[length(written) + 1L]] <- ff
      if (length(ans$missing_lines)) {
        miss_path <- file.path(temp_dir, paste0("sequences/", id, "_missing_sequences.txt"))
        dir.create(dirname(miss_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(ans$missing_lines, miss_path)
      }
    }
  }

  # Text artifacts
  for (id in plan$text_components) {
    extractor <- pitax_export_extractor(id)
    items <- extractor(snap, scope_ids)
    related <- Filter(function(f) identical(f$component_id, id) && identical(f$kind, "text"), plan$files)
    for (ff in related) {
      item <- items[[ff$record_ids[1]]]
      if (is.null(item)) stop("Required text export missing for ", ff$path, call. = FALSE)
      out_path <- file.path(temp_dir, ff$path)
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
      writeLines(as.character(item$text %||% ""), out_path)
      written[[length(written) + 1L]] <- ff
    }
  }

  # Raw AB1 copies
  for (id in plan$raw_components) {
    extractor <- pitax_export_extractor(id)
    df <- extractor(snap, scope_ids)
    related <- Filter(function(f) identical(f$component_id, id), plan$files)
    for (ff in related) {
      row <- df[df$Source_ID == ff$record_ids[1], , drop = FALSE]
      if (!nrow(row)) stop("Required raw AB1 file missing for ", ff$path, call. = FALSE)
      src <- as.character(row$datapath[1])
      if (!file.exists(src)) stop("Original AB1 bytes unavailable for ", ff$path, call. = FALSE)
      out_path <- file.path(temp_dir, ff$path)
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(src, out_path, overwrite = TRUE)
      if (!isTRUE(ok)) stop("Failed to copy raw AB1 for ", ff$path, call. = FALSE)
      written[[length(written) + 1L]] <- ff
    }
  }

  # Plots (optional failures)
  if (isTRUE(plan$selection$include_plots)) {
    for (id in plan$plot_components) {
      related <- Filter(function(f) identical(f$component_id, id) && identical(f$kind, "plot"), plan$files)
      for (ff in related) {
        out_path <- file.path(temp_dir, ff$path)
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
        tg <- ff$record_ids[1]
        err <- tryCatch({
          if (identical(id, "qc_metrics_plots")) {
            result <- snap$results[[tg]]
            if (is.null(result)) stop("Missing result for ", tg)
            settings <- pitax_export_settings_for_result(result, snap$settings)
            write_qc_plot_png(result, settings, out_path)
          } else if (identical(id, "qc_chromatograms")) {
            result <- snap$results[[tg]]
            if (is.null(result)) stop("Missing result for ", tg)
            settings <- pitax_export_settings_for_result(result, snap$settings)
            pitax_export_write_chromatogram_png(result, settings, out_path)
          } else if (identical(id, "taxonomy_score_plot")) {
            hits <- snap$taxonomy_hits
            if (is.data.frame(hits) && "original_name" %in% names(hits)) {
              hits <- hits[hits$original_name == tg, , drop = FALSE]
            }
            pitax_export_write_taxonomy_score_png(hits, out_path, tg)
          } else if (identical(id, "multilocus_plot")) {
            ev <- pitax_export_extract_multilocus_evidence(snap, scope_ids)
            if (is.data.frame(ev) && "Isolate" %in% names(ev)) ev <- ev[ev$Isolate == tg, , drop = FALSE]
            pitax_export_write_multilocus_png(ev, out_path, tg)
          } else {
            stop("Unknown plot component")
          }
          written[[length(written) + 1L]] <- ff
          NULL
        }, error = function(e) conditionMessage(e))
        if (!is.null(err)) {
          failures[[ff$path]] <- err
        }
      }
    }
  }

  run_info_path <- file.path(temp_dir, "RUN_INFO.txt")
  pitax_export_write_run_info(run_info_path, snap, plan, written, failures)
  written[[length(written) + 1L]] <- pitax_export_planned_file("RUN_INFO.txt", "run_info", kind = "text")

  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("Export produced no files.", call. = FALSE)
  zip::zipr(dest_zip, files, root = temp_dir)

  invisible(list(
    zip = dest_zip,
    written = written,
    failures = failures,
    warnings = c(plan$warnings, if (length(failures)) paste0("Optional plot failures: ", length(failures)) else character())
  ))
}

pitax_export_package_filename <- function(snap, selection = list()) {
  stem <- "PITAX_export"
  if (is.list(snap$architecture) && is.data.frame(snap$architecture$loci) && "Locus" %in% names(snap$architecture$loci)) {
    locus_names <- unique(trimws(as.character(snap$architecture$loci$Locus)))
    locus_names <- locus_names[nzchar(locus_names)]
    if (length(locus_names) == 1L) stem <- clean_fasta_name(locus_names)
    if (length(locus_names) > 1L) stem <- "PITAX_multi_locus"
  } else if (is.list(snap$settings) && !is.null(snap$settings$target)) {
    stem <- clean_fasta_name(as.character(snap$settings$target)[1])
  }
  stage <- as.character(selection$stage_id %||% snap$source_stage %||% "export")[1]
  paste0(stem, "_export_", stage, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
}
