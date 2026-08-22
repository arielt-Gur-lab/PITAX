  # ---------------- Trimming ----------------
  observeEvent(input$run_trimming, {
    req(input$ab1_files)
    rv$read_assignments <- collect_assignment_editor()
    name_error <- stage2_identity_error(rv$read_assignments)
    if (!is.null(name_error)) {
      showNotification(paste("Rename error:", name_error), type = "error", duration = 10)
      return()
    }
    sync_assignment_state()
    assignment_error <- stage2_validate_assignments(rv$read_assignments, current_upload_source_ids(), rv$assay_profiles)
    if (!is.null(assignment_error)) {
      showNotification(paste("Read assignment error:", assignment_error), type = "error", duration = 10)
      return()
    }
    locus_error <- stage3_run_locus_error(rv$read_assignments)
    if (!is.null(locus_error)) {
      showNotification(locus_error, type = "error", duration = 10)
      return()
    }
    rv$read_assignments <- stage2_coerce_assignments(rv$read_assignments)
    rv$architecture <- tryCatch(
      stage2_build_architecture(rv$read_assignments, assay_profiles = rv$assay_profiles),
      error = function(e) {
        showNotification(paste("Could not build project architecture:", conditionMessage(e)), type = "error", duration = 10)
        NULL
      }
    )
    if (is.null(rv$architecture)) return()
    settings <- current_settings_from_inputs()
    rv$project_defaults <- assay_project_defaults_from_legacy_settings(settings)
    rv$settings <- settings
    all_results <- list(); summaries <- list(); files <- input$ab1_files

    withProgress(message="Processing AB1 files", value=0, {
      for (i in seq_len(nrow(files))) {
        sample_id <- sub("\\.ab1$", "", files$name[i], ignore.case=TRUE)
        assignment_idx <- match(sample_id, rv$read_assignments$Source_ID)
        assignment <- rv$read_assignments[assignment_idx, , drop = FALSE]
        profile_idx <- match(assignment$Assay_ID[1], rv$assay_profiles$Assay_ID)
        read_settings <- assay_resolve_read_settings(
          rv$assay_profiles[profile_idx, , drop = FALSE],
          rv$project_defaults,
          assignment$Direction[1]
        )
        incProgress(1/nrow(files), detail=paste("Processing", files$name[i]))
        result <- tryCatch(
          trim_one_ab1(files$datapath[i], sample_id, read_settings),
          error=function(e) structure(list(error=conditionMessage(e)), class="ab1_error")
        )
        if (inherits(result,"ab1_error")) {
          summaries[[sample_id]] <- make_failure_summary(sample_id, read_settings, result$error)
        } else {
          result$read_assignment <- as.list(assignment[1, , drop = FALSE])
          result$processing_settings <- read_settings
          result <- ensure_curation_state(result)
          result <- curation_rebuild(result, read_settings)
          all_results[[sample_id]] <- result
          summaries[[sample_id]] <- result$summary
        }
      }
    })

    rv$results <- all_results
    rv$summary <- Reduce(rbind_fill, summaries); rownames(rv$summary) <- NULL
    rv$consensus_set <- stage3_empty_consensus_set()
    sync_qc_sample_choices()
    session$sendCustomMessage("showLoader", list(text = "Opening Trim & QC workspace…"))
    updateTabsetPanel(session,"pipeline_step",selected="qc")
  })

  qc_display_name <- function(original_name) {
    original_name <- as.character(original_name)[1]
    if (is.null(rv$rename) || !nrow(rv$rename)) return(original_name)
    idx <- match(original_name, rv$rename$Original_name)
    if (is.na(idx) || !nzchar(trimws(rv$rename$New_name[idx]))) original_name else as.character(rv$rename$New_name[idx])
  }

  sync_qc_sample_choices <- function(preferred = NULL) {
    keys <- names(rv$results)
    if (!length(keys)) {
      updateSelectInput(session, "inspect_sample", choices = character(), selected = character())
      return(invisible(NULL))
    }
    labels <- vapply(keys, qc_display_name, character(1))
    choices <- stats::setNames(keys, labels)
    current <- if (!is.null(preferred) && preferred %in% keys) preferred else isolate(input$inspect_sample)
    selected <- if (!is.null(current) && length(current) == 1L && current %in% keys) current else if (length(keys)) keys[1] else character()
    updateSelectInput(session, "inspect_sample", choices = choices, selected = selected)
  }
