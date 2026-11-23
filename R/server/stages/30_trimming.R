  # ---------------- Trimming ----------------
  observeEvent(input$run_trimming, {
    files <- input$ab1_files
    if (is.null(files) || !is.data.frame(files) || !nrow(files)) {
      showNotification(
        if (length(rv$results)) "No newly uploaded AB1 files. Use Re-trim saved reads to re-run trimming without uploading or assigning again." else "Upload at least one AB1 file and complete Assign before starting trimming.",
        type = "warning",
        duration = 8
      )
      return()
    }
    if (!is.data.frame(rv$read_assignments) || !nrow(rv$read_assignments)) {
      showNotification(
        "Assign isolate, assay and direction to every uploaded read before starting trimming.",
        type = "warning",
        duration = 8
      )
      return()
    }
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
    all_results <- list(); summaries <- list()

    withProgress(message = paste0("Trimming ", nrow(files), " AB1 file(s)"), value = 0, {
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
        incProgress(
          1/nrow(files),
          detail = paste0(i, "/", nrow(files), " | ", files$name[i])
        )
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
    workflow_mark_unlocked("upload", "settings", "rename", "qc")
    workflow_mark_completed("upload", "settings", "rename", "qc")
    session$sendCustomMessage("showLoader", list(text = "Opening Trim & QC workspace..."))
    updateTabsetPanel(session,"pipeline_step",selected="qc")
  })

  # Re-trim from immutable raw trace data embedded in a saved project.
  # Unlike Start trimming, this does NOT use input$ab1_files and never calls
  # collect_assignment_editor() or sync_assignment_state() (which can replace
  # saved Final_Name values with defaults after a project load).
  saved_retrim_preview <- reactiveVal(NULL)

  open_saved_retrim_dialog <- function(scope = c("all", "selected")) {
    scope <- match.arg(scope)
    results <- isolate(rv$results)
    if (!is.list(results) || !length(results)) {
      showNotification("Load a project with previously processed chromatograms before re-trimming.", type = "warning", duration = 8)
      return(invisible(FALSE))
    }
    keys <- names(results)
    if (identical(scope, "selected")) {
      selected <- as.character(isolate(input$inspect_sample))[1]
      if (!length(selected) || is.na(selected) || !selected %in% keys) {
        showNotification("Choose a saved read in Sequence inspection.", type = "warning")
        return(invisible(FALSE))
      }
      keys <- selected
    }
    preview <- tryCatch(pitax_retrim_saved_batch(results, keys), error = function(e) e)
    if (inherits(preview, "error")) {
      showNotification(
        paste0("Could not prepare re-trim: ", conditionMessage(preview)),
        type = "error", duration = 14
      )
      return(invisible(FALSE))
    }
    saved_retrim_preview(preview)
    changes <- sum(preview$summary$Sequence_changed)
    showModal(modalDialog(
      title = paste0("Confirm re-trim of ", nrow(preview$summary), " saved read(s)"),
      div(class = "status-note",
          tags$strong(paste0(nrow(preview$summary), " read(s) ready; ", changes, " resulting sequence(s) will change.")),
          tags$p("Re-trim uses the chromatogram and original base calls ALREADY SAVED with this project. No AB1 upload, Assay setup or Assign step is needed. Existing read identities, locus/assay assignments and names are preserved."),
          tags$p("Existing manual base substitutions remain, but old manual TRIM boundaries are replaced by the new automatic bounds. Their previous values are recorded in re-trim history. Existing undo/redo stacks are reset for updated reads; the audit log is kept."),
          if (changes > 0) tags$p(class = "status-warning",
             "Changed sequences will invalidate their BLAST and taxonomy evidence. The affected consensus/analysis sequences will require a rebuild before further BLAST. Unchanged sequences retain their evidence.") else NULL,
          tags$p("Save a copy of the current project first if you want an independent backup. This action cannot be undone using the ordinary curation Undo button.")),
      DTOutput("saved_retrim_preview_table"),
      easyClose = FALSE, size = "l",
      footer = tagList(
        actionButton("cancel_saved_retrim", "Cancel"),
        actionButton("confirm_saved_retrim", "Confirm re-trim", class = "btn-warning", icon = icon("refresh"))
      )
    ))
    invisible(TRUE)
  }

  observeEvent(input$open_saved_retrim, {
    scope <- if (identical(input$saved_retrim_scope, "selected")) "selected" else "all"
    open_saved_retrim_dialog(scope)
  }, ignoreInit = TRUE)

  observeEvent(input$open_saved_retrim_from_assign, {
    open_saved_retrim_dialog("all")
  }, ignoreInit = TRUE)

  output$saved_retrim_preview_table <- renderDT({
    preview <- saved_retrim_preview()
    req(preview)
    display <- preview$summary
    display$Sequence_changed <- ifelse(display$Sequence_changed, "Yes", "No")
    names(display) <- c("Source read", "Old bounds", "New bounds", "Old length",
                        "New length", "Sequence changed")
    datatable(display, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
              selection = "none")
  })

  observeEvent(input$cancel_saved_retrim, {
    saved_retrim_preview(NULL)
    removeModal()
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_saved_retrim, {
    preview <- saved_retrim_preview()
    if (is.null(preview)) return()
    current <- isolate(rv$results)
    ids <- names(preview$results)
    unchanged <- all(vapply(ids, function(id) {
      is.list(current[[id]]) && identical(
        pitax_saved_retrim_before(current[[id]]), preview$original[[id]]
      )
    }, logical(1)))
    if (!unchanged) {
      saved_retrim_preview(NULL)
      removeModal()
      showNotification("The saved read state changed after preview. Re-open the re-trim preview to avoid overwriting newer edits.", type = "error", duration = 12)
      return()
    }

    # Commit all validated candidate reads together. Validation happens
    # before any mutation, so legacy/missing-trace reads cannot cause a
    # half-retrimmed project.
    for (id in ids) current[[id]] <- preview$results[[id]]
    changed_sources <- preview$changed_sources
    affected_analysis <- character()
    if (length(changed_sources) && is.list(rv$consensus_set) &&
        is.list(rv$consensus_set$records) && length(rv$consensus_set$records)) {
      affected_analysis <- names(rv$consensus_set$records)[vapply(
        rv$consensus_set$records, function(record) {
          any(as.character(record$source_read_ids) %in% changed_sources)
        }, logical(1)
      )]
    }
    rv$results <- current
    sync_summary_from_results()
    # Preserve the existing consensus history but prevent the stale analysis
    # from being used. stage3_consensus_gate_error detects changed sources.
    # BLAST/taxonomy are invalidated only for affected analysis IDs.
    for (key in unique(c(changed_sources, affected_analysis))) {
      invalidate_downstream_for_sample(key, "Saved chromatogram re-trim changed the sequence")
    }
    saved_retrim_preview(NULL)
    removeModal()
    sync_qc_sample_choices()
    workflow_mark_unlocked("qc")
    workflow_mark_completed("qc")
    rv$project_status_text <- paste0(
      "Unsaved re-trim from stored chromatograms: ", length(ids), " processed; ",
      length(changed_sources), " sequence(s) changed. Assay and Assign were preserved."
    )
    showNotification(
      paste0("Re-trim complete: ", length(ids), " saved reads processed, ",
             length(changed_sources), " sequence(s) changed. ",
             if (length(changed_sources)) "Rebuild affected analysis sequences before BLAST. " else "",
             "Save the updated project."),
      type = "message", duration = 12
    )
    updateTabsetPanel(session, "pipeline_step", selected = "qc")
  }, ignoreInit = TRUE)

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
