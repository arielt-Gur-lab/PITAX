  # ---------------- Upload ----------------
  assignment_state_signature <- function(assignments) {
    assignments <- stage2_coerce_assignments(assignments)
    if (!nrow(assignments)) return("")
    cols <- c("Source_ID", "Isolate", "Locus", "Direction", "Final_Name")
    paste(apply(assignments[, cols, drop = FALSE], 1, paste, collapse = "\r"), collapse = "\n")
  }

  output$uploaded_files_table <- renderDT({
    req(input$ab1_files)
    datatable(data.frame(File=input$ab1_files$name, Size_KB=round(input$ab1_files$size/1024,1)),
              rownames=FALSE, options=list(pageLength=15, dom="tip"))
  })

  observeEvent(input$ab1_files, {
    fresh <- initialize_current_read_assignments()
    previous <- stage2_coerce_assignments(rv$read_assignments)
    if (nrow(previous) && nrow(fresh)) {
      carry <- intersect(c("Isolate", "Locus", "Direction", "Primer", "Notes"), names(previous))
      for (i in seq_len(nrow(fresh))) {
        j <- match(fresh$Source_ID[i], previous$Source_ID)
        if (!is.na(j)) fresh[i, carry] <- previous[j, carry]
      }
    }
    rv$read_assignments <- fresh
    sync_assignment_state()
    rv$project_migration_log <- ""
  })

  sync_project_mode_navigation <- function(mode) {
    simple <- identical(as.character(mode)[1], "simple")
    updateActionButton(session, "to_consensus", label = if (simple) "Continue to Export" else "Continue to Consensus")
    updateActionButton(session, "back_consensus_from_export", label = if (simple) "Back to QC" else "Back to Consensus")
    if (simple) {
      hideTab(inputId = "pipeline_step", target = "consensus", session = session)
    } else {
      showTab(inputId = "pipeline_step", target = "consensus", session = session)
    }
    invisible(NULL)
  }

  session$onFlushed(function() sync_project_mode_navigation(isolate(rv$project_mode)), once = TRUE)

  observeEvent(input$project_mode, {
    mode <- as.character(input$project_mode)[1]
    if (!mode %in% c("simple", "paired_consensus") || identical(mode, rv$project_mode)) return()
    rv$project_mode <- mode
    if (is.list(rv$consensus_set) && length(rv$consensus_set$records)) {
      old_ids <- names(rv$consensus_set$records)
      rv$consensus_set <- stage3_empty_consensus_set()
      for (consensus_id in old_ids) invalidate_downstream_for_sample(consensus_id, "Project read model changed")
    }
    rv$project_status_text <- paste0("Unsaved project read model: ", if (mode == "simple") "simple independent reads" else "Forward/Reverse pairing", ".")
    sync_project_mode_navigation(mode)
  }, ignoreInit = TRUE)

  sync_assignment_state <- function() {
    if (!nrow(rv$read_assignments)) {
      identity_changed <- nzchar(rv$assignment_signature)
      rv$rename <- stage2_default_rename_map(rv$read_assignments)
      rv$architecture <- NULL
      rv$assignment_signature <- ""
      if (identity_changed && is.list(rv$consensus_set) && length(rv$consensus_set$records)) rv$consensus_set <- stage3_empty_consensus_set()
      return(invisible(NULL))
    }
    rv$read_assignments <- stage2_sync_generated_names(
      rv$read_assignments,
      forward_primer = input$forward_primer,
      reverse_primer = input$reverse_primer
    )
    new_signature <- assignment_state_signature(rv$read_assignments)
    identity_changed <- nzchar(rv$assignment_signature) && !identical(new_signature, rv$assignment_signature)
    if (identity_changed && is.list(rv$consensus_set) && length(rv$consensus_set$records)) {
      stale_consensus_ids <- names(rv$consensus_set$records)
      rv$consensus_set <- stage3_empty_consensus_set()
      for (consensus_id in stale_consensus_ids) invalidate_downstream_for_sample(consensus_id, "Read identity changed")
    }
    rv$assignment_signature <- new_signature
    rv$rename <- stage2_default_rename_map(rv$read_assignments)
    assignment_error <- stage2_identity_error(rv$read_assignments)
    rv$architecture <- if (is.null(assignment_error)) stage2_build_architecture(rv$read_assignments) else NULL
    if (is.null(assignment_error) && length(rv$results)) {
      for (source_id in intersect(names(rv$results), rv$read_assignments$Source_ID)) {
        i <- match(source_id, rv$read_assignments$Source_ID)
        rv$results[[source_id]]$read_assignment <- as.list(rv$read_assignments[i, , drop = FALSE])
        if (is.list(rv$results[[source_id]]$processing_settings)) {
          rv$results[[source_id]]$processing_settings$target <- rv$read_assignments$Locus[i]
          rv$results[[source_id]]$processing_settings$sequencing_primer <- rv$read_assignments$Direction[i]
          if (rv$read_assignments$Direction[i] == "Forward" && nzchar(rv$read_assignments$Primer[i])) rv$results[[source_id]]$processing_settings$forward_primer <- rv$read_assignments$Primer[i]
          if (rv$read_assignments$Direction[i] == "Reverse" && nzchar(rv$read_assignments$Primer[i])) rv$results[[source_id]]$processing_settings$reverse_primer <- rv$read_assignments$Primer[i]
        }
      }
    }
    invisible(assignment_error)
  }

  observeEvent(list(input$target, input$sequencing_primer, input$forward_primer, input$reverse_primer), {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return()
    sync_assignment_state()
  }, ignoreInit = TRUE)

  architecture_summary_ui <- function() {
    error <- stage2_identity_error(rv$read_assignments)
    if (!is.null(error)) return(div(class = "status-warning", "Architecture preview will appear after every read has explicit Isolate, Gene / locus and Forward / Reverse fields."))
    architecture <- tryCatch(stage2_build_architecture(rv$read_assignments), error = function(e) NULL)
    if (is.null(architecture)) return(NULL)
    sm <- stage2_architecture_summary(architecture)
    div(
      class = "compact-hint",
      paste0(
        "Architecture preview: ", sm$Isolates, " isolate(s) · ", sm$Loci, " locus/loci · ", sm$Reads, " read(s) · ",
        if (identical(rv$project_mode, "simple")) paste0(sm$Reads, " independent analysis sequence(s).") else paste0(sm$Paired_loci, " Forward/Reverse pair(s) · ", sm$Single_read_loci, " single-read locus/loci.")
      )
    )
  }
  output$architecture_summary <- renderUI(architecture_summary_ui())

  observeEvent(input$to_settings, {
    if (is.null(input$ab1_files) || nrow(input$ab1_files)==0) {
      showNotification("Please upload at least one AB1 file.", type="error"); return()
    }
    updateTabsetPanel(session, "pipeline_step", selected="settings")
  })
  observeEvent(input$back_upload, updateTabsetPanel(session,"pipeline_step",selected="upload"))
  observeEvent(input$back_settings_from_rename, updateTabsetPanel(session,"pipeline_step",selected="settings"))
  observeEvent(input$back_rename_from_qc, updateTabsetPanel(session,"pipeline_step",selected="rename"))
  observeEvent(input$back_qc_from_consensus, updateTabsetPanel(session,"pipeline_step",selected="qc"))
  observeEvent(input$back_consensus_from_export, updateTabsetPanel(session,"pipeline_step",selected=if (identical(rv$project_mode, "simple")) "qc" else "consensus"))
  observeEvent(input$back_export, updateTabsetPanel(session,"pipeline_step",selected="export"))
  observeEvent(input$back_blast, updateTabsetPanel(session,"pipeline_step",selected="blast"))

  current_settings_from_inputs <- function() {
    list(
      target=input$target,
      forward_primer=input$forward_primer,
      forward_primer_seq=sanitize_dna(input$forward_primer_seq),
      reverse_primer=input$reverse_primer,
      reverse_primer_seq=sanitize_dna(input$reverse_primer_seq),
      sequencing_primer=input$sequencing_primer,
      enable_primer_mapping=isTRUE(input$enable_primer_mapping),
      expected_amplicon_len=as.integer(input$expected_amplicon_len),
      absolute_max_base_index=as.integer(input$absolute_max_base_index),
      window=as.integer(input$window),
      min_peak_ratio=as.numeric(input$min_peak_ratio),
      min_relative_signal=as.numeric(input$min_relative_signal),
      min_len_before_collapse=as.integer(input$min_len_before_collapse),
      bad_run_windows=as.integer(input$bad_run_windows),
      min_usable_len=as.integer(input$min_usable_len),
      ambiguous_peak_strong_ratio=1.25,
      ambiguous_peak_moderate_ratio=1.75,
      ambiguous_peak_min_relative_signal=0.20,
      auto_correct_min_alt_to_called=1.80,
      auto_correct_min_alt_to_third=2.00,
      auto_correct_max_peak_offset=2L,
      auto_correct_min_relative_signal=0.50
    )
  }

  observeEvent(input$to_rename, {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) {
      showNotification("Please upload at least one AB1 file.", type = "error")
      return()
    }
    source_ids <- current_upload_source_ids()
    if (!nrow(rv$read_assignments) || !identical(sort(rv$read_assignments$Source_ID), sort(source_ids))) {
      rv$read_assignments <- initialize_current_read_assignments()
    }
    rv$settings <- current_settings_from_inputs()
    sync_assignment_state()
    updateTabsetPanel(session, "pipeline_step", selected = "rename")
  })

