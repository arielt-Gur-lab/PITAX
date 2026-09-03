  # ---------------- Upload ----------------
  assignment_state_signature <- function(assignments) {
    assignments <- stage2_coerce_assignments(assignments)
    if (!nrow(assignments)) return("")
    cols <- c("Source_ID", "Isolate", "Assay_ID", "Locus", "Direction", "Final_Name")
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
    updateActionButton(session, "to_consensus", label = if (simple) "Continue to NCBI BLAST" else "Continue to Consensus")
    updateActionButton(session, "back_from_export", label = if (simple) "Back to Trim & QC" else "Back to Consensus")
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
    if (identical(mode, "simple")) {
      rv$workflow_unlocked <- setdiff(as.character(rv$workflow_unlocked), "consensus")
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
      reverse_primer = input$reverse_primer,
      assay_profiles = rv$assay_profiles
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
    if (is.null(assignment_error)) assignment_error <- stage2_validate_assignments(rv$read_assignments, assay_profiles = rv$assay_profiles)
    rv$architecture <- if (is.null(assignment_error)) stage2_build_architecture(rv$read_assignments, assay_profiles = rv$assay_profiles) else NULL
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

  selected_assay_id <- reactiveVal("")

  sync_assay_editor_choices <- function(preferred = NULL) {
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    if (!nrow(profiles)) {
      profiles <- assay_default_profiles()
      rv$assay_profiles <- profiles
    }
    labels <- paste0(profiles$Assay_Name, " | ", profiles$Locus_Display_Name)
    choices <- setNames(profiles$Assay_ID, labels)
    preferred <- if (!is.null(preferred) && preferred %in% profiles$Assay_ID) preferred else if (nzchar(selected_assay_id()) && selected_assay_id() %in% profiles$Assay_ID) selected_assay_id() else profiles$Assay_ID[1]
    selected_assay_id(preferred)
    updateSelectInput(session, "assay_editor_select", choices = choices, selected = preferred)
    invisible(preferred)
  }

  load_assay_editor_inputs <- function(assay_id) {
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    idx <- match(assay_id, profiles$Assay_ID)
    if (is.na(idx)) return(invisible(NULL))
    row <- profiles[idx, , drop = FALSE]
    updateTextInput(session, "assay_name", value = row$Assay_Name[1])
    updateSelectInput(session, "target", selected = row$Locus_ID[1])
    updateTextInput(session, "forward_primer", value = row$Forward_Primer_Name[1])
    updateTextInput(session, "reverse_primer", value = row$Reverse_Primer_Name[1])
    updateTextInput(session, "forward_primer_seq", value = row$Forward_Primer_Sequence[1])
    updateTextInput(session, "reverse_primer_seq", value = row$Reverse_Primer_Sequence[1])
    updateNumericInput(session, "expected_amplicon_len", value = row$Expected_Amplicon_Length[1])
    updateNumericInput(session, "absolute_max_base_index", value = row$Maximum_Sequence_Position[1])
    invisible(NULL)
  }

  commit_active_assay_from_inputs <- function() {
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    if (!nrow(profiles)) profiles <- assay_default_profiles()
    result <- assay_try_apply_editor_inputs(profiles, selected_assay_id(), list(
      assay_name = input$assay_name,
      target = input$target,
      forward_primer = input$forward_primer,
      reverse_primer = input$reverse_primer,
      forward_primer_seq = input$forward_primer_seq,
      reverse_primer_seq = input$reverse_primer_seq,
      expected_amplicon_len = input$expected_amplicon_len,
      absolute_max_base_index = input$absolute_max_base_index
    ))
    if (!isTRUE(result$committed)) return(invisible(NULL))
    idx <- match(result$assay_id, result$profiles$Assay_ID)
    if (is.na(idx) || !length(idx)) return(invisible(NULL))
    rv$assay_profiles <- result$profiles
    rv$project_defaults <- assay_project_defaults_from_legacy_settings(current_settings_from_inputs())
    rv$settings <- assay_resolve_read_settings(result$profiles[idx, , drop = FALSE], rv$project_defaults, input$sequencing_primer)
    selected_assay_id(result$assay_id)
    invisible(result$assay_id)
  }

  observe({
    rv$assay_profiles
    isolate(sync_assay_editor_choices(selected_assay_id()))
  })

  observeEvent(input$assay_editor_select, {
    req(nzchar(input$assay_editor_select))
    if (nzchar(selected_assay_id()) && !identical(selected_assay_id(), input$assay_editor_select)) {
      commit_active_assay_from_inputs()
    }
    selected_assay_id(input$assay_editor_select)
    load_assay_editor_inputs(input$assay_editor_select)
    sync_assay_editor_choices(input$assay_editor_select)
  }, ignoreInit = TRUE)

  observeEvent(list(input$assay_name, input$target, input$forward_primer, input$reverse_primer,
                    input$forward_primer_seq, input$reverse_primer_seq,
                    input$expected_amplicon_len, input$absolute_max_base_index,
                    input$window, input$min_peak_ratio, input$min_relative_signal,
                    input$min_len_before_collapse, input$bad_run_windows, input$min_usable_len), {
    if (!nzchar(selected_assay_id())) return()
    if (is.null(input$target) || !nzchar(as.character(input$target)[1])) return()
    commit_active_assay_from_inputs()
    sync_assay_editor_choices(selected_assay_id())
    if (!is.null(input$ab1_files) && nrow(input$ab1_files)) sync_assignment_state()
  }, ignoreInit = TRUE)

  observeEvent(input$assay_add_profile, {
    commit_active_assay_from_inputs()
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    new_row <- assay_profile_from_legacy_settings(list(target = "ITS", assay_name = "New assay"), assay_id = assay_make_id("ITS", profiles$Assay_ID))
    new_row$Assay_Name[1] <- "New assay"
    rv$assay_profiles <- rbind(profiles, new_row)
    selected_assay_id(new_row$Assay_ID[1])
    sync_assay_editor_choices(new_row$Assay_ID[1])
    load_assay_editor_inputs(new_row$Assay_ID[1])
  })

  observeEvent(input$assay_remove_profile, {
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    if (nrow(profiles) <= 1L) {
      showNotification("At least one assay profile is required.", type = "warning")
      return()
    }
    assay_id <- selected_assay_id()
    if (any(rv$read_assignments$Assay_ID == assay_id)) {
      showNotification("Cannot remove an assay that is still assigned to reads.", type = "error", duration = 8)
      return()
    }
    profiles <- profiles[profiles$Assay_ID != assay_id, , drop = FALSE]
    rv$assay_profiles <- profiles
    selected_assay_id(profiles$Assay_ID[1])
    sync_assay_editor_choices(profiles$Assay_ID[1])
    load_assay_editor_inputs(profiles$Assay_ID[1])
  })

  session$onFlushed(function() {
    isolate({
      sync_assay_editor_choices()
      load_assay_editor_inputs(selected_assay_id())
    })
  }, once = TRUE)

  architecture_summary_ui <- function() {
    error <- stage2_identity_error(rv$read_assignments)
    if (is.null(error)) error <- stage2_validate_assignments(rv$read_assignments, assay_profiles = rv$assay_profiles)
    if (!is.null(error)) return(div(class = "status-warning", "Architecture preview will appear after every read has explicit Isolate, Gene / locus and Forward / Reverse fields."))
    architecture <- tryCatch(stage2_build_architecture(rv$read_assignments, assay_profiles = rv$assay_profiles), error = function(e) NULL)
    if (is.null(architecture)) return(NULL)
    sm <- stage2_architecture_summary(architecture)
    div(
      class = "compact-hint",
      paste0(
        "Architecture preview: ", sm$Isolates, " isolate(s) | ", sm$Loci, " locus/loci | ", sm$Reads, " read(s) | ",
        if (identical(rv$project_mode, "simple")) paste0(sm$Reads, " independent analysis sequence(s).") else paste0(sm$Paired_loci, " Forward/Reverse pair(s) | ", sm$Single_read_loci, " single-read locus/loci.")
      )
    )
  }
  output$architecture_summary <- renderUI(architecture_summary_ui())

  observeEvent(input$to_settings, {
    workflow_mark_unlocked("upload", "settings")
    updateTabsetPanel(session, "pipeline_step", selected="settings")
  })
  observeEvent(input$back_upload, updateTabsetPanel(session,"pipeline_step",selected="upload"))
  observeEvent(input$back_settings_from_rename, updateTabsetPanel(session,"pipeline_step",selected="settings"))
  observeEvent(input$back_rename_from_qc, updateTabsetPanel(session,"pipeline_step",selected="rename"))
  observeEvent(input$back_qc_from_consensus, updateTabsetPanel(session,"pipeline_step",selected="qc"))
  observeEvent(input$back_from_export, {
    updateTabsetPanel(session,"pipeline_step",selected=if (identical(rv$project_mode, "simple")) "qc" else "consensus")
  })
  observeEvent(input$back_to_process, {
    updateTabsetPanel(session,"pipeline_step",selected=if (identical(rv$project_mode, "simple")) "qc" else "consensus")
  })
  observeEvent(input$back_blast, updateTabsetPanel(session,"pipeline_step",selected="blast"))
  observeEvent(input$open_export_output, {
    workflow_goto("export")
  })

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
    source_ids <- current_upload_source_ids()
    if (length(source_ids) && (!nrow(rv$read_assignments) || !identical(sort(rv$read_assignments$Source_ID), sort(source_ids)))) {
      rv$read_assignments <- initialize_current_read_assignments()
    }
    rv$settings <- current_settings_from_inputs()
    commit_active_assay_from_inputs()
    sync_assignment_state()
    workflow_mark_unlocked("upload", "settings", "rename")
    # Explicit Assay completion - never inferred from default ITS settings on load.
    workflow_mark_completed("settings")
    if (length(source_ids)) workflow_mark_completed("upload")
    updateTabsetPanel(session, "pipeline_step", selected = "rename")
  })
