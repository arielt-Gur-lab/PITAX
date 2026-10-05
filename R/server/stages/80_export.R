# ============================================================
# Shared export dialog + retained analysis-record helpers
# ============================================================

  # Analysis-record helpers remain available to BLAST/taxonomy screens.
  read_export_records <- reactive({
    req(rv$results, rv$rename)
    records <- list()
    for (original_name in names(rv$results)) {
      r <- rv$results[[original_name]]
      idx <- match(original_name, rv$rename$Original_name)
      r$final_name <- if (!is.na(idx)) rv$rename$New_name[idx] else original_name
      records[[original_name]] <- r
    }
    records
  })
  read_export_summary_df <- reactive({
    req(rv$summary, rv$rename)
    df <- rv$summary
    df$final_name <- rv$rename$New_name[match(df$sample_id, rv$rename$Original_name)]
    df
  })
  export_records <- reactive({
    req(is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
    stage3_analysis_records(rv$consensus_set)
  })
  export_summary_df <- reactive({
    req(is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
    stage3_analysis_summary(rv$consensus_set)
  })

  project_export_stem <- function() {
    if (is.list(rv$architecture) && is.data.frame(rv$architecture$loci) && "Locus" %in% names(rv$architecture$loci)) {
      locus_names <- unique(trimws(as.character(rv$architecture$loci$Locus)))
      locus_names <- locus_names[nzchar(locus_names)]
      if (length(locus_names) > 1L) return("PITAX_multi_locus")
      if (length(locus_names) == 1L) return(clean_fasta_name(locus_names))
    }
    if (!is.null(rv$settings)) return(clean_fasta_name(stage2_scalar_text(rv$settings$target, "PITAX_project")))
    "PITAX_project"
  }

  observeEvent(input$to_export, {
    e <- stage3_consensus_gate_error(rv$consensus_set, rv$results)
    if (!is.null(e)) {
      showNotification(e, type = "error", duration = 10)
      return()
    }
    workflow_mark_unlocked(also_export = TRUE)
    open_shared_export_dialog("export", pitax_export_component_ids(), "project")
  })

  output$export_summary <- renderUI({
    if (is.null(rv$summary) || is.null(rv$settings)) {
      return(div(class = "compact-hint", "Open Export from any stage. Available components depend on project progress."))
    }
    gate <- tryCatch(stage3_consensus_gate_error(rv$consensus_set, rv$results), error = function(e) conditionMessage(e))
    if (!is.null(gate)) {
      return(tagList(
        p("Analysis sequences are not ready for a full results package."),
        p(gate),
        actionButton("open_export_from_stage", "Open shared Export", icon = icon("download"), class = "btn-primary")
      ))
    }
    consensus_summary <- rv$consensus_set$summary
    loci <- unique(as.character(consensus_summary$Locus))
    architecture_summary <- stage2_architecture_summary(rv$architecture)
    tagList(
      p(strong("Gene / locus: "), paste(loci, collapse = ", ")),
      p(strong("Source architecture: "), paste0(architecture_summary$Isolates, " isolate(s), ", architecture_summary$Reads, " read(s)")),
      p(strong("Source reads processed: "), nrow(rv$summary)),
      p(strong("Isolate-level sequences: "), sum(consensus_summary$Length > 0, na.rm = TRUE)),
      actionButton("export_preset_export_package_summary", "Open shared Export package builder", icon = icon("download"), class = "btn-primary")
    )
  })

  export_dialog_state <- reactiveValues(
    snap = NULL,
    plan = NULL,
    stage_id = "export",
    last_result_message = ""
  )

  pitax_collect_export_selection_from_input <- function(stage_id) {
    scope <- if (!is.null(input$export_scope_mode)) as.character(input$export_scope_mode)[1] else {
      if (identical(stage_id, "export")) "project" else "stage"
    }
    ids <- character()
    for (cid in pitax_export_component_ids()) {
      val <- input[[paste0("export_comp_", cid)]]
      if (isTRUE(val)) ids <- c(ids, cid)
    }
    if (!length(ids)) ids <- pitax_export_default_component_ids(stage_id, scope)
    sel <- pitax_export_default_selection(
      stage_id = stage_id,
      scope = scope,
      sample_mode = if (!is.null(input$export_sample_mode)) as.character(input$export_sample_mode)[1] else "all",
      selected_samples = if (!is.null(input$export_selected_samples)) as.character(input$export_selected_samples) else character(),
      isolates = if (!is.null(input$export_filter_isolates)) as.character(input$export_filter_isolates) else character(),
      loci = if (!is.null(input$export_filter_loci)) as.character(input$export_filter_loci) else character(),
      assays = if (!is.null(input$export_filter_assays)) as.character(input$export_filter_assays) else character()
    )
    sel$component_ids <- ids
    sel$excel_layout <- if (!is.null(input$export_excel_layout)) as.character(input$export_excel_layout)[1] else "sheet_per_component"
    sel$sequence_layout <- if (!is.null(input$export_sequence_layout)) as.character(input$export_sequence_layout)[1] else "combined"
    sel$sequence_header <- if (!is.null(input$export_sequence_header)) as.character(input$export_sequence_header)[1] else "short"
    sel$include_plots <- if (is.null(input$export_include_plots)) TRUE else isTRUE(input$export_include_plots)
    sel
  }

  open_shared_export_dialog <- function(stage_id = NULL, preset_components = NULL, scope = NULL) {
    if (is.null(stage_id) || !nzchar(as.character(stage_id)[1])) {
      stage_id <- if (!is.null(input$pipeline_step)) as.character(input$pipeline_step)[1] else "export"
    }
    stage_id <- as.character(stage_id)[1]
    if (identical(stage_id, "help")) stage_id <- "export"
    snap <- pitax_export_snapshot_from_rv(rv, input, source_stage = stage_id)
    scope_mode <- if (!is.null(scope)) scope else if (identical(stage_id, "export")) "project" else "stage"
    sel <- pitax_export_default_selection(stage_id, scope_mode)
    if (!is.null(preset_components) && length(preset_components)) {
      sel$component_ids <- unique(as.character(preset_components))
      if (!is.null(scope)) sel$scope <- scope
    }
    plan <- pitax_export_build_plan(snap, sel)
    export_dialog_state$snap <- snap
    export_dialog_state$plan <- plan
    export_dialog_state$stage_id <- stage_id

    asg <- snap$read_assignments
    sample_choices <- if (is.data.frame(asg) && nrow(asg)) as.character(asg$Source_ID) else names(snap$results)
    isolate_choices <- if (is.data.frame(asg) && nrow(asg)) unique(as.character(asg$Isolate)) else character()
    locus_choices <- if (is.data.frame(asg) && nrow(asg)) unique(as.character(asg$Locus)) else character()
    assay_choices <- if (is.data.frame(asg) && nrow(asg) && "Assay_ID" %in% names(asg)) unique(as.character(asg$Assay_ID)) else character()
    isolate_choices <- isolate_choices[nzchar(isolate_choices)]
    locus_choices <- locus_choices[nzchar(locus_choices)]
    assay_choices <- assay_choices[nzchar(assay_choices)]

    showModal(pitax_export_dialog_ui(
      stage_id = stage_id,
      availability = plan$availability,
      default_ids = sel$component_ids,
      preview = plan$preview,
      file_plan = plan$files,
      warnings = plan$warnings,
      sample_choices = sample_choices,
      isolate_choices = isolate_choices,
      locus_choices = locus_choices,
      assay_choices = assay_choices,
      current_sample = snap$current_sample
    ))
  }

  observeEvent(input$open_export_output, {
    open_shared_export_dialog(if (!is.null(input$pipeline_step)) as.character(input$pipeline_step)[1] else "export")
  }, ignoreInit = TRUE)

  observeEvent(input$open_export_from_stage, {
    open_shared_export_dialog(if (!is.null(input$pipeline_step)) as.character(input$pipeline_step)[1] else "export")
  }, ignoreInit = TRUE)

  observeEvent(input$export_preset_assign_checkpoint, {
    open_shared_export_dialog("rename", c("assign_identity_table", "assign_rename_map", "assign_architecture"), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_qc_checkpoint, {
    open_shared_export_dialog("qc", c(
      "qc_summary", "qc_peak_flags", "qc_active_curation", "qc_trimmed_sequences",
      "qc_metrics_plots", "assign_identity_table", "assign_rename_map"
    ), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_qc_ab1_run, {
    open_shared_export_dialog("qc", "qc_ab1_run_evidence", "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_qc_ab1_detail, {
    open_shared_export_dialog("qc", "qc_ab1_base_evidence", "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_consensus, {
    open_shared_export_dialog("consensus", c(
      "consensus_summary", "consensus_analysis_sequences", "consensus_column_evidence",
      "consensus_pairwise_alignments", "consensus_conflicts", "consensus_active_curation"
    ), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_consensus_fasta, {
    open_shared_export_dialog("consensus", "consensus_analysis_sequences", "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_export_package, {
    open_shared_export_dialog("export", pitax_export_component_ids(), "project")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_export_package_summary, {
    open_shared_export_dialog("export", pitax_export_component_ids(), "project")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_blast_jobs, {
    open_shared_export_dialog("blast", c("blast_jobs", "blast_hits"), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_blast_selected, {
    open_shared_export_dialog("blast", "blast_query_sequences", "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_taxonomy, {
    open_shared_export_dialog("taxonomy", c(
      "taxonomy_team_summary", "taxonomy_summary", "taxonomy_hits", "taxonomy_species_evidence"
    ), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_analysis_fasta_tile, {
    open_shared_export_dialog("export", "consensus_analysis_sequences", "project")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_taxonomy_checkpoint, {
    open_shared_export_dialog("taxonomy", c(
      "taxonomy_team_summary", "taxonomy_summary", "taxonomy_hits", "taxonomy_species_evidence", "taxonomy_score_plot"
    ), "stage")
  }, ignoreInit = TRUE)
  observeEvent(input$export_preset_multilocus, {
    open_shared_export_dialog("multilocus", c(
      "multilocus_profiles", "multilocus_evidence", "multilocus_sources", "multilocus_sequences"
    ), "stage")
  }, ignoreInit = TRUE)

  output$download_export_package <- downloadHandler(
    filename = function() {
      snap <- export_dialog_state$snap
      if (is.null(snap)) return("PITAX_export.zip")
      sel <- tryCatch(
        pitax_collect_export_selection_from_input(export_dialog_state$stage_id),
        error = function(e) list(stage_id = export_dialog_state$stage_id)
      )
      pitax_export_package_filename(snap, sel)
    },
    content = function(file) {
      snap <- export_dialog_state$snap
      if (is.null(snap)) stop("Export snapshot is missing. Re-open the Export dialog.")
      sel <- pitax_collect_export_selection_from_input(export_dialog_state$stage_id)
      plan <- pitax_export_build_plan(snap, sel)
      export_dialog_state$plan <- plan
      result <- pitax_export_execute(snap, plan, file)
      if (length(result$failures)) {
        export_dialog_state$last_result_message <- paste(
          "Export finished with optional plot failures:",
          paste(names(result$failures), collapse = ", ")
        )
      } else {
        export_dialog_state$last_result_message <- "Export package created."
      }
    }
  )
