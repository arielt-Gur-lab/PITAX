# ============================================================
# Free workflow navigation (stepper is orientation, not a gate)
# ============================================================

  # rv$workflow_unlocked is initialized in 00_state.R; kept for progress hints /
  # project load compatibility. Navigation itself is not gated on it.
  # rv$workflow_completed holds explicit "user finished this step" flags (e.g. Assay
  # after Continue to Assign). Default assay settings alone never mark done.

  workflow_all_step_ids <- function(project_mode = isolate(rv$project_mode)) {
    ids <- vapply(pitax_workflow_catalog(project_mode), function(s) s$id, character(1))
    unique(c(ids, "export", "help"))
  }

  workflow_mark_unlocked <- function(..., also_export = FALSE) {
    steps <- unique(c(as.character(unlist(list(...))), "upload"))
    steps <- steps[nzchar(steps)]
    current <- isolate(as.character(rv$workflow_unlocked))
    if (!length(current) || all(!nzchar(current))) current <- "upload"
    rv$workflow_unlocked <- unique(c(current, steps))
    if (isTRUE(also_export)) {
      rv$workflow_unlocked <- unique(c(as.character(isolate(rv$workflow_unlocked)), "export"))
    }
    invisible(isolate(rv$workflow_unlocked))
  }

  workflow_mark_completed <- function(...) {
    steps <- unique(as.character(unlist(list(...))))
    steps <- steps[nzchar(steps)]
    current <- isolate(as.character(rv$workflow_completed))
    if (!length(current) || all(!nzchar(current))) current <- character(0)
    rv$workflow_completed <- unique(c(current, steps))
    invisible(isolate(rv$workflow_completed))
  }

  workflow_can_visit <- function(step) {
    step <- as.character(step)[1]
    if (!nzchar(step)) return(FALSE)
    if (identical(step, "help")) return(TRUE)
    if (identical(step, "export")) return(TRUE)
    if (identical(step, "consensus") && identical(isolate(rv$project_mode), "simple")) return(FALSE)
    step %in% workflow_all_step_ids()
  }

  # Green chip markers: real success only — never default assay settings / visit order.
  workflow_completed_steps <- function() {
    done <- as.character(rv$workflow_completed)
    if (!length(done) || all(!nzchar(done))) done <- character(0)

    has_files <- !is.null(input$ab1_files) && is.data.frame(input$ab1_files) && nrow(input$ab1_files) > 0
    has_assignments <- is.data.frame(rv$read_assignments) && nrow(rv$read_assignments) > 0
    has_results <- is.list(rv$results) && length(rv$results) > 0

    if (has_files) done <- c(done, "upload")
    # Assay is NEVER inferred from rv$settings alone (defaults fill that on load).
    if (has_results) done <- c(done, "settings", "rename", "qc")

    assign_ok <- FALSE
    if (has_files && has_assignments) {
      assign_ok <- is.null(tryCatch(
        stage2_identity_error(rv$read_assignments),
        error = function(e) "assignment identity check failed"
      ))
      if (assign_ok) {
        assign_ok <- is.null(tryCatch(
          stage2_validate_assignments(rv$read_assignments, assay_profiles = rv$assay_profiles),
          error = function(e) "assignment validation failed"
        ))
      }
    }
    if (assign_ok) done <- c(done, "rename")

    consensus_ready <- is.list(rv$consensus_set) &&
      length(rv$consensus_set$records) > 0 &&
      is.null(tryCatch(
        stage3_consensus_gate_error(rv$consensus_set, rv$results),
        error = function(e) "consensus gate failed"
      ))
    if (identical(rv$project_mode, "paired_consensus") && consensus_ready) {
      done <- c(done, "consensus")
    }
    if (consensus_ready) done <- c(done, "export")

    if (is.data.frame(rv$blast_hits) && nrow(rv$blast_hits) > 0) done <- c(done, "blast")
    if (is.data.frame(rv$taxonomy_summary) && nrow(rv$taxonomy_summary) > 0) done <- c(done, "taxonomy")

    ml <- rv$multilocus_profile
    if (is.list(ml) && (
      (is.data.frame(ml$profiles) && nrow(ml$profiles) > 0) ||
      (is.data.frame(ml$evidence) && nrow(ml$evidence) > 0)
    )) {
      done <- c(done, "multilocus")
    }

    unique(done)
  }

  workflow_hide_loader <- function() {
    session$sendCustomMessage("hideLoader", list())
  }

  workflow_goto <- function(step, notify = TRUE) {
    step <- as.character(step)[1]
    if (identical(step, "help")) {
      updateTabsetPanel(session, "pipeline_step", selected = "help")
      workflow_hide_loader()
      return(TRUE)
    }
    if (identical(step, "consensus") && identical(isolate(rv$project_mode), "simple")) {
      if (isTRUE(notify)) {
        showNotification("Consensus is only used in Paired Forward/Reverse mode.", type = "message", duration = 5)
      }
      workflow_hide_loader()
      return(FALSE)
    }
    if (!workflow_can_visit(step)) {
      if (isTRUE(notify)) {
        showNotification("Unknown workflow step.", type = "warning", duration = 5)
      }
      workflow_hide_loader()
      return(FALSE)
    }
    updateTabsetPanel(session, "pipeline_step", selected = step)
    workflow_hide_loader()
    TRUE
  }

  output$workflow_stepper <- renderUI({
    mode <- rv$project_mode
    catalog_ids <- vapply(pitax_workflow_catalog(mode), function(s) s$id, character(1))
    current <- if (!is.null(input$pipeline_step)) input$pipeline_step else "upload"
    # Touch reactive inputs used for completion so the stepper refreshes with real state.
    invisible(input$ab1_files)
    invisible(rv$workflow_completed)
    invisible(rv$read_assignments)
    invisible(rv$results)
    invisible(rv$consensus_set)
    invisible(rv$blast_hits)
    invisible(rv$taxonomy_summary)
    invisible(rv$multilocus_profile)
    pitax_workflow_stepper_ui(
      current = current,
      unlocked = c(catalog_ids, "export", "help"),
      project_mode = mode,
      completed = workflow_completed_steps(),
      free_nav = TRUE
    )
  })

  # Chip navigation via Shiny actionButtons (ids workflow_goto_*). No JS overlay.
  workflow_chip_step_ids <- c(
    "upload", "settings", "rename", "qc", "consensus",
    "blast", "taxonomy", "multilocus", "export"
  )
  lapply(workflow_chip_step_ids, function(step_id) {
    local({
      sid <- step_id
      observeEvent(input[[paste0("workflow_goto_", sid)]], {
        workflow_goto(sid, notify = TRUE)
      }, ignoreInit = TRUE)
    })
  })

  # Static Help button (not inside renderUI) — always works.
  observeEvent(input$workflow_open_help, {
    updateTabsetPanel(session, "pipeline_step", selected = "help")
    workflow_hide_loader()
  }, ignoreInit = TRUE)

  observe({
    help_on <- identical(as.character(input$pipeline_step)[1], "help")
    session$sendCustomMessage("setHelpActive", list(active = isTRUE(help_on)))
  })

  output$workflow_stage_heading <- renderUI({
    step <- if (!is.null(input$pipeline_step)) as.character(input$pipeline_step)[1] else "upload"
    pitax_workflow_stage_heading(step, rv$project_mode)
  })

  output$workflow_stage_actions <- renderUI({
    step <- if (!is.null(input$pipeline_step)) as.character(input$pipeline_step)[1] else "upload"
    simple <- identical(as.character(rv$project_mode)[1], "simple")
    to_consensus_label <- if (simple) "Continue to NCBI BLAST" else "Continue to Consensus"
    back_export_label <- if (simple) "Back to Trim & QC" else "Back to Consensus"

    actions <- switch(
      step,
      "upload" = tagList(
        div(class = "stage-topbar-spacer"),
        actionButton("to_settings", "Continue to Assay", icon = icon("arrow-right"), class = "btn-primary")
      ),
      "settings" = tagList(
        actionButton("back_upload", "Back", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("to_rename", "Continue to Assign", icon = icon("arrow-right"), class = "btn-primary")
      ),
      "rename" = tagList(
        actionButton("back_settings_from_rename", "Back to Assay", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("run_trimming", "Start trimming", icon = icon("play"), class = "btn-primary")
      ),
      "qc" = tagList(
        actionButton("back_rename_from_qc", "Back to Assign", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("to_consensus", to_consensus_label, icon = icon("arrow-right"), class = "btn-primary")
      ),
      "consensus" = tagList(
        actionButton("back_qc_from_consensus", "Back to QC", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("build_consensus", "Build / rebuild", icon = icon("cogs"), class = "btn-success"),
        actionButton("to_blast", "Continue to NCBI BLAST", icon = icon("arrow-right"), class = "btn-primary")
      ),
      "export" = tagList(
        actionButton("back_from_export", back_export_label, icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        span(class = "compact-hint", "Requires finished PROCESS: Trim & QC (Simple) or Consensus (Paired) with ready analysis sequences.")
      ),
      "blast" = tagList(
        actionButton("back_to_process", "Back to Process", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("to_taxonomy", "Continue to Taxonomic Summary", icon = icon("arrow-right"), class = "btn-primary"),
        actionButton("reset_pipeline", "Start new run")
      ),
      "taxonomy" = tagList(
        actionButton("back_blast", "Back to NCBI BLAST", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("run_taxonomy", "Analyze selected", icon = icon("play"), class = "btn-primary"),
        actionButton("run_taxonomy_all", "Analyze all retrieved", icon = icon("tasks"), class = "btn-success"),
        actionButton("to_multilocus", "Continue to Multi-locus", icon = icon("arrow-right")),
        actionButton("reset_pipeline_tax", "Start new run")
      ),
      "multilocus" = tagList(
        actionButton("back_taxonomy_from_multilocus", "Back to Taxonomic Summary", icon = icon("arrow-left")),
        div(class = "stage-topbar-spacer"),
        actionButton("build_multilocus_profile", "Build / rebuild profile", icon = icon("cogs"), class = "btn-success")
      ),
      "help" = tagList(
        div(class = "stage-topbar-spacer"),
        span(class = "compact-hint", "Documentation and scientific notes for the current PITAX workflow.")
      ),
      tagList(div(class = "stage-topbar-spacer"))
    )

    div(class = "stage-topbar workflow-stage-actions", actions)
  })
