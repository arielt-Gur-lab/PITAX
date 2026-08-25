# ============================================================
# Free workflow navigation (stepper is orientation, not a gate)
# ============================================================

  # rv$workflow_unlocked is initialized in 00_state.R; kept for progress hints /
  # project load compatibility. Navigation itself is not gated on it.

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

  workflow_can_visit <- function(step) {
    step <- as.character(step)[1]
    if (!nzchar(step)) return(FALSE)
    if (identical(step, "help")) return(TRUE)
    if (identical(step, "export")) return(TRUE)
    if (identical(step, "consensus") && identical(isolate(rv$project_mode), "simple")) return(FALSE)
    step %in% workflow_all_step_ids()
  }

  workflow_goto <- function(step, notify = TRUE) {
    step <- as.character(step)[1]
    if (identical(step, "help")) {
      updateTabsetPanel(session, "pipeline_step", selected = "help")
      session$sendCustomMessage("hideLoader", list())
      return(TRUE)
    }
    if (identical(step, "consensus") && identical(isolate(rv$project_mode), "simple")) {
      if (isTRUE(notify)) {
        showNotification("Consensus is only used in Paired Forward/Reverse mode.", type = "message", duration = 5)
      }
      return(FALSE)
    }
    if (!workflow_can_visit(step)) {
      if (isTRUE(notify)) {
        showNotification("Unknown workflow step.", type = "warning", duration = 5)
      }
      return(FALSE)
    }
    updateTabsetPanel(session, "pipeline_step", selected = step)
    TRUE
  }

  output$workflow_stepper <- renderUI({
    mode <- rv$project_mode
    catalog_ids <- vapply(pitax_workflow_catalog(mode), function(s) s$id, character(1))
    current <- if (!is.null(input$pipeline_step)) input$pipeline_step else "upload"
    cur_idx <- match(current, catalog_ids)
    completed <- if (!is.na(cur_idx) && cur_idx > 1) catalog_ids[seq_len(cur_idx - 1L)] else character()
    pitax_workflow_stepper_ui(
      current = current,
      unlocked = c(catalog_ids, "export", "help"),
      project_mode = mode,
      completed = completed,
      free_nav = TRUE
    )
  })

  # Static Help button (not inside renderUI) — always works.
  observeEvent(input$workflow_open_help, {
    updateTabsetPanel(session, "pipeline_step", selected = "help")
    session$sendCustomMessage("hideLoader", list())
  }, ignoreInit = TRUE)

  observe({
    help_on <- identical(as.character(input$pipeline_step)[1], "help")
    session$sendCustomMessage("setHelpActive", list(active = isTRUE(help_on)))
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
        span(class = "compact-hint", "Export is available from OUTPUT whenever analysis sequences are ready.")
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

  observeEvent(input$workflow_nav_step, {
    step <- as.character(input$workflow_nav_step)[1]
    if (!nzchar(step)) return()
    if (identical(step, "help")) {
      # Prefer the static Help button path; still honour chip clicks if present.
      updateTabsetPanel(session, "pipeline_step", selected = "help")
      session$sendCustomMessage("hideLoader", list())
      return()
    }
    workflow_goto(step, notify = TRUE)
  }, ignoreInit = TRUE)
