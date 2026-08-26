# ============================================================
# UI helpers
# ============================================================

panel_box <- function(...) div(class = "panel-box", ...)
section_title <- function(x) div(class = "section-title", x)
info_tip <- function(text) tags$span(class = "info-tip", title = text, "(i)")

stage_heading <- function(icon_name, title, subtitle, badge = NULL) {
  div(
    class = "stage-heading",
    div(class = "stage-heading-icon", icon(icon_name)),
    div(
      class = "stage-heading-copy",
      h3(title),
      div(class = "stage-heading-subtitle", subtitle)
    ),
    if (!is.null(badge)) div(class = "stage-heading-badge", badge) else NULL
  )
}

pitax_workflow_stage_heading <- function(step, project_mode = "simple") {
  simple <- identical(as.character(project_mode)[1], "simple")
  meta <- switch(
    as.character(step)[1],
    "upload" = list("upload", "Upload chromatograms", "Keep the sequencer barcode and source filename unchanged. Biological identity is assigned later in Assign.", "SETUP | 1"),
    "settings" = list("sliders", "Assay setup", "Define one or more assay profiles and shared project trimming defaults. Trimming starts only after Assign.", "SETUP | 2"),
    "rename" = list("tags", "Assign read identity", "Assign isolate, assay and Forward/Reverse direction here, after Upload and Assay and before trimming. The upload barcode remains unchanged.", "SETUP | 3"),
    "qc" = list("bar-chart", "Trimming results, QC & curation", "Review the completed trim, inspect assigned chromatograms, and document manual sequence curation.", "PROCESS | 4"),
    "consensus" = list("random", "Forward/Reverse consensus", "Build the auditable analysis sequence used for BLAST according to the Paired project read model.", "PROCESS | 5"),
    "export" = list("download", "Export analysis sequences", "Export unlocks after PROCESS is complete: Trim & QC in Simple mode, or Consensus in Paired Forward/Reverse mode, once analysis sequences are ready.", "OUTPUT"),
    "blast" = list("search", "NCBI BLAST workspace", "Submit analysis sequences, retrieve accession-level hits, and keep each RID linked to the active sequence revision.", if (simple) "IDENTIFY | 5" else "IDENTIFY | 6"),
    "taxonomy" = list("sitemap", "Taxonomic interpretation", "Identify the best molecular match, inspect close alternatives and report the most conservative supported taxonomic level.", if (simple) "IDENTIFY | 6" else "IDENTIFY | 7"),
    "multilocus" = list("th", "Multi-locus isolate profile", "Integrate Isolate + Locus evidence from the current multi-locus project and/or imported projects, without flat voting.", if (simple) "INTEGRATE | 7" else "INTEGRATE | 8"),
    "help" = list("question-circle", "Help / About", "Documentation for the laboratory workflow, BLAST/taxonomy interpretation logic, and the scientific sources used to guide the application.", "DOCS"),
    list("flask", "PITAX", "Taxonomic identification workspace.", NULL)
  )
  stage_heading(meta[[1]], meta[[2]], meta[[3]], meta[[4]])
}

card_title <- function(title, tip = NULL, icon_name = NULL) {
  div(
    class = "card-title-row",
    if (!is.null(icon_name)) div(class = "card-title-icon", icon(icon_name)) else NULL,
    div(class = "card-title-text", title),
    if (!is.null(tip)) info_tip(tip) else NULL
  )
}

rbind_fill <- function(a, b) {
  if (is.null(a) || !is.data.frame(a) || !nrow(a)) return(b)
  if (is.null(b) || !is.data.frame(b) || !nrow(b)) return(a)
  cols <- union(names(a), names(b))
  for (nm in setdiff(cols, names(a))) a[[nm]] <- NA
  for (nm in setdiff(cols, names(b))) b[[nm]] <- NA
  a <- a[, cols, drop = FALSE]
  b <- b[, cols, drop = FALSE]
  rbind(a, b)
}

stage_topbar <- function(...) {
  div(class = "stage-topbar", ...)
}

# Kept as a no-op stub so older help/docs references do not crash if called.
pipeline_stage_footer <- function(current_step = NULL) {
  NULL
}

pitax_workflow_catalog <- function(project_mode = "simple") {
  paired <- identical(as.character(project_mode)[1], "paired_consensus")
  steps <- list(
    list(id = "upload", label = "Upload", number = 1L, category = "SETUP"),
    list(id = "settings", label = "Assay", number = 2L, category = "SETUP"),
    list(id = "rename", label = "Assign", number = 3L, category = "SETUP"),
    list(id = "qc", label = "Trim & QC", number = 4L, category = "PROCESS")
  )
  if (paired) {
    steps[[length(steps) + 1L]] <- list(id = "consensus", label = "Consensus", number = 5L, category = "PROCESS")
    blast_n <- 6L
    tax_n <- 7L
    multi_n <- 8L
  } else {
    blast_n <- 5L
    tax_n <- 6L
    multi_n <- 7L
  }
  steps[[length(steps) + 1L]] <- list(id = "blast", label = "BLAST", number = blast_n, category = "IDENTIFY")
  steps[[length(steps) + 1L]] <- list(id = "taxonomy", label = "Taxonomy", number = tax_n, category = "IDENTIFY")
  steps[[length(steps) + 1L]] <- list(id = "multilocus", label = "Multi-locus", number = multi_n, category = "INTEGRATE")
  steps
}

pitax_workflow_category_order <- c("SETUP", "PROCESS", "IDENTIFY", "INTEGRATE")

pitax_workflow_stepper_ui <- function(current, unlocked, project_mode = "simple", completed = character(), free_nav = TRUE) {
  current <- as.character(current)[1]
  unlocked <- unique(as.character(unlocked))
  completed <- unique(as.character(completed))
  catalog <- pitax_workflow_catalog(project_mode)
  current_category <- {
    hit <- Filter(function(s) identical(s$id, current), catalog)
    if (length(hit)) hit[[1]]$category else NA_character_
  }

  render_step_chip <- function(step) {
    is_current <- identical(step$id, current)
    is_unlocked <- isTRUE(free_nav) || step$id %in% unlocked || is_current
    # Green = real completion only (caller supplies completed); never "visited earlier".
    is_done <- !is_current && step$id %in% completed
    state <- if (is_current) {
      "current"
    } else if (is_done) {
      "done"
    } else if (is_unlocked) {
      "available"
    } else {
      "locked"
    }
    # Shiny actionButton - no custom JS / overlay path for navigation.
    actionButton(
      inputId = paste0("workflow_goto_", step$id),
      label = tagList(
        span(class = "workflow-chip-num", as.character(step$number)),
        span(class = "workflow-chip-label", step$label)
      ),
      class = paste("workflow-chip", state),
      `aria-current` = if (is_current) "step" else NULL
    )
  }

  render_category_column <- function(cat) {
    cat_steps <- Filter(function(s) identical(s$category, cat), catalog)
    if (!length(cat_steps)) return(NULL)
    chips <- lapply(cat_steps, render_step_chip)
    div(
      class = paste("workflow-col", if (identical(cat, current_category)) "is-active" else NULL),
      div(class = "workflow-col-label", cat),
      div(class = "workflow-col-steps", chips)
    )
  }

  export_current <- identical(current, "export")
  export_ok <- isTRUE(free_nav) || "export" %in% unlocked || export_current
  export_done <- !export_current && "export" %in% completed
  export_state <- if (export_current) {
    "current"
  } else if (export_done) {
    "done"
  } else if (export_ok) {
    "available"
  } else {
    "locked"
  }

  output_col <- div(
    class = paste("workflow-col", "workflow-col-output", if (export_current) "is-active" else NULL),
    div(class = "workflow-col-label", "OUTPUT"),
    div(
      class = "workflow-col-steps",
      actionButton(
        inputId = "workflow_goto_export",
        label = tagList(
          icon("download"),
          span(class = "workflow-chip-label", "Export")
        ),
        class = paste("workflow-chip", export_state)
      )
    )
  )

  columns <- c(
    lapply(pitax_workflow_category_order, render_category_column),
    list(output_col)
  )

  # Help lives as a static actionButton in the app header (under the version badge),
  # so stepper re-renders cannot break it.
  div(
    class = "workflow-nav",
    div(class = "workflow-nav-grid", columns)
  )
}
