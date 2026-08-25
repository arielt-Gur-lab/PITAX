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
  step_ids <- vapply(catalog, function(s) s$id, character(1))
  cur_idx <- match(current, step_ids)
  current_category <- {
    hit <- Filter(function(s) identical(s$id, current), catalog)
    if (length(hit)) hit[[1]]$category else NA_character_
  }

  render_step_chip <- function(step) {
    is_current <- identical(step$id, current)
    is_unlocked <- isTRUE(free_nav) || step$id %in% unlocked || is_current
    step_idx <- match(step$id, step_ids)
    is_done <- (!is_current &&
                  (step$id %in% completed ||
                     (isTRUE(free_nav) && !is.na(cur_idx) && !is.na(step_idx) && step_idx < cur_idx)))
    state <- if (is_current) {
      "current"
    } else if (is_done) {
      "done"
    } else if (is_unlocked) {
      "available"
    } else {
      "locked"
    }
    tags$button(
      id = paste0("workflow_goto_", step$id),
      class = paste("workflow-chip", state),
      type = "button",
      disabled = if (!is_unlocked) TRUE else NULL,
      `aria-disabled` = if (!is_unlocked) "true" else "false",
      `aria-current` = if (is_current) "step" else NULL,
      `data-step` = step$id,
      span(class = "workflow-chip-num", as.character(step$number)),
      span(class = "workflow-chip-label", step$label)
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

  output_col <- div(
    class = paste("workflow-col", "workflow-col-output", if (export_current) "is-active" else NULL),
    div(class = "workflow-col-label", "OUTPUT"),
    div(
      class = "workflow-col-steps",
      tags$button(
        id = "workflow_goto_export",
        class = paste(
          "workflow-chip",
          if (export_current) "current" else if (export_ok) "available" else "locked"
        ),
        type = "button",
        disabled = if (!export_ok) TRUE else NULL,
        `aria-disabled` = if (!export_ok) "true" else "false",
        `data-step` = "export",
        icon("download"),
        span(class = "workflow-chip-label", "Export")
      )
    )
  )

  columns <- c(
    lapply(pitax_workflow_category_order, render_category_column),
    list(output_col)
  )

  # Help lives as a static actionButton beside this uiOutput (app_ui.R),
  # so stepper re-renders cannot break it.
  div(
    class = "workflow-nav",
    div(class = "workflow-nav-grid", columns)
  )
}
