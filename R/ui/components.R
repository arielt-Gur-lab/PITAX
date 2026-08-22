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

pipeline_stage_footer <- function(current_step) {
  labels <- c(
    "Upload",
    "Assay settings",
    "Rename & assign",
    "Trim & QC",
    "Analysis sequence",
    "Export",
    "NCBI BLAST",
    "Taxonomic summary",
    "Multi-locus profile"
  )

  pieces <- list()
  for (i in seq_along(labels)) {
    state <- if (i < current_step) "done" else if (i == current_step) "current" else "future"
    pieces[[length(pieces) + 1]] <- div(
      class = paste("schema-stage", state),
      div(class = "schema-dot"),
      div(class = "schema-label", paste0(i, ". ", labels[[i]]))
    )
    if (i < length(labels)) {
      pieces[[length(pieces) + 1]] <- div(
        class = paste("schema-line", if (i < current_step) "done" else "")
      )
    }
  }

  div(
    class = "pipeline-schema-wrap",
    div(class = "pipeline-schema-title", "Workflow position"),
    div(class = "pipeline-schema", pieces)
  )
}

