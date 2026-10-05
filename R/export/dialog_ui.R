# ============================================================
# Shared export dialog UI helpers
# ============================================================

pitax_export_dialog_ui <- function(stage_id = "export", availability = list(),
                                   default_ids = character(), preview = NULL,
                                   file_plan = list(), warnings = character(),
                                   sample_choices = character(),
                                   isolate_choices = character(),
                                   locus_choices = character(),
                                   assay_choices = character(),
                                   current_sample = NULL) {
  status_badge <- function(status) {
    cls <- switch(
      as.character(status)[1],
      "available" = "export-status-available",
      "partial" = "export-status-partial",
      "not_created" = "export-status-missing",
      "stale" = "export-status-stale",
      "blocked" = "export-status-blocked",
      "export-status-missing"
    )
    tags$span(class = paste("export-status-badge", cls), as.character(status))
  }

  rows <- list()
  for (id in pitax_export_component_ids()) {
    comp <- pitax_export_get_component(id)
    av <- availability[[id]]
    if (is.null(av)) next
    # Stage-scoped view still lists only stage components; project scope lists all.
    selectable <- av$status %in% pitax_export_selectable_statuses()
    rows[[length(rows) + 1L]] <- tags$div(
      class = "export-component-row",
      tags$label(
        class = "export-component-check",
        checkboxInput(
          inputId = paste0("export_comp_", id),
          label = NULL,
          value = isTRUE(id %in% default_ids) && selectable,
          width = "auto"
        )
      ),
      tags$div(
        class = "export-component-meta",
        tags$div(class = "export-component-title", comp$label, status_badge(av$status)),
        tags$div(class = "export-component-reason", av$reason),
        tags$div(class = "export-component-count", paste0("Records: ", av$record_count))
      )
    )
    # Disable non-selectable checkboxes via JS-friendly class
    if (!selectable) {
      rows[[length(rows)]] <- tagAppendAttributes(rows[[length(rows)]], class = "is-disabled")
    }
  }

  modalDialog(
    title = paste0("Export | source stage: ", stage_id),
    size = "l",
    easyClose = TRUE,
    footer = tagList(
      modalButton("Cancel"),
      downloadButton("download_export_package", "Download ZIP package", class = "btn-primary")
    ),
    div(
      class = "export-dialog",
      if (length(warnings)) div(class = "status-warning", lapply(warnings, tags$p)),
      radioButtons(
        "export_scope_mode", "Data scope",
        choices = c(
          "Current stage data" = "stage",
          "All available project data" = "project"
        ),
        selected = if (identical(stage_id, "export")) "project" else "stage",
        inline = TRUE
      ),
      div(
        class = "form-grid-2",
        selectInput(
          "export_sample_mode", "Samples",
          choices = c("All samples" = "all", "Current sample" = "current", "Selected samples" = "selected"),
          selected = "all"
        ),
        selectizeInput(
          "export_selected_samples", "Sample selection",
          choices = sample_choices,
          selected = if (!is.null(current_sample)) current_sample else NULL,
          multiple = TRUE
        )
      ),
      div(
        class = "form-grid-3",
        selectizeInput("export_filter_isolates", "Isolates", choices = isolate_choices, multiple = TRUE),
        selectizeInput("export_filter_loci", "Loci", choices = locus_choices, multiple = TRUE),
        selectizeInput("export_filter_assays", "Assays", choices = assay_choices, multiple = TRUE)
      ),
      div(
        class = "form-grid-2",
        selectInput(
          "export_excel_layout", "Excel layout",
          choices = c(
            "One sheet per data type" = "sheet_per_component",
            "Single sheet for one table" = "single_sheet",
            "One sheet per sample" = "sheet_per_sample",
            "One workbook per isolate" = "workbook_per_isolate"
          ),
          selected = "sheet_per_component"
        ),
        selectInput(
          "export_sequence_layout", "Sequence files",
          choices = c(
            "Combined FASTA" = "combined",
            "One file per record" = "per_record",
            "One file per isolate" = "per_isolate",
            "One file per locus" = "per_locus"
          ),
          selected = "combined"
        )
      ),
      div(
        class = "form-grid-2",
        selectInput(
          "export_sequence_header", "FASTA headers",
          choices = c("Short headers" = "short", "Headers with metadata" = "rich"),
          selected = "short"
        ),
        checkboxInput("export_include_plots", "Include selected plots as PNG", value = TRUE)
      ),
      div(class = "subsection-title", "Components"),
      div(class = "export-component-list", rows),
      div(class = "subsection-title", "Preview"),
      if (is.data.frame(preview) && nrow(preview)) {
        tags$table(
          class = "table table-condensed export-preview-table",
          tags$thead(tags$tr(lapply(names(preview), tags$th))),
          tags$tbody(lapply(seq_len(nrow(preview)), function(i) {
            tags$tr(lapply(as.character(unlist(preview[i, , drop = TRUE])), tags$td))
          }))
        )
      } else {
        div(class = "compact-hint", "No exportable components are selected.")
      },
      div(class = "subsection-title", "Planned files"),
      tags$ul(lapply(file_plan, function(f) tags$li(f$path)))
    )
  )
}
