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
                                   current_sample = NULL, header_config = pitax_header_default_config()) {
  header_config <- pitax_header_coerce_config(header_config)
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
      # Build/validate on the active Shiny session before browser navigation.
      # The actual download link is triggered only after a ZIP has been written.
      actionButton("export_prepare_download", "Download ZIP package", class = "btn-primary",
                   icon = icon("download")),
      tags$span(style = "display:none;",
                downloadButton("download_export_package", "Prepared ZIP", class = "btn-default"))
    ),
    div(
      class = "export-dialog",
      if (length(warnings)) div(class = "status-warning", lapply(warnings, tags$p)),
      div(class = "compact-hint", textOutput("export_download_status")),
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
            "One FASTA per read / record (full final name)" = "per_record",
            "Combined FASTA" = "combined",
            "One FASTA per isolate (grouped)" = "per_isolate",
            "One FASTA per locus (grouped)" = "per_locus"
          ),
          selected = "per_record"
        )
      ),
      div(class = "compact-hint", "One FASTA per read uses the full assigned Final_Name (for example 320_TEF_R.fasta). Grouped FASTA files use the isolate or locus as their filename, with full record names inside the FASTA."),
      div(
        class = "form-grid-2",
        selectInput(
          "export_sequence_header", "FASTA headers",
          choices = c("Short headers" = "short", "Headers with metadata" = "rich",
                      "Custom FASTA headers (template + per-record fields)" = "custom"),
          selected = "short"
        ),
        checkboxInput("export_include_plots", "Include selected plots as PNG", value = TRUE)
      ),
      conditionalPanel(
        condition = "input.export_sequence_header === 'custom'",
        div(class = "fasta-header-builder",
          div(class = "subsection-title", "Custom FASTA Header Builder"),
          p("Build reusable one-line FASTA headers using {Field} tokens. This changes only exported headers, not the original read/analysis identity or filename."),
          div(class = "form-grid-2",
            selectInput("export_header_saved_template", "Saved templates",
                        choices = names(header_config$templates), selected = header_config$active_template),
            textInput("export_header_template_name", "Template name", value = header_config$active_template)
          ),
          div(class = "button-row",
            actionButton("export_header_load_template", "Load selected template", icon = icon("folder-open")),
            actionButton("export_header_save_template", "Save / update template", icon = icon("save")),
            downloadButton("download_fasta_header_templates", "Export templates CSV"),
            fileInput("import_fasta_header_templates", "Import templates CSV", accept = ".csv")
          ),
          textAreaInput("export_header_template", "Header template",
                        value = header_config$draft_template, width = "100%", rows = 2),
          div(class = "compact-hint",
              "Built-in: {Final_Name}, {Isolate}, {Locus}, {Direction}, {Source_ID}, {Record_ID}. ",
              "Custom defaults: {Organism}, {Experiment}, {Seq_Date}. ",
              "Include {Final_Name} for unique identifiers. Date of sequencing is user-supplied, never inferred from export time."),
          textOutput("export_header_validation"),
          div(class = "compact-hint",
              "The table lists available metadata columns, not mandatory values. ",
              "Only tokens actually present in the active template are required. ",
              "Remove experiment=\"{Experiment}\" from the template to omit it entirely."),
          div(class = "subsection-title", "Per-record metadata | double-click a metadata cell to edit"),
          div(class = "button-row",
            textInput("export_header_new_field", "New metadata field (e.g. Fungicide)", value = ""),
            actionButton("export_header_add_field", "Add field", icon = icon("plus"))
          ),
          div(class = "form-grid-3",
            selectInput("export_header_bulk_field", "Field to set", choices = header_config$fields,
                        selected = if ("Experiment" %in% header_config$fields) "Experiment" else header_config$fields[1]),
            textInput("export_header_bulk_value", "Value for selected records", value = ""),
            selectInput("export_header_bulk_scope", "Apply to", 
                        choices = c("Selected table rows" = "selected", "All records" = "all"), selected = "selected")
          ),
          actionButton("export_header_apply_bulk", "Apply metadata value", icon = icon("check")),
          div(class = "compact-hint", "Select rows in the table for bulk edits; search/filter the table first. Values are stored with immutable record keys in the saved PITAX project."),
          div(class = "fasta-header-table-scroll", DTOutput("export_header_values_table")),
          div(class = "subsection-title", "Live header preview | review errors before download"),
          div(class = "fasta-header-table-scroll", DTOutput("export_header_preview_table"))
        )
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
