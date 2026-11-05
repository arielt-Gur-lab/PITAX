# Custom FASTA header builder: only active in the shared export dialog.
# Values are linked to immutable read/analysis IDs rather than parsing filenames.
  export_header_state <- reactiveValues(config = pitax_header_default_config())

  header_builder_config <- function() {
    config <- pitax_header_coerce_config(export_header_state$config)
    if (!is.null(input$export_header_template)) {
      config$draft_template <- as.character(input$export_header_template)[1]
    }
    config
  }
  header_builder_store <- function(config) {
    config <- pitax_header_coerce_config(config)
    export_header_state$config <- config
    rv$fasta_header_config <- config
    invisible(config)
  }

  output$export_header_values_table <- DT::renderDT({
    req(export_dialog_state$snap)
    config <- pitax_header_coerce_config(export_header_state$config)
    df <- pitax_header_edit_table(export_dialog_state$snap, config)
    df <- df[, c("Record_Key", "Final_Name", "Isolate", "Locus", "Direction", config$fields), drop = FALSE]
    DT::datatable(
      df, rownames = FALSE,
      selection = list(mode = "multiple", selected = NULL, target = "row"),
      editable = list(target = "cell", disable = list(columns = 0:4)),
      options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100),
                     scrollX = FALSE, autoWidth = TRUE, dom = "ftip")
    )
  }, server = FALSE)

  observeEvent(input$export_header_values_table_cell_edit, {
    edit <- input$export_header_values_table_cell_edit
    config <- pitax_header_coerce_config(export_header_state$config)
    df <- pitax_header_edit_table(export_dialog_state$snap, config)
    columns <- c("Record_Key", "Final_Name", "Isolate", "Locus", "Direction", config$fields)
    i <- suppressWarnings(as.integer(edit$row)); j <- suppressWarnings(as.integer(edit$col)) + 1L
    if (length(i) != 1L || length(j) != 1L || is.na(i) || is.na(j) ||
        i < 1L || i > nrow(df) || j < 6L || j > length(columns)) return()
    ans <- tryCatch(pitax_header_set_value(config, df$Record_Key[i], columns[j], edit$value),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_add_field, {
    config <- pitax_header_coerce_config(export_header_state$config)
    field <- tryCatch(pitax_header_validate_field(input$export_header_new_field, config$fields),
                      error = function(e) e)
    if (inherits(field, "error")) {
      showNotification(conditionMessage(field), type = "error", duration = 8)
      return()
    }
    if (length(config$fields) >= 40L) {
      showNotification("Maximum 40 custom fields per project.", type = "error")
      return()
    }
    config$fields <- c(config$fields, field)
    config$values[[field]] <- rep("", nrow(config$values))
    header_builder_store(config)
    updateSelectInput(session, "export_header_bulk_field", choices = config$fields, selected = field)
    updateTextInput(session, "export_header_new_field", value = "")
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_apply_bulk, {
    config <- pitax_header_coerce_config(export_header_state$config)
    df <- pitax_header_edit_table(export_dialog_state$snap, config)
    field <- as.character(input$export_header_bulk_field %||% "")[1]
    scope <- as.character(input$export_header_bulk_scope %||% "selected")[1]
    selected <- as.integer(input$export_header_values_table_rows_selected %||% integer())
    if (identical(scope, "selected")) {
      selected <- selected[!is.na(selected) & selected >= 1L & selected <= nrow(df)]
      if (!length(selected)) {
        showNotification("Select rows in the metadata table, or choose All records.", type = "warning", duration = 8)
        return()
      }
      keys <- df$Record_Key[selected]
    } else keys <- df$Record_Key
    ans <- tryCatch(pitax_header_set_value(config, keys, field, input$export_header_bulk_value),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans)
    showNotification(paste("Set", field, "for", length(keys), "record(s)."), type = "message", duration = 4)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_template, {
    config <- pitax_header_coerce_config(export_header_state$config)
    template <- as.character(input$export_header_template)[1]
    if (is.na(template)) return()
    config$draft_template <- template
    header_builder_store(config)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_save_template, {
    config <- header_builder_config()
    name <- trimws(as.character(input$export_header_template_name %||% "")[1])
    if (is.na(name) || !nzchar(name) || nchar(name) > 80L) {
      showNotification("Enter a template name (1-80 characters).", type = "error")
      return()
    }
    error <- tryCatch({pitax_header_validate_template(config$draft_template, config$fields); NULL},
                      error = function(e) conditionMessage(e))
    if (!is.null(error)) {
      showNotification(error, type = "error", duration = 8)
      return()
    }
    config$templates[[name]] <- config$draft_template
    config$active_template <- name
    header_builder_store(config)
    updateSelectInput(session, "export_header_saved_template", choices = names(config$templates), selected = name)
    showNotification(paste("Saved template:", name), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_load_template, {
    config <- pitax_header_coerce_config(export_header_state$config)
    name <- as.character(input$export_header_saved_template %||% "")[1]
    if (!name %in% names(config$templates)) return()
    config$active_template <- name
    config$draft_template <- config$templates[[name]]
    header_builder_store(config)
    updateTextInput(session, "export_header_template_name", value = name)
    updateTextAreaInput(session, "export_header_template", value = config$draft_template)
  }, ignoreInit = TRUE)

  output$export_header_preview_table <- DT::renderDT({
    req(export_dialog_state$snap)
    config <- pitax_header_coerce_config(export_header_state$config)
    if (!is.null(input$export_header_template)) config$draft_template <- as.character(input$export_header_template)[1]
    df <- pitax_header_catalog(export_dialog_state$snap)
    preview <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, , drop = FALSE]
      key <- as.character(row$Record_Key)
      kind <- sub("::.*$", "", key)
      record_id <- sub("^[^:]+::", "", key)
      comp <- switch(kind, read = "qc_trimmed_sequences",
                     analysis = "consensus_analysis_sequences", multilocus = "multilocus_sequences")
      rec <- list(final_name = row$Final_Name, isolate = row$Isolate,
                  locus = row$Locus, direction = row$Direction)
      result <- tryCatch(pitax_header_render(rec, comp, record_id, config), error = function(e) e)
      data.frame(Record = as.character(row$Final_Name),
                 Header = if (inherits(result, "error")) "" else paste0(">", result),
                 Status = if (inherits(result, "error")) conditionMessage(result) else "Ready",
                 stringsAsFactors = FALSE)
    })
    out <- if (length(preview)) do.call(rbind, preview) else
      data.frame(Record = "", Header = "", Status = "No assigned reads", stringsAsFactors = FALSE)
    # Within each sequence type, flag identifiers that would collide in FASTA.
    if (nrow(df) && nrow(out) == nrow(df)) {
      categories <- sub("::.*$", "", df$Record_Key)
      identifiers <- ifelse(nzchar(out$Header),
                           sub("[[:space:]].*$", "", sub("^>", "", out$Header)), "")
      for (category in unique(categories)) {
        subset <- which(categories == category & nzchar(identifiers) & out$Status == "Ready")
        if (length(subset) > 1L) {
          duplicates <- duplicated(tolower(identifiers[subset])) |
                        duplicated(tolower(identifiers[subset]), fromLast = TRUE)
          out$Status[subset[duplicates]] <- "Duplicate primary FASTA ID; include {Final_Name}"
        }
      }
    }
    DT::datatable(out, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = FALSE, autoWidth = TRUE,
                                 dom = "ftip", columnDefs = list(list(targets = 1, className = "dt-left"))))
  }, server = FALSE)

  output$export_header_validation <- renderText({
    req(export_dialog_state$snap)
    config <- header_builder_config()
    error <- tryCatch({pitax_header_validate_template(config$draft_template, config$fields); NULL},
                      error = function(e) conditionMessage(e))
    if (!is.null(error)) paste("Template error:", error) else
      "Preview updates after editing. Blank required variables are errors: fill metadata before downloading custom FASTA."
  })

  output$download_fasta_header_templates <- downloadHandler(
    filename = function() "PITAX_FASTA_header_templates.csv",
    content = function(file) {
      config <- pitax_header_coerce_config(export_header_state$config)
      utils::write.csv(data.frame(Name = names(config$templates),
                                  Template = unlist(config$templates, use.names = FALSE),
                                  Fields = paste(config$fields, collapse = ";"),
                                  stringsAsFactors = FALSE),
                       file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  observeEvent(input$import_fasta_header_templates, {
    config <- pitax_header_coerce_config(export_header_state$config)
    uploaded <- input$import_fasta_header_templates
    req(uploaded$datapath)
    table <- tryCatch(utils::read.csv(uploaded$datapath, stringsAsFactors = FALSE,
                                      check.names = FALSE, fileEncoding = "UTF-8"),
                      error = function(e) e)
    if (inherits(table, "error") || !all(c("Name", "Template") %in% names(table))) {
      showNotification("Invalid template CSV; expected Name and Template columns.", type = "error")
      return()
    }
    # Portable template CSV carries field definitions; values remain project-specific.
    if ("Fields" %in% names(table)) {
      imported_fields <- unique(unlist(strsplit(paste(table$Fields, collapse = ";"), ";", fixed = TRUE)))
      for (field in imported_fields) {
        field <- trimws(as.character(field))
        if (!nzchar(field) || field %in% config$fields || length(config$fields) >= 40L) next
        valid <- tryCatch(pitax_header_validate_field(field, config$fields), error = function(e) NULL)
        if (is.null(valid)) next
        config$fields <- c(config$fields, valid)
        config$values[[valid]] <- rep("", nrow(config$values))
      }
    }
    for (i in seq_len(nrow(table))) {
      nm <- trimws(as.character(table$Name[i]))
      tmpl <- as.character(table$Template[i])
      if (is.na(nm) || !nzchar(nm) || nchar(nm) > 80L) next
      if (is.na(tmpl) || !nzchar(tmpl) || grepl("[\r\n]", tmpl)) next
      config$templates[[nm]] <- tmpl
    }
    header_builder_store(config)
    updateSelectInput(session, "export_header_bulk_field", choices = config$fields,
                      selected = config$fields[1])
    updateSelectInput(session, "export_header_saved_template", choices = names(config$templates),
                      selected = config$active_template)
    showNotification("Imported FASTA header templates and field names. Populate metadata values before export.", type = "message")
  }, ignoreInit = TRUE)
