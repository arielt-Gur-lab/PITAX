# Structured Custom FASTA Header Builder (profile-based).
# Single authoritative state: export_header_state$config / rv$fasta_header_config.

  export_header_state <- reactiveValues(config = pitax_header_default_config())
  header_browser_restore_done <- reactiveVal(FALSE)

  header_builder_config <- function() {
    pitax_header_coerce_config(export_header_state$config)
  }

  header_builder_store <- function(config, mark_dirty = NULL) {
    config <- pitax_header_coerce_config(config)
    if (!is.null(mark_dirty)) config$dirty <- isTRUE(mark_dirty)
    export_header_state$config <- config
    rv$fasta_header_config <- config
    header_builder_push_ui(config)
    invisible(config)
  }

  header_builder_push_ui <- function(config = NULL) {
    config <- pitax_header_coerce_config(config %||% header_builder_config())
    # Keep hidden textarea / selects synchronized for diagnostics and legacy hooks.
    updateTextAreaInput(session, "export_header_template", value = config$draft_template)
    updateTextInput(session, "export_header_template_name", value = config$draft$name)
    updateSelectInput(session, "export_header_tag_style", selected = config$draft$tag_style)
    profile_choices <- setNames(
      vapply(config$profiles, function(p) p$id, character(1)),
      vapply(config$profiles, function(p) p$name, character(1))
    )
    updateSelectInput(session, "export_header_saved_template",
                      choices = profile_choices, selected = config$active_profile_id)
    active <- pitax_header_active_metadata_fields(config)
    labels <- setNames(active, vapply(active, function(id) {
      pitax_header_field_def(id, config)$display
    }, character(1)))
    selected_bulk <- isolate(input$export_header_bulk_field %||% "")
    if (!length(selected_bulk) || !selected_bulk %in% active) {
      selected_bulk <- if (length(active)) active[1] else character()
    }
    updateSelectInput(session, "export_header_bulk_field",
                      choices = labels, selected = selected_bulk)
    custom_ids <- names(config$custom_fields)
    updateSelectInput(session, "export_header_existing_field",
                      choices = custom_ids,
                      selected = if (length(custom_ids)) custom_ids[1] else character())
    session$sendCustomMessage("pitaxHeaderBuilderState", pitax_header_builder_state_json(config))
    invisible(config)
  }

  # Push builder state whenever the export dialog snapshot is available / config changes.
  observe({
    req(export_dialog_state$snap)
    header_builder_push_ui(header_builder_config())
    cat <- pitax_header_catalog(export_dialog_state$snap)
    isolates <- unique(cat$Isolate[nzchar(cat$Isolate)])
    loci <- unique(cat$Locus[nzchar(cat$Locus)])
    updateSelectizeInput(session, "export_header_bulk_isolate",
                         choices = isolates, server = FALSE)
    updateSelectizeInput(session, "export_header_bulk_locus",
                         choices = loci, server = FALSE)
  })

  output$export_header_values_table <- DT::renderDT({
    req(export_dialog_state$snap)
    config <- header_builder_config()
    df <- pitax_header_active_edit_table(export_dialog_state$snap, config)
    export_header_state$edit_columns <- names(df)
    # Display labels for active field columns (identity columns stay as-is).
    for (fid in pitax_header_active_metadata_fields(config)) {
      if (fid %in% names(df)) {
        names(df)[names(df) == fid] <- pitax_header_field_def(fid, config)$display
      }
    }
    freeze <- which(export_header_state$edit_columns %in%
                      c("Record_Key", "Final_Name", "Isolate", "Locus", "Direction")) - 1L
    DT::datatable(
      df, rownames = FALSE,
      selection = list(mode = "multiple", selected = NULL, target = "row"),
      editable = list(target = "cell", disable = list(columns = freeze)),
      options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100),
                     scrollX = FALSE, autoWidth = TRUE, dom = "ftip")
    )
  }, server = FALSE)

  observeEvent(input$export_header_values_table_cell_edit, {
    edit <- input$export_header_values_table_cell_edit
    config <- header_builder_config()
    df <- pitax_header_active_edit_table(export_dialog_state$snap, config)
    cols <- export_header_state$edit_columns %||% names(df)
    i <- suppressWarnings(as.integer(edit$row))
    j <- suppressWarnings(as.integer(edit$col)) + 1L
    if (length(i) != 1L || length(j) != 1L || is.na(i) || is.na(j) ||
        i < 1L || i > nrow(df) || j < 6L || j > length(cols)) return()
    field_id <- cols[j]
    if (!field_id %in% pitax_header_active_metadata_fields(config)) return()
    ans <- tryCatch(pitax_header_set_value(config, df$Record_Key[i], field_id, edit$value),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
  }, ignoreInit = TRUE)

  output$export_header_registered_fields <- renderText({
    config <- header_builder_config()
    registry <- pitax_header_registered_field_table(config)
    if (!nrow(registry)) return("No registered metadata fields.")
    paste0("Registered metadata: ",
           paste(paste0(registry$Field,
                        ifelse(registry$Status == "In active header",
                               " [active]", " [available]")), collapse = ", "))
  })

  output$export_header_active_template <- renderText({
    config <- header_builder_config()
    paste0("Active export profile: ", config$draft$name, " | >", config$draft_template)
  })

  output$export_header_validation <- renderText({
    config <- header_builder_config()
    issues <- pitax_header_validate_profile_structure(config)
    active <- pitax_header_active_metadata_fields(config)
    required <- intersect(active, config$draft$required)
    parts <- character()
    if (length(issues)) parts <- c(parts, paste(issues, collapse = " "))
    if (length(required)) {
      parts <- c(parts, paste0("Required fields: ", paste(required, collapse = ", "), "."))
    } else {
      parts <- c(parts, "No required custom fields. Empty optional annotations are omitted.")
    }
    paste(parts, collapse = " ")
  })

  # ---- Drag/drop builder events from JS ----
  observeEvent(input$pitax_header_builder_add, {
    payload <- input$pitax_header_builder_add
    field_id <- payload$field_id %||% payload[["field_id"]]
    ans <- tryCatch(pitax_header_profile_add_field(header_builder_config(), field_id),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    # If assay triggered auto-add of locus, tell the user.
    if (identical(as.character(field_id)[1], "assay") &&
        "locus" %in% ans$draft$ordered_fields) {
      showNotification("Assay requires locus; locus was added to the profile.",
                       type = "message", duration = 7)
    }
    header_builder_store(ans, mark_dirty = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_header_builder_remove, {
    payload <- input$pitax_header_builder_remove
    field_id <- payload$field_id %||% payload[["field_id"]]
    ans <- tryCatch(pitax_header_profile_remove_field(header_builder_config(), field_id),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_header_builder_order, {
    payload <- input$pitax_header_builder_order
    fields <- unlist(payload$fields %||% list(), use.names = FALSE)
    ans <- tryCatch(pitax_header_profile_reorder_fields(header_builder_config(), fields),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      header_builder_push_ui()
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_tag_style, {
    config <- header_builder_config()
    style <- as.character(input$export_header_tag_style)[1]
    if (!style %in% c("bare", "bracket")) return()
    if (identical(config$draft$tag_style, style)) return()
    config$draft$tag_style <- style
    config$dirty <- TRUE
    header_builder_store(config, mark_dirty = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_add_field, {
    result <- tryCatch(
      pitax_header_create_custom_field(header_builder_config(),
                                       input$export_header_new_field,
                                       activate = TRUE),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      showNotification(conditionMessage(result), type = "error", duration = 8)
      return()
    }
    header_builder_store(result$config, mark_dirty = TRUE)
    updateTextInput(session, "export_header_new_field", value = "")
    showNotification(
      paste0(if (result$created) "Created" else "Reused",
             " custom field '", result$field, "' and added it to the active profile."),
      type = "message", duration = 7
    )
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_delete_field, {
    field <- input$export_header_existing_field %||% ""
    ans <- tryCatch(pitax_header_delete_custom_field(header_builder_config(), field),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 12)
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
    showNotification(paste0("Deleted custom field '", field, "'."), type = "warning", duration = 8)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_apply_bulk, {
    config <- header_builder_config()
    df <- pitax_header_active_edit_table(export_dialog_state$snap, config)
    field <- as.character(input$export_header_bulk_field %||% "")[1]
    scope <- as.character(input$export_header_bulk_scope %||% "selected")[1]
    value <- input$export_header_bulk_value
    keys <- character()
    ans <- tryCatch({
      if (identical(scope, "selected")) {
        selected <- as.integer(input$export_header_values_table_rows_selected %||% integer())
        selected <- selected[!is.na(selected) & selected >= 1L & selected <= nrow(df)]
        if (!length(selected)) stop("Select rows in the metadata table, or choose another scope.", call. = FALSE)
        keys <- df$Record_Key[selected]
        pitax_header_set_value(config, keys, field, value)
      } else if (identical(scope, "all")) {
        keys <- df$Record_Key
        pitax_header_set_value(config, keys, field, value)
      } else if (identical(scope, "isolate")) {
        iso <- as.character(input$export_header_bulk_isolate %||% "")[1]
        if (!nzchar(iso)) stop("Choose an isolate filter.", call. = FALSE)
        pitax_header_apply_isolate_value(config, export_dialog_state$snap, iso, field, value)
      } else if (identical(scope, "isolate_locus")) {
        iso <- as.character(input$export_header_bulk_isolate %||% "")[1]
        loc <- as.character(input$export_header_bulk_locus %||% "")[1]
        if (!nzchar(iso) || !nzchar(loc)) stop("Choose both isolate and locus filters.", call. = FALSE)
        pitax_header_apply_isolate_value(config, export_dialog_state$snap, iso, field, value, locus = loc)
      } else stop("Unknown bulk scope.", call. = FALSE)
    }, error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
    showNotification("Metadata value applied.", type = "message", duration = 4)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_save_template, {
    config <- header_builder_config()
    name <- trimws(as.character(input$export_header_template_name %||% "")[1])
    ans <- tryCatch(pitax_header_save_profile(config, name), error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    ans$persistence$project_included <- TRUE
    ans$persistence$last_status <- "Saved in current project state (use Save project for .sangerproject)"
    header_builder_store(ans, mark_dirty = FALSE)
    showNotification(
      paste0("Profile '", ans$draft$name, "' saved in the current project state. ",
             "Use 'Save to this browser' for browser persistence, or Save project for .sangerproject."),
      type = "message", duration = 10
    )
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_load_template, {
    config <- header_builder_config()
    id <- as.character(input$export_header_saved_template %||% "")[1]
    ans <- tryCatch(pitax_header_load_profile(config, id), error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = FALSE)
    showNotification(paste0("Loaded profile '", ans$draft$name, "'."), type = "message", duration = 5)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_duplicate_profile, {
    ans <- tryCatch(pitax_header_duplicate_profile(header_builder_config()), error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = FALSE)
    showNotification(paste0("Duplicated as '", ans$draft$name, "'."), type = "message", duration = 5)
  }, ignoreInit = TRUE)

  observeEvent(input$export_header_delete_profile, {
    id <- as.character(input$export_header_saved_template %||% "")[1]
    ans <- tryCatch(pitax_header_delete_profile(header_builder_config(), id), error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 8)
      return()
    }
    header_builder_store(ans, mark_dirty = FALSE)
    showNotification("Profile deleted.", type = "warning", duration = 5)
  }, ignoreInit = TRUE)

  # ---- Browser persistence (profiles only) ----
  observeEvent(input$export_header_save_browser, {
    config <- header_builder_config()
    # Ensure named draft is stored before browser save.
    config <- tryCatch(pitax_header_save_profile(config, input$export_header_template_name),
                       error = function(e) e)
    if (inherits(config, "error")) {
      showNotification(conditionMessage(config), type = "error", duration = 8)
      return()
    }
    header_builder_store(config, mark_dirty = FALSE)
    payload <- pitax_header_profiles_portable(config)
    payload$active_name <- config$draft$name
    session$sendCustomMessage("pitaxStoreHeaderProfiles", payload)
    showNotification("Writing profiles to this browser… waiting for verification.",
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_header_browser_save_result, {
    res <- input$pitax_header_browser_save_result
    config <- header_builder_config()
    if (isTRUE(res$ok)) {
      config$persistence$browser_saved <- TRUE
      config$persistence$last_browser_save_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      config$persistence$last_status <- paste0("Verified save in this browser at ",
                                               config$persistence$last_browser_save_at)
      header_builder_store(config, mark_dirty = FALSE)
      showNotification(
        paste0("Verified: profile layout '", res$name %||% config$draft$name,
               "' is saved in this browser. Record metadata remains in the project only."),
        type = "message", duration = 10
      )
    } else {
      config$persistence$browser_saved <- FALSE
      config$persistence$last_status <- as.character(res$error %||% "Browser save failed")[1]
      header_builder_store(config)
      showNotification(
        paste0(config$persistence$last_status,
               " Download the portable JSON profile as a backup."),
        type = "error", duration = 12
      )
    }
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_header_browser_profiles_restore, {
    payload <- input$pitax_header_browser_profiles_restore
    data <- payload$data
    if (is.null(data) || (!is.list(data))) {
      header_browser_restore_done(TRUE)
      return()
    }
    # Shiny may deliver JS arrays as lists; normalize profiles to list-of-lists.
    if (!is.null(data$profiles) && is.data.frame(data$profiles)) {
      data$profiles <- lapply(seq_len(nrow(data$profiles)), function(i) {
        as.list(data$profiles[i, , drop = FALSE])
      })
    }
    if (is.list(data$profiles) && length(data$profiles) && !is.null(names(data$profiles)) &&
        !"name" %in% names(data$profiles[[1]]) && !is.null(data$profiles$name)) {
      # Column-oriented list from JSON simplification.
      n <- length(data$profiles$name)
      data$profiles <- lapply(seq_len(n), function(i) {
        list(
          id = data$profiles$id[[i]],
          name = data$profiles$name[[i]],
          schema_version = data$profiles$schema_version[[i]] %||% 2L,
          primary_identifier = data$profiles$primary_identifier[[i]] %||% "final_name",
          tag_style = data$profiles$tag_style[[i]] %||% "bare",
          ordered_fields = unlist(data$profiles$ordered_fields[[i]], use.names = FALSE),
          required = unlist(data$profiles$required[[i]] %||% list(), use.names = FALSE)
        )
      })
    }
    # Fresh Shiny sessions adopt browser layouts. A loaded .sangerproject sets
    # persistence$from_project so project profiles keep precedence.
    current <- header_builder_config()
    precedence <- if (isTRUE(current$persistence$from_project)) "project" else "browser"
    config <- pitax_header_merge_portable_profiles(current, data, precedence = precedence)
    config$persistence$browser_saved <- TRUE
    config$persistence$last_status <- if (identical(precedence, "project")) {
      "Merged browser profiles (loaded project kept precedence on conflicts)"
    } else {
      "Restored named profiles from this browser"
    }
    header_builder_store(config)
    header_browser_restore_done(TRUE)
    if (length(config$migration_warnings)) {
      showNotification(paste(config$migration_warnings, collapse = " | "),
                       type = "warning", duration = 12)
    }
  }, ignoreInit = TRUE)

  output$export_header_preview_table <- DT::renderDT({
    req(export_dialog_state$snap)
    config <- header_builder_config()
    selection <- pitax_collect_export_selection_from_input(export_dialog_state$stage_id)
    snap <- export_dialog_state$snap
    plan <- tryCatch(pitax_export_build_plan(snap, selection), error = function(e) e)
    if (inherits(plan, "error")) {
      return(DT::datatable(data.frame(Status = conditionMessage(plan)),
                           rownames = FALSE, options = list(dom = "t")))
    }
    rows <- list()
    scope_ids <- pitax_export_filter_scope_ids(snap, plan$selection)
    for (comp in plan$sequence_components) {
      recs <- pitax_export_extractor(comp)(snap, scope_ids)
      for (rid in names(recs)) {
        rec <- recs[[rid]]
        if (isTRUE(rec$missing) || !nzchar(as.character(rec$sequence %||% ""))) next
        rendered <- tryCatch(pitax_header_render(rec, comp, rid, config), error = function(e) e)
        rows[[length(rows) + 1L]] <- data.frame(
          Component = comp,
          Record = as.character(rec$final_name %||% rid),
          Header = if (inherits(rendered, "error")) "" else paste0(">", rendered),
          Status = if (inherits(rendered, "error")) conditionMessage(rendered) else "Ready",
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(rows)) {
      return(DT::datatable(data.frame(Status = "Select a FASTA sequence component to preview."),
                           rownames = FALSE, options = list(dom = "t")))
    }
    out <- do.call(rbind, rows)
    for (component_id in unique(out$Component)) {
      indexes <- which(out$Component == component_id & out$Status == "Ready")
      if (length(indexes) < 2L) next
      ids <- tolower(sub("[[:space:]].*$", "", substring(out$Header[indexes], 2L)))
      dup <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
      out$Status[indexes[dup]] <- "Duplicate FASTA identifier"
    }
    DT::datatable(out, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = FALSE, autoWidth = TRUE,
                                 dom = "ftip",
                                 columnDefs = list(list(targets = 2, className = "dt-left"))))
  }, server = FALSE)

  # Portable JSON export/import
  output$download_fasta_header_profiles_json <- downloadHandler(
    filename = function() "PITAX_FASTA_header_profiles.json",
    content = function(file) {
      config <- header_builder_config()
      portable <- pitax_header_profiles_portable(config)
      json <- if (requireNamespace("jsonlite", quietly = TRUE)) {
        jsonlite::toJSON(portable, auto_unbox = TRUE, pretty = TRUE, null = "null")
      } else {
        # Minimal fallback encoder for the portable subset.
        paste0('{"schema":', portable$schema,
               ',"active_profile_id":"', portable$active_profile_id, '","profiles":[]}')
      }
      writeLines(json, file, useBytes = TRUE)
      config$persistence$portable_exported <- TRUE
      config$persistence$last_status <- "Portable JSON profile exported"
      header_builder_store(config)
    }
  )

  observeEvent(input$import_fasta_header_profiles_json, {
    uploaded <- input$import_fasta_header_profiles_json
    req(uploaded$datapath)
    txt <- paste(readLines(uploaded$datapath, warn = FALSE), collapse = "\n")
    ans <- tryCatch(pitax_header_import_portable_json(header_builder_config(), txt),
                    error = function(e) e)
    if (inherits(ans, "error")) {
      showNotification(conditionMessage(ans), type = "error", duration = 10)
      return()
    }
    header_builder_store(ans, mark_dirty = TRUE)
    showNotification("Imported portable header profiles (layouts only). Populate metadata values before export.",
                     type = "message", duration = 8)
  }, ignoreInit = TRUE)

  output$download_fasta_header_templates <- downloadHandler(
    filename = function() "PITAX_FASTA_header_templates.csv",
    content = function(file) {
      config <- header_builder_config()
      utils::write.csv(
        data.frame(
          Name = vapply(config$profiles, function(p) p$name, character(1)),
          Template = vapply(config$profiles, function(p) pitax_header_profile_to_template(p, config), character(1)),
          Fields = paste(unique(c(
            unlist(lapply(config$profiles, function(p) p$ordered_fields), use.names = FALSE),
            names(config$custom_fields)
          )), collapse = ";"),
          stringsAsFactors = FALSE
        ),
        file, row.names = FALSE, fileEncoding = "UTF-8"
      )
    }
  )

  observeEvent(input$import_fasta_header_templates, {
    config <- header_builder_config()
    uploaded <- input$import_fasta_header_templates
    req(uploaded$datapath)
    table <- tryCatch(utils::read.csv(uploaded$datapath, stringsAsFactors = FALSE,
                                      check.names = FALSE, fileEncoding = "UTF-8"),
                      error = function(e) e)
    if (inherits(table, "error") || !all(c("Name", "Template") %in% names(table))) {
      showNotification("Invalid template CSV; expected Name and Template columns.", type = "error")
      return()
    }
    # Convert each legacy row into a structured profile via migration helpers.
    legacy <- list(
      fields = character(),
      values = config$values,
      templates = setNames(as.list(as.character(table$Template)), trimws(as.character(table$Name))),
      active_template = trimws(as.character(table$Name[1])),
      draft_template = as.character(table$Template[1])
    )
    if ("Fields" %in% names(table)) {
      imported_fields <- unique(unlist(strsplit(paste(table$Fields, collapse = ";"), ";", fixed = TRUE)))
      legacy$fields <- trimws(as.character(imported_fields))
      legacy$fields <- legacy$fields[nzchar(legacy$fields)]
    }
    migrated <- pitax_header_migrate_legacy_config(legacy)
    # Merge into current config without dropping existing values.
    for (p in migrated$profiles) {
      if (is.na(pitax_header_find_profile_index(config, p$name))) {
        config$profiles <- c(config$profiles, list(p))
      }
    }
    for (nm in names(migrated$custom_fields)) {
      if (is.null(config$custom_fields[[nm]])) config$custom_fields[[nm]] <- migrated$custom_fields[[nm]]
    }
    config <- pitax_header_refresh_derived(config)
    header_builder_store(config, mark_dirty = TRUE)
    showNotification("Imported legacy CSV templates as structured profiles. Review before export.",
                     type = "message", duration = 8)
  }, ignoreInit = TRUE)
