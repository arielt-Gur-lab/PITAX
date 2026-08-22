  # ---------------- Rename key ----------------
  rename_key_data <- reactive({
    req(input$rename_key_file)
    f <- input$rename_key_file; ext <- tolower(tools::file_ext(f$name))
    key <- if (ext=="xlsx") readxl::read_excel(f$datapath) else if (ext=="csv") read.csv(f$datapath,stringsAsFactors=FALSE,check.names=FALSE) else stop("Assignment key must be XLSX or CSV.")
    key <- as.data.frame(key,stringsAsFactors=FALSE)
    normalized <- tolower(gsub("[^A-Za-z0-9]+", "_", names(key)))
    names(key) <- normalized
    if ("gene" %in% names(key) && !"locus" %in% names(key)) names(key)[names(key) == "gene"] <- "locus"
    required <- c("old_id", "isolate", "direction")
    if (!all(required %in% names(key)) || !any(c("assay_id", "locus") %in% names(key))) stop("Assignment key must contain old_id, isolate, direction, and assay_id or locus (gene).")
    keep <- unique(c(required, intersect(c("assay_id", "locus"), names(key))))
    key <- key[, keep, drop = FALSE]
    for (nm in names(key)) key[[nm]] <- trimws(as.character(key[[nm]]))
    key$direction <- vapply(key$direction, stage2_normalize_direction, character(1))
    key[key$old_id != "", , drop = FALSE]
  })

  output$rename_key_status <- renderUI({
    if (is.null(input$rename_key_file)) return(div(class = "compact-hint", "Key columns: old_id, isolate, direction, plus assay_id (preferred) or locus/gene."))
    div(class = "compact-hint", paste("Selected:", input$rename_key_file$name))
  })
  output$download_assignment_key_template <- downloadHandler(
    filename = function() "PITAX_assignment_key_template.csv",
    content = function(file) file.copy(file.path("templates", "assignment_key_template.csv"), file, overwrite = TRUE)
  )

  observeEvent(input$apply_rename_key, {
    req(nrow(rv$read_assignments))
    rv$read_assignments <- collect_assignment_editor()
    key <- tryCatch(rename_key_data(), error=function(e){showNotification(conditionMessage(e),type="error");NULL})
    if (is.null(key)) return()
    matched <- 0L
    for(i in seq_len(nrow(rv$read_assignments))) {
      original <- rv$read_assignments$Source_ID[i]
      idx <- which(key$old_id==original)
      if(!length(idx)) idx <- which(startsWith(original,key$old_id))
      if(length(idx)==1) {
        rv$read_assignments$Isolate[i] <- key$isolate[idx]
        assay_id <- if ("assay_id" %in% names(key)) key$assay_id[idx] else ""
        if (!nzchar(assay_id) && "locus" %in% names(key)) {
          locus_id <- pitax_normalize_locus_id(key$locus[idx])
          candidates <- rv$assay_profiles$Assay_ID[rv$assay_profiles$Locus_ID == locus_id]
          if (length(candidates) == 1L) assay_id <- candidates[1]
        }
        rv$read_assignments$Assay_ID[i] <- assay_id
        rv$read_assignments$Direction[i] <- key$direction[idx]
        matched <- matched+1L
      }
    }
    sync_assignment_state()
    showNotification(paste(matched,"of",nrow(rv$read_assignments),"reads matched; identity fields and generated names were updated."),type="message")
  })

  assignment_input_id <- function(prefix, read_id) paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", read_id))

  collect_assignment_editor <- function() {
    edited <- stage2_coerce_assignments(rv$read_assignments)
    if (!nrow(edited)) return(edited)
    for (i in seq_len(nrow(edited))) {
      read_id <- edited$Read_ID[i]
      isolate_value <- input[[assignment_input_id("assign_isolate", read_id)]]
      assay_value <- input[[assignment_input_id("assign_assay", read_id)]]
      direction_value <- input[[assignment_input_id("assign_direction", read_id)]]
      if (!is.null(isolate_value)) edited$Isolate[i] <- trimws(as.character(isolate_value)[1])
      if (!is.null(assay_value)) edited$Assay_ID[i] <- trimws(as.character(assay_value)[1])
      if (!is.null(direction_value)) edited$Direction[i] <- stage2_normalize_direction(direction_value)
    }
    edited
  }

  selected_assignment_rows <- function() {
    if (!nrow(rv$read_assignments)) return(integer())
    selected <- vapply(seq_len(nrow(rv$read_assignments)), function(i) {
      read_id <- rv$read_assignments$Read_ID[i]
      isTRUE(input[[assignment_input_id("assign_select", read_id)]])
    }, logical(1))
    rows <- which(selected)
    if (length(rows)) rows else seq_len(nrow(rv$read_assignments))
  }

  output$assignment_editor <- renderUI({
    assignments <- stage2_coerce_assignments(rv$read_assignments)
    if (!nrow(assignments)) return(div(class = "status-note", "Upload AB1 files to create the Rename table."))
    header <- tags$thead(tags$tr(
      tags$th("Use"), tags$th("Upload barcode / source"), tags$th("Final read / FASTA name"),
      tags$th("Isolate"), tags$th("Assay profile"), tags$th("Locus"), tags$th("Direction")
    ))
    body_rows <- lapply(seq_len(nrow(assignments)), function(i) {
      read_id <- assignments$Read_ID[i]
      direction <- stage2_normalize_direction(assignments$Direction[i])
      tags$tr(
        tags$td(tags$input(id = assignment_input_id("assign_select", read_id), type = "checkbox", class = "shiny-input-checkbox")),
        tags$td(class = "assignment-source", assignments$Source_ID[i]),
        tags$td(class = "assignment-final", if (nzchar(assignments$Final_Name[i])) assignments$Final_Name[i] else "- generated after Apply -"),
        tags$td(tags$input(id = assignment_input_id("assign_isolate", read_id), type = "text", class = "form-control", value = assignments$Isolate[i], autocomplete = "off")),
        tags$td(tags$select(
          id = assignment_input_id("assign_assay", read_id), class = "form-control assignment-assay",
          tags$option(value = "", disabled = NA, selected = if (!assignments$Assay_ID[i] %in% rv$assay_profiles$Assay_ID) NA, "Select..."),
          lapply(seq_len(nrow(rv$assay_profiles)), function(j) tags$option(
            value = rv$assay_profiles$Assay_ID[j],
            selected = if (assignments$Assay_ID[i] == rv$assay_profiles$Assay_ID[j]) NA,
            paste0(rv$assay_profiles$Assay_Name[j], " | ", rv$assay_profiles$Locus_Display_Name[j])
          ))
        )),
        tags$td(class = "assignment-locus", if (nzchar(assignments$Locus[i])) pitax_locus_display_name(assignments$Locus[i], assignments$Locus[i]) else "-"),
        tags$td(tags$select(
          id = assignment_input_id("assign_direction", read_id), class = "form-control assignment-direction",
          tags$option(value = "", disabled = NA, selected = if (!direction %in% c("Forward", "Reverse")) NA, "Select..."),
          tags$option(value = "Forward", selected = if (direction == "Forward") NA, "Forward"),
          tags$option(value = "Reverse", selected = if (direction == "Reverse") NA, "Reverse")
        ))
      )
    })
    div(class = "assignment-editor-wrap", tags$table(class = "assignment-editor", header, tags$tbody(body_rows)))
  })

  observeEvent(input$save_assignment_edits, {
    req(nrow(rv$read_assignments))
    rv$read_assignments <- collect_assignment_editor()
    error <- sync_assignment_state()
    if (is.null(error)) {
      showNotification("Read identity and generated FASTA names were updated.", type = "message")
    } else {
      showNotification(paste("Saved draft assignments:", error), type = "warning", duration = 8)
    }
  })

  observeEvent(input$apply_assignment_batch, {
    req(nrow(rv$read_assignments))
    rows <- selected_assignment_rows()
    rv$read_assignments <- collect_assignment_editor()
    isolate_values <- as.character(rv$read_assignments$Isolate[rows])
    if (nzchar(input$batch_isolate_find)) isolate_values <- gsub(input$batch_isolate_find, input$batch_isolate_replace, isolate_values, fixed = TRUE)
    if (nzchar(input$batch_isolate_prefix)) isolate_values <- paste0(input$batch_isolate_prefix, isolate_values)
    if (nzchar(input$batch_isolate_suffix)) isolate_values <- paste0(isolate_values, input$batch_isolate_suffix)
    rv$read_assignments$Isolate[rows] <- isolate_values
    if (!is.null(input$batch_assay) && nzchar(trimws(input$batch_assay))) rv$read_assignments$Assay_ID[rows] <- trimws(input$batch_assay)
    if (nzchar(input$batch_direction)) rv$read_assignments$Direction[rows] <- input$batch_direction
    sync_assignment_state()
    showNotification(paste("Updated", length(rows), "read assignment(s)."), type = "message")
  })

  observeEvent(input$reset_assignments, {
    req(nrow(rv$read_assignments))
    rows <- selected_assignment_rows()
    rv$read_assignments <- collect_assignment_editor()
    rv$read_assignments$Isolate[rows] <- ""
    rv$read_assignments$Assay_ID[rows] <- ""
    rv$read_assignments$Locus[rows] <- ""
    rv$read_assignments$Direction[rows] <- "Unknown"
    rv$read_assignments$Primer[rows] <- ""
    sync_assignment_state()
  })

  output$batch_assay_control <- renderUI({
    profiles <- assay_coerce_profiles(rv$assay_profiles)
    choices <- c("No change" = "", stats::setNames(profiles$Assay_ID, paste0(profiles$Assay_Name, " | ", profiles$Locus_Display_Name)))
    selectInput("batch_assay", "Set assay profile", choices = choices, selected = "")
  })

  rename_error <- reactive({
    error <- stage2_identity_error(rv$read_assignments)
    if (is.null(error)) error <- stage2_validate_assignments(rv$read_assignments, assay_profiles = rv$assay_profiles)
    error
  })
  assignment_validation_ui <- function() {
    e <- stage2_identity_error(rv$read_assignments)
    if (is.null(e)) e <- stage2_validate_assignments(rv$read_assignments, assay_profiles = rv$assay_profiles)
    if(is.null(e)) div(class="status-ok","OK: Explicit identity fields are complete; final read names were generated by PITAX.") else div(class="status-error",paste0("Warning: ",e))
  }
  output$rename_validation <- renderUI(assignment_validation_ui())
