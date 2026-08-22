  # ---------------- Stage 3 consensus ----------------
  output$consensus_mode_note <- renderUI({
    if (identical(rv$project_mode, "simple")) {
      div(class = "status-note", strong("Project mode: Simple reads. "), "Each read remains an independent analysis sequence. Reverse reads are reverse-complemented, but PITAX does not search for or merge a matching Forward read.")
    } else {
      div(class = "status-note", strong("Project mode: Paired Forward/Reverse. "), "Reads sharing the same isolate and gene are paired only when one is Forward and one is Reverse. A lone read remains a valid single-read representative.")
    }
  })

  sync_consensus_choices <- function(preferred = NULL) {
    ids <- if (is.list(rv$consensus_set) && is.list(rv$consensus_set$records)) names(rv$consensus_set$records) else character()
    if (!length(ids)) {
      updateSelectInput(session, "consensus_sample", choices = character(), selected = character())
      return(invisible(NULL))
    }
    labels <- vapply(rv$consensus_set$records, function(x) stage3_scalar_text(x$final_name), character(1))
    choices <- setNames(ids, labels)
    selected <- if (!is.null(preferred) && preferred %in% ids) preferred else if (!is.null(input$consensus_sample) && input$consensus_sample %in% ids) input$consensus_sample else ids[1]
    updateSelectInput(session, "consensus_sample", choices = choices, selected = selected)
    invisible(selected)
  }

  build_analysis_sequences <- function(notify = TRUE) {
    if (!length(rv$results)) {
      if (notify) showNotification("No processed reads are available.", type = "error")
      return(FALSE)
    }
    locus_error <- stage3_run_locus_error(rv$read_assignments)
    if (!is.null(locus_error)) {
      if (notify) showNotification(locus_error, type = "error", duration = 10)
      return(FALSE)
    }
    built <- tryCatch(
      stage3_build_consensus_set(
        rv$read_assignments, rv$results,
        min_overlap = input$consensus_min_overlap,
        min_identity = input$consensus_min_identity,
        quality_delta = input$consensus_quality_delta,
        strong_quality = input$consensus_strong_quality,
        project_mode = rv$project_mode
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "consensus_build_error")
    )
    if (inherits(built, "consensus_build_error")) {
      if (notify) showNotification(paste("Consensus build failed:", built$error), type = "error", duration = 12)
      return(FALSE)
    }
    old_records <- tryCatch(stage3_analysis_records(rv$consensus_set), error = function(e) list())
    new_records <- stage3_analysis_records(built)
    old_signature <- if (length(old_records)) vapply(old_records, function(x) paste0(x$final_name, "\r", x$seq), character(1)) else character()
    new_signature <- if (length(new_records)) vapply(new_records, function(x) paste0(x$final_name, "\r", x$seq), character(1)) else character()
    if (!identical(old_signature, new_signature) && (nrow(rv$blast_jobs) || nrow(rv$blast_hits) || nrow(rv$taxonomy_summary))) {
      if (nrow(rv$blast_jobs)) rv$blast_jobs$status <- "STALE"
      rv$blast_raw <- list(); rv$blast_hits <- data.frame(); rv$blast_ids <- data.frame()
      rv$taxonomy_summary <- data.frame(); rv$taxonomy_hits <- data.frame(); rv$taxonomy_counts <- data.frame()
      rv$blast_batch_status_text <- "Stage 3 analysis sequences changed; previous BLAST evidence is stale and must be rerun."
      rv$taxonomy_status_text <- "Stage 3 analysis sequences changed; previous taxonomic interpretation was removed. Re-run BLAST, then Taxonomy."
    }
    rv$consensus_set <- built
    sync_consensus_choices()
    rv$project_status_text <- paste0("Unsaved Stage 3 build: ", nrow(built$summary), " analysis sequence(s).")
    gate_error <- stage3_consensus_gate_error(rv$consensus_set, rv$results)
    if (notify) {
      showNotification(
        if (is.null(gate_error)) "Stage 3 sequence gate is green." else gate_error,
        type = if (is.null(gate_error)) "message" else "warning", duration = 10
      )
    }
    is.null(gate_error)
  }

  observeEvent(input$to_consensus, {
    if (!length(rv$results)) {
      showNotification("Run trimming and review QC before building analysis sequences.", type = "error")
      return()
    }
    if (identical(rv$project_mode, "simple")) {
      if (isTRUE(build_analysis_sequences(notify = TRUE))) {
        updateTabsetPanel(session, "pipeline_step", selected = "export")
      }
      return()
    }
    locus_error <- stage3_run_locus_error(rv$read_assignments)
    if (!is.null(locus_error)) {
      showNotification(locus_error, type = "error", duration = 10)
      return()
    }
    updateTabsetPanel(session, "pipeline_step", selected = "consensus")
  })

  observeEvent(input$build_consensus, {
    build_analysis_sequences(notify = TRUE)
  })

  output$consensus_gate_status <- renderUI({
    error <- stage3_consensus_gate_error(rv$consensus_set, rv$results)
    if (is.null(error)) {
      div(class = "status-ok", "✓ Every expected analysis sequence is current and review-complete.")
    } else {
      div(class = "status-warning", paste0("⚠ ", error))
    }
  })

  output$consensus_summary_table <- renderDT({
    df <- if (is.list(rv$consensus_set) && is.data.frame(rv$consensus_set$summary)) rv$consensus_set$summary else stage3_empty_summary()
    if (!nrow(df)) return(datatable(data.frame(Message = "Build analysis sequences to populate this table."), rownames = FALSE, options = list(dom = "t")))
    visible <- c("Final_Name", "Isolate", "Locus", "Status", "Source_Reads", "Length", "Overlap", "Identity_percent", "Review_positions", "Revision")
    df <- df[, intersect(visible, names(df)), drop = FALSE]
    datatable(
      df, rownames = FALSE, selection = "single", class = "compact stripe hover consensus-summary-stable",
      options = list(
        pageLength = 20, scrollX = TRUE, autoWidth = TRUE, deferRender = TRUE, dom = "tip",
        columnDefs = list(
          list(targets = 0, width = "190px"), list(targets = 1, width = "115px"),
          list(targets = 2, width = "85px"), list(targets = 3, width = "150px"),
          list(targets = 4:9, width = "82px", className = "dt-center")
        )
      )
    )
  })

  observeEvent(input$consensus_summary_table_rows_selected, {
    row <- input$consensus_summary_table_rows_selected
    if (length(row) == 1L && row >= 1L && row <= nrow(rv$consensus_set$summary)) {
      updateSelectInput(session, "consensus_sample", selected = rv$consensus_set$summary$Consensus_ID[row])
    }
  })

  selected_consensus <- reactive({
    req(is.list(rv$consensus_set), length(rv$consensus_set$records), input$consensus_sample)
    record <- rv$consensus_set$records[[as.character(input$consensus_sample)]]
    req(is.list(record))
    record
  })

  sync_consensus_review_choices <- function(preferred = NULL) {
    id <- stage3_scalar_text(input$consensus_sample)
    record <- if (nzchar(id) && is.list(rv$consensus_set$records[[id]])) rv$consensus_set$records[[id]] else NULL
    review <- if (is.list(record)) stage3_reviewable_evidence(record) else data.frame()
    if (!nrow(review)) {
      updateSelectInput(session, "consensus_review_position", choices = character(), selected = character())
      return(invisible(NULL))
    }
    values <- as.character(review$Alignment_Column)
    labels <- paste0(
      "Position ", review$Consensus_Position, " · F:", review$Forward_Base,
      " / R:", review$Reverse_Base, ifelse(review$Needs_Review, " · unresolved", " · reviewed")
    )
    selected <- if (!is.null(preferred) && as.character(preferred) %in% values) as.character(preferred) else values[1]
    updateSelectInput(session, "consensus_review_position", choices = setNames(values, labels), selected = selected)
    invisible(selected)
  }

  observeEvent(input$consensus_sample, sync_consensus_review_choices(), ignoreInit = FALSE)

  selected_consensus_review_row <- reactive({
    record <- selected_consensus()
    column <- suppressWarnings(as.integer(input$consensus_review_position))
    req(is.finite(column), is.data.frame(record$evidence), nrow(record$evidence))
    idx <- match(column, record$evidence$Alignment_Column)
    req(!is.na(idx))
    record$evidence[idx, , drop = FALSE]
  })

  output$consensus_review_call_ui <- renderUI({
    row <- selected_consensus_review_row()
    candidate_values <- c(as.character(row$Forward_Base), as.character(row$Reverse_Base), as.character(row$Automatic_Call))
    candidate_labels <- c(paste0("Forward · ", row$Forward_Base), paste0("Reverse · ", row$Reverse_Base), paste0("IUPAC / automatic · ", row$Automatic_Call))
    candidates <- setNames(candidate_values, candidate_labels)
    candidates <- candidates[candidates != "-" & nzchar(candidates)]
    candidates <- candidates[!duplicated(unname(candidates))]
    radioButtons("consensus_review_call", "Accepted consensus call", choices = candidates,
                 selected = as.character(row$Consensus_Call), inline = FALSE)
  })

  output$consensus_review_state <- renderUI({
    record <- selected_consensus()
    record <- stage3_ensure_record_curation(record)
    reviewable <- stage3_reviewable_evidence(record)
    unresolved <- if (nrow(record$evidence)) sum(record$evidence$Needs_Review, na.rm = TRUE) else 0L
    div(class = if (unresolved) "status-warning" else "status-ok",
        paste0("Revision ", record$curation$revision, " · ", unresolved, " unresolved · ", nrow(reviewable), " reviewable position(s)."))
  })

  consensus_review_plot <- function(direction) {
    record <- selected_consensus()
    row <- selected_consensus_review_row()
    source_id <- if (identical(direction, "Forward")) stage3_scalar_text(record$forward_read) else stage3_scalar_text(record$reverse_read)
    raw_position <- if (identical(direction, "Forward")) row$Forward_Raw_Position else row$Reverse_Raw_Position
    result <- if (nzchar(source_id)) rv$results[[source_id]] else NULL
    title <- paste0(direction, " evidence", if (nzchar(source_id)) paste0(" · ", source_id) else "")
    make_chromatogram_focus_plot(result, settings_for_result(result), raw_position, title)
  }

  output$consensus_forward_conflict_plot <- plotly::renderPlotly(consensus_review_plot("Forward"))
  output$consensus_reverse_conflict_plot <- plotly::renderPlotly(consensus_review_plot("Reverse"))

  commit_consensus_record <- function(id, record, previous_sequence, action_label) {
    rv$consensus_set$records[[id]] <- record
    rv$consensus_set <- stage3_refresh_consensus_summary(rv$consensus_set)
    sequence_changed <- !identical(previous_sequence, stage3_scalar_text(record$sequence))
    invalidate_downstream_for_sample(
      id,
      paste0("Consensus ", action_label, if (sequence_changed) " changed the sequence" else " changed the review state", "; revision ", record$curation$revision)
    )
    rv$project_status_text <- paste0("Unsaved consensus ", action_label, ": ", record$final_name, " · revision ", record$curation$revision, ".")
    sync_consensus_review_choices(input$consensus_review_position)
  }

  observeEvent(input$apply_consensus_review, {
    id <- stage3_scalar_text(input$consensus_sample)
    req(nzchar(id), input$consensus_review_position, input$consensus_review_call)
    old <- rv$consensus_set$records[[id]]
    previous_sequence <- stage3_scalar_text(old$sequence)
    updated <- tryCatch(
      stage3_apply_consensus_call(old, input$consensus_review_position, input$consensus_review_call,
                                  method = "Forward / Reverse / IUPAC decision", note = input$consensus_review_note),
      error = function(e) e
    )
    if (inherits(updated, "error")) {
      showNotification(conditionMessage(updated), type = "error", duration = 10)
      return()
    }
    commit_consensus_record(id, updated, previous_sequence, "review")
    updateTextInput(session, "consensus_review_note", value = "")
    showNotification("Consensus decision saved with a new revision.", type = "message")
  })

  observeEvent(input$consensus_review_undo, {
    id <- stage3_scalar_text(input$consensus_sample); req(nzchar(id))
    old <- rv$consensus_set$records[[id]]; previous_sequence <- stage3_scalar_text(old$sequence)
    updated <- tryCatch(stage3_consensus_undo(old), error = function(e) e)
    if (inherits(updated, "error")) { showNotification(conditionMessage(updated), type = "warning"); return() }
    commit_consensus_record(id, updated, previous_sequence, "undo")
  })

  observeEvent(input$consensus_review_redo, {
    id <- stage3_scalar_text(input$consensus_sample); req(nzchar(id))
    old <- rv$consensus_set$records[[id]]; previous_sequence <- stage3_scalar_text(old$sequence)
    updated <- tryCatch(stage3_consensus_redo(old), error = function(e) e)
    if (inherits(updated, "error")) { showNotification(conditionMessage(updated), type = "warning"); return() }
    commit_consensus_record(id, updated, previous_sequence, "redo")
  })

  observeEvent(input$consensus_review_history, {
    record <- stage3_ensure_record_curation(selected_consensus())
    output$consensus_review_history_table <- renderDT({
      datatable(record$curation$audit_log, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE, autoWidth = TRUE, order = list(list(0, "desc"))))
    })
    showModal(modalDialog(title = paste0("Consensus revision history · ", record$final_name),
                          DTOutput("consensus_review_history_table"), size = "l", easyClose = TRUE, footer = modalButton("Close")))
  })

  output$consensus_selected_status <- renderUI({
    record <- selected_consensus()
    class_name <- if (record$status %in% c("READY", "SINGLE_READ", "INDEPENDENT_READ")) "status-ok" else if (record$status == "REVIEW_REQUIRED") "status-warning" else "status-error"
    label <- switch(record$status,
      READY = "Ready · Forward and Reverse overlap passed",
      SINGLE_READ = "Valid single-read representative · not a two-read consensus",
      INDEPENDENT_READ = "Ready · independent read, oriented without Forward/Reverse matching",
      REVIEW_REQUIRED = "Review required · unresolved conflicts are retained as IUPAC calls",
      NO_RELIABLE_OVERLAP = "Blocked · no reliable Forward/Reverse overlap",
      SOURCE_MISSING = "Blocked · a processed source read is missing",
      record$status
    )
    div(class = class_name, label)
  })

  output$consensus_selected_metrics <- renderDT({
    record <- selected_consensus()
    df <- data.frame(
      Metric = c("Isolate", "Locus", "Forward read", "Reverse read", "Sequence length", "Overlap", "Overlap identity", "Mismatches", "Indels", "Review positions"),
      Value = c(record$isolate, record$locus, stage3_scalar_text(record$forward_read, "—"), stage3_scalar_text(record$reverse_read, "—"),
                nchar(stage3_scalar_text(record$sequence)), record$metrics$overlap,
                if (is.finite(record$metrics$identity_percent)) paste0(record$metrics$identity_percent, "%") else "—",
                record$metrics$mismatches, record$metrics$indels, record$metrics$review_positions),
      stringsAsFactors = FALSE
    )
    datatable(df, rownames = FALSE, options = list(dom = "t", ordering = FALSE))
  })

  output$consensus_alignment_text <- renderText({
    record <- selected_consensus()
    if (is.null(record$alignment)) return("Single read: no pairwise overlap alignment is required.")
    f <- record$alignment$forward_aligned
    r <- record$alignment$reverse_aligned
    n <- nchar(f)
    starts <- seq.int(1L, max(1L, n), by = 80L)
    blocks <- vapply(starts, function(st) {
      en <- min(n, st + 79L)
      fs <- substr(f, st, en); rs <- substr(r, st, en)
      fc <- strsplit(fs, "", fixed = TRUE)[[1]]; rc <- strsplit(rs, "", fixed = TRUE)[[1]]
      middle <- paste0(ifelse(fc == rc & fc != "-", "|", ifelse(fc == "-" | rc == "-", " ", ".")), collapse = "")
      paste0("F  ", fs, "\n   ", middle, "\nR  ", rs, "\n")
    }, character(1))
    paste(blocks, collapse = "\n")
  })

  output$consensus_sequence_text <- renderText({
    record <- selected_consensus()
    sequence <- stage3_scalar_text(record$sequence)
    if (!nzchar(sequence)) return("No downstream sequence is emitted while this record is blocked.")
    wrap_sequence(sequence, 80)
  })

  output$consensus_evidence_table <- renderDT({
    record <- selected_consensus()
    evidence <- record$evidence
    if (!is.data.frame(evidence) || !nrow(evidence)) return(datatable(data.frame(Message = "No per-column evidence is available."), rownames = FALSE, options = list(dom = "t")))
    datatable(evidence, rownames = FALSE, filter = "top", options = list(pageLength = 25, scrollX = TRUE, autoWidth = TRUE))
  })

  output$download_consensus_fasta <- downloadHandler(
    filename = function() paste0(project_export_stem(), "_isolate_level.fasta"),
    content = function(file) {
      req(is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
      writeLines(make_fasta(stage3_analysis_records(rv$consensus_set), FALSE), file)
    }
  )

  output$download_consensus_checkpoint <- downloadHandler(
    filename = function() paste0(project_export_stem(), "_checkpoint_C_consensus.zip"),
    content = function(file) {
      req(is.list(rv$consensus_set), length(rv$consensus_set$records))
      write_consensus_checkpoint_zip(file, rv$consensus_set, rv$read_assignments, rv$settings)
    }
  )

