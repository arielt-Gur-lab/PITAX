  # ---------------- PITAX 3.0 Stage 1: AB1 evidence audit ----------------
  output$ab1_evidence_run_table <- renderDT({
    req(rv$results)
    df <- ab1_evidence_run_summary(rv$results)
    if (!nrow(df)) {
      return(datatable(data.frame(Message = "No AB1 evidence has been captured in this run."),
                       rownames = FALSE, selection = "none", options = list(dom = "t")))
    }
    display <- df
    display$Sample <- vapply(display$Sample, qc_display_name, character(1))
    names(display) <- c(
      "Sample", "Evidence", "Quality tag", "Quality coverage (%)",
      "Auto trim start", "Auto trim end", "Auto trim length", "Median quality | auto trim",
      "Q>=20 | auto trim (%)", "Q>=30 | auto trim (%)",
      "PCON window start", "PCON window end", "PCON window length", "PCON window shift",
      "Median quality | PCON window", "Q>=20 | PCON window (%)", "Q>=30 | PCON window (%)",
      "Primary position source",
      "Primary position coverage (%)", "Primary positions different (%)", "Median |position Delta|",
      "Legacy map", "Canonical A/C/G/T map", "Maps match",
      "Legacy call is max (%)", "Canonical call is max (%)",
      "Canonical call is max | auto trim (%)", "Median called/alternative ratio | auto trim"
    )
    datatable(
      display,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      options = list(pageLength = 12, scrollX = TRUE, autoWidth = TRUE, dom = "tip")
    )
  })

  # The run-audit table and the Sample dropdown now share one source of truth.
  # Clicking an audit row selects the same sample throughout the QC workspace.
  observeEvent(input$ab1_evidence_run_table_rows_selected, {
    idx <- input$ab1_evidence_run_table_rows_selected
    if (is.null(idx) || length(idx) != 1L || idx < 1L) return()
    df <- ab1_evidence_run_summary(rv$results)
    if (idx > nrow(df)) return()
    key <- pitax_result_key_for_sample(rv$results, df$Sample[idx])
    if (nzchar(key)) updateSelectInput(session, "inspect_sample", selected = key)
  }, ignoreInit = TRUE)

  output$ab1_evidence_selected_note <- renderUI({
    key <- selected_sample_key()
    r <- selected_result()
    sid <- pitax_result_sample_id(r, key)
    display_name <- qc_display_name(key)
    sample_label <- if (!identical(display_name, sid)) paste0(display_name, " (original: ", sid, ")") else sid
    ev <- r$ab1_evidence
    if (is.null(ev)) {
      return(div(class = "status-note",
                 tags$strong(paste0("Selected sample: ", sample_label, ". ")),
                 "This sample has no Stage 1 evidence object. Reprocess the original AB1 with the current build to create the audit."))
    }
    if (!is.null(ev$error) && nzchar(as.character(ev$error)[1])) {
      return(div(class = "status-error", tags$strong(paste0("Selected sample: ", sample_label, ". ")),
                 paste0("Evidence audit error: ", as.character(ev$error)[1],
                        ". The established trimming result was preserved.")))
    }
    sm <- ab1_evidence_result_summary(r)
    tagList(
      div(class = "compact-hint", tags$strong(paste0("Selected sample: ", sample_label))),
      div(class = "peak-flag-summary",
          div(class = "peak-flag-pill", tags$strong(sm$Quality_tag[1]), " quality tag"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Median_quality_auto_trim[1]), sm$Median_quality_auto_trim[1], "NA")), " median quality in auto trim"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Q20_auto_trim_percent[1]), paste0(sm$Q20_auto_trim_percent[1], "%"), "NA")), " Q>=20 in auto trim"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Primary_position_difference_percent[1]), paste0(sm$Primary_position_difference_percent[1], "%"), "NA")), " PLOC vs legacy positions differ")),
      div(class = "compact-hint",
          "The current build retains the active legacy auto trim and a same-length PCON-only comparison. The comparison is observational and does not change the processed sequence.")
    )
  })

  output$ab1_trim_comparison_table <- renderDT({
    d <- pitax_trim_window_comparison(selected_result())
    names(d) <- c(
      "Window", "Status", "Start", "End", "Length", "Quality coverage (%)",
      "Median quality", "Q>=20 (%)", "Q>=30 (%)", "Start shift"
    )
    datatable(d, rownames = FALSE, selection = "none", options = list(dom = "t", scrollX = TRUE, autoWidth = TRUE))
  })

  output$ab1_evidence_detail_table <- renderDT({
    key <- selected_sample_key()
    r <- selected_result()
    ev <- r$ab1_evidence
    if (is.null(ev) || !is.data.frame(ev$detail) || !nrow(ev$detail)) {
      return(datatable(data.frame(Message = "No per-base Stage 1 evidence available for this sample."),
                       rownames = FALSE, selection = "none", options = list(dom = "t")))
    }
    d <- pitax_add_trim_membership(ev$detail, r)
    d$Sample_ID <- pitax_result_sample_id(r, key)
    d$Final_Name <- qc_display_name(key)
    keep <- c(
      "Sample_ID", "Final_Name", "Position", "Base", "In_auto_trim", "In_quality_proposed_window", "Basecaller_quality",
      "Legacy_primary_peak_pos", "Raw_ABIF_primary_peak_pos", "Primary_peak_pos_delta",
      "Legacy_called_signal", "Legacy_called_is_max",
      "Canonical_called_signal", "Canonical_best_alt_signal", "Canonical_called_to_alt_ratio",
      "Canonical_best_channel", "Canonical_called_is_max"
    )
    d <- d[, intersect(keep, names(d)), drop = FALSE]
    numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    for (nm in numeric_cols) d[[nm]] <- round(d[[nm]], 3)
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE, dom = "tip"))
  })

  output$download_ab1_evidence_run <- downloadHandler(
    filename = function() paste0("PITAX_v3_stage1_AB1_run_audit_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      df <- ab1_evidence_run_summary(rv$results)
      df$Final_Name <- vapply(df$Sample, qc_display_name, character(1))
      df <- df[, c("Sample", "Final_Name", setdiff(names(df), c("Sample", "Final_Name"))), drop = FALSE]
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )

  output$download_ab1_evidence_detail <- downloadHandler(
    filename = function() {
      key <- selected_sample_key()
      r <- selected_result()
      pitax_assert_export_identity(r, key)
      paste0(clean_fasta_name(qc_display_name(key)), "_PITAX_v3_stage1_AB1_base_audit.csv")
    },
    content = function(file) {
      key <- selected_sample_key()
      r <- selected_result()
      d <- pitax_evidence_detail_export(r, key)
      d$Final_Name <- qc_display_name(key)
      d <- d[, c("Sample_ID", "Final_Name", setdiff(names(d), c("Sample_ID", "Final_Name"))), drop = FALSE]
      utils::write.csv(d, file, row.names = FALSE, na = "")
    }
  )
  output$primer_match_table <- renderDT({
    req(isTRUE(input$enable_primer_mapping))
    datatable(primer_match_table(selected_result(), selected_processing_settings()), rownames = FALSE,
              options = list(dom = "t", scrollX = TRUE, autoWidth = TRUE))
  })
  output$primer_alignment <- renderText({
    req(isTRUE(input$enable_primer_mapping))
    primer_alignment_text(selected_result(), selected_processing_settings())
  })
  output$amplicon_overview <- renderPlot({
    draw_amplicon_overview(selected_result(), selected_processing_settings())
  })
  output$expected_amplicon_note <- renderUI({
    settings <- selected_processing_settings()
    span(class = "peak-flag-pill", title = "Assay metadata only; not a fixed coordinate on the Sanger read.",
         paste0("Expected amplicon: ", settings$expected_amplicon_len, " bp | Locus: ", settings$target))
  })

  all_current_peak_flags <- reactive({
    req(selected_result())
    scope <- if (!is.null(input$peak_flag_scope) && input$peak_flag_scope %in% c("trimmed", "raw")) input$peak_flag_scope else "trimmed"
    df <- ambiguous_peak_flags(selected_result(), scope = scope, params = ambiguous_peak_params_from_settings(selected_processing_settings()))
    if (!nrow(df)) return(df)
    r <- ensure_curation_state(selected_result())
    reviewed <- as.integer(r$curation$reviewed_positions)
    df$Review_status <- ifelse(df$Position %in% reviewed, "Reviewed", "Active")
    df
  })

  current_peak_flags <- reactive({
    df <- all_current_peak_flags()
    if (!nrow(df)) return(df)
    if (!isTRUE(input$show_reviewed_flags)) df <- df[df$Review_status != "Reviewed", , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  center_chromatogram_at <- function(pos, half_window = 20) {
    n <- length(selected_result()$peak_pos)
    pos <- suppressWarnings(as.numeric(pos))
    if (!is.finite(pos) || n < 1) return(invisible(NULL))
    range <- c(max(1, pos - half_window), min(n, pos + half_window))
    proxy <- plotly::plotlyProxy("chromatogram_plot", session)
    plotly::plotlyProxyInvoke(proxy, "relayout", list(xaxis.range = range))
    invisible(NULL)
  }

  output$ambiguous_peak_summary <- renderUI({
    all_df <- all_current_peak_flags()
    active <- if (nrow(all_df)) all_df[all_df$Review_status == "Active", , drop = FALSE] else all_df
    total <- nrow(active)
    strong_count <- if (total) sum(active$Severity == "Strong", na.rm = TRUE) else 0L
    moderate_count <- if (total) sum(active$Severity == "Moderate", na.rm = TRUE) else 0L
    auto_count <- if (total && "Auto_correct_candidate" %in% names(active)) sum(active$Auto_correct_candidate %in% TRUE, na.rm = TRUE) else 0L
    reviewed_count <- if (nrow(all_df)) sum(all_df$Review_status == "Reviewed", na.rm = TRUE) else 0L
    r <- ensure_curation_state(selected_result())
    div(class = "peak-flag-summary",
        div(class = "peak-flag-pill", tags$strong(total), " active"),
        div(class = "peak-flag-pill", tags$strong(strong_count), " strong"),
        div(class = "peak-flag-pill", tags$strong(moderate_count), " moderate"),
        if (auto_count > 0) div(class = "peak-flag-pill", tags$strong(auto_count), " high-confidence proposals") else NULL,
        if (reviewed_count > 0) div(class = "peak-flag-pill", tags$strong(reviewed_count), " reviewed") else NULL,
        if (r$curation$revision > 0) div(class = "manual-edit-badge", paste0("Curation revision ", r$curation$revision)) else NULL)
  })

  output$ambiguous_peak_table <- renderDT({
    df <- current_peak_flags()
    if (!nrow(df)) {
      return(datatable(
        data.frame(Message = "No active ambiguous channel-competition positions in the selected scope."),
        rownames = FALSE, selection = "none", options = list(dom = "t")
      ))
    }
    show <- df[, c("Position", "Call", "Competing_channel", "Peak_ratio", "Competitor_percent",
                   "Competitor_peak_offset", "Severity", "Flag", "Auto_correct_candidate", "Review_status"), drop = FALSE]
    show$Auto_correct_candidate <- ifelse(show$Auto_correct_candidate %in% TRUE, "Yes", "")
    names(show) <- c("Position", "Call", "Competing channel", "Peak ratio", "Competitor / called (%)",
                     "Peak offset", "Severity", "Flag", "Auto", "Status")
    datatable(
      show,
      rownames = FALSE,
      selection = list(mode = "single", selected = NULL, target = "row"),
      filter = "none",
      callback = JS("
        table.on('contextmenu', 'tbody tr', function(e) {
          e.preventDefault();
          var row = table.row(this).data();
          if (!row || row.length < 1) return;
          Shiny.setInputValue('peak_flag_context', {position: row[0], nonce: Math.random()}, {priority: 'event'});
        });
      "),
      options = list(pageLength = 10, lengthChange = FALSE, scrollX = TRUE, autoWidth = TRUE, dom = "tip", order = list(list(0, "asc")))
    )
  })

  observeEvent(input$ambiguous_peak_table_rows_selected, {
    idx <- input$ambiguous_peak_table_rows_selected
    df <- current_peak_flags()
    if (!length(idx) || !nrow(df) || idx[1] < 1 || idx[1] > nrow(df)) return()
    center_chromatogram_at(df$Position[idx[1]])
  }, ignoreInit = TRUE)

  flag_evidence_text <- function(flag_row) {
    if (is.null(flag_row) || !nrow(flag_row)) return("")
    paste0(
      "Call ", flag_row$Call[1], "; competitor ", flag_row$Competing_channel[1],
      "; called signal ", flag_row$Called_signal[1],
      "; competitor signal ", flag_row$Competitor_signal[1],
      "; peak ratio ", flag_row$Peak_ratio[1],
      "; competitor/called ", flag_row$Competitor_percent[1], "%",
      if ("Competitor_peak_offset" %in% names(flag_row) && is.finite(flag_row$Competitor_peak_offset[1])) paste0("; peak offset ", flag_row$Competitor_peak_offset[1]) else "",
      "; flag: ", flag_row$Flag[1]
    )
  }

  show_peak_curation_menu <- function(flag_row) {
    if (is.null(flag_row) || !nrow(flag_row)) return()
    rv$context_peak_flag <- flag_row[1, , drop = FALSE]
    pos <- as.integer(flag_row$Position[1])
    call <- as.character(flag_row$Call[1])
    comp <- as.character(flag_row$Competing_channel[1])
    r <- ensure_curation_state(selected_result())
    st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
    in_trim <- is.finite(st) && is.finite(en) && pos >= st && pos <= en
    trim_len <- if (in_trim) en - st + 1L else NA_integer_
    edge_zone <- if (in_trim) max(8L, min(40L, ceiling(trim_len * 0.15))) else NA_integer_
    recommendation <- ""
    if (in_trim) {
      dl <- pos - st; dr <- en - pos
      if (dl <= edge_zone || dr <= edge_zone) recommendation <- if (dl <= dr) "Recommended: trim the left edge through this position." else "Recommended: trim the right edge from this position."
    }
    comp_candidates <- unique(strsplit(gsub("[^ACGT/]", "", comp), "/", fixed = TRUE)[[1]])
    comp_candidates <- comp_candidates[comp_candidates %in% c("A","C","G","T")]
    base_choices <- unique(c(comp_candidates, c("A","C","G","T","N")))
    base_choices <- setNames(base_choices, base_choices)
    pair_code <- if (call %in% c("A","C","G","T") && length(comp_candidates) && comp_candidates[1] != call) iupac_for_pair(call, comp_candidates[1]) else NA_character_

    showModal(modalDialog(
      title = paste0("Review chromatogram position ", pos),
      easyClose = TRUE, size = "m",
      div(class = "about-callout",
          tags$strong(paste0(call, "  <->  ", comp)), tags$br(),
          span(flag_evidence_text(flag_row))),
      if (nzchar(recommendation)) div(class = "status-note", recommendation) else NULL,
      selectInput("manual_base_choice", "Change current base to", choices = base_choices, selected = if (length(comp_candidates)) comp_candidates[1] else "N"),
      div(class = "curation-actions",
          actionButton("context_change_base_preview", "Review base change", class = "btn-primary"),
          if (!is.na(pair_code)) actionButton("context_iupac_preview", paste0("Set ambiguity ", pair_code)) else NULL,
          actionButton("context_mark_reviewed", "Keep call / mark reviewed")),
      if (in_trim) tagList(
        tags$hr(),
        div(class = "curation-actions",
            actionButton("context_trim_left", paste0("Trim left through ", pos)),
            actionButton("context_trim_right", paste0("Trim right from ", pos)))
      ) else NULL,
      footer = modalButton("Close")
    ))
  }

  observeEvent(input$peak_flag_context, {
    pos <- suppressWarnings(as.integer(input$peak_flag_context$position))
    df <- all_current_peak_flags()
    row <- df[df$Position == pos, , drop = FALSE]
    if (!nrow(row)) return()
    center_chromatogram_at(pos)
    show_peak_curation_menu(row[1, , drop = FALSE])
  }, ignoreInit = TRUE)

  show_curation_confirmation <- function(type, position, value = NULL) {
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(position)
    flag_row <- rv$context_peak_flag
    evidence <- flag_evidence_text(flag_row)
    if (type == "base") {
      calls <- curated_raw_calls(r); before <- if (pos >= 1 && pos <= length(calls)) calls[pos] else ""
      after <- toupper(as.character(value))
      title <- "Confirm base edit"
      body <- tagList(
        h4(paste0("Position ", pos, ": ", before, " -> ", after)),
        p(evidence),
        div(class = "status-note", "The raw AB1 trace and original base call are preserved. This edit changes only the curated sequence and will be recorded in the audit log."))
    } else {
      st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
      side <- as.character(value)
      new_len <- if (side == "left") en - pos else pos - st
      title <- "Confirm manual trimming"
      body <- tagList(
        h4(if (side == "left") paste0("Remove positions ", st, "-", pos) else paste0("Remove positions ", pos, "-", en)),
        p(paste0("Current length: ", en - st + 1L, " bp | New length: ", max(0, new_len), " bp")),
        p(evidence),
        div(class = "status-note", "The automatic trimming boundaries remain stored separately. This manual boundary change is reversible and logged."))
    }
    rv$pending_curation <- list(type = type, position = pos, value = value, evidence = evidence)
    showModal(modalDialog(
      title = title, body,
      footer = tagList(modalButton("Cancel"), actionButton("confirm_curation_action", "Confirm change", class = "btn-danger")),
      easyClose = FALSE
    ))
  }

  observeEvent(input$context_change_base_preview, {
    req(rv$context_peak_flag, input$manual_base_choice)
    show_curation_confirmation("base", rv$context_peak_flag$Position[1], input$manual_base_choice)
  }, ignoreInit = TRUE)

  observeEvent(input$context_iupac_preview, {
    req(rv$context_peak_flag)
    call <- as.character(rv$context_peak_flag$Call[1])
    comp <- strsplit(as.character(rv$context_peak_flag$Competing_channel[1]), "/", fixed = TRUE)[[1]][1]
    code <- iupac_for_pair(call, comp)
    if (is.na(code)) { showNotification("No two-base IUPAC ambiguity code is available for this pair.", type = "warning"); return() }
    show_curation_confirmation("base", rv$context_peak_flag$Position[1], code)
  }, ignoreInit = TRUE)

  observeEvent(input$context_trim_left, {
    req(rv$context_peak_flag)
    show_curation_confirmation("trim", rv$context_peak_flag$Position[1], "left")
  }, ignoreInit = TRUE)

  observeEvent(input$context_trim_right, {
    req(rv$context_peak_flag)
    show_curation_confirmation("trim", rv$context_peak_flag$Position[1], "right")
  }, ignoreInit = TRUE)

  observeEvent(input$context_mark_reviewed, {
    req(rv$context_peak_flag, input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(rv$context_peak_flag$Position[1])
    snap <- curation_review_snapshot(r, pos)
    row <- data.frame(Action = "Mark reviewed / keep call", Position = pos,
                      Before = as.character(rv$context_peak_flag$Call[1]), After = as.character(rv$context_peak_flag$Call[1]),
                      Method = "Manual review", Evidence = flag_evidence_text(rv$context_peak_flag), Details = "No sequence change", stringsAsFactors = FALSE)
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), paste0("Reviewed position ", pos))
    commit_curated_result(input$inspect_sample, r2, paste0("Reviewed position ", pos))
    removeModal()
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_curation_action, {
    req(rv$pending_curation, input$inspect_sample)
    pnd <- rv$pending_curation
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(pnd$position)
    calls <- curated_raw_calls(r)
    before_call <- if (pos >= 1 && pos <= length(calls)) calls[pos] else ""
    if (identical(pnd$type, "base")) {
      new_base <- toupper(as.character(pnd$value))
      if (identical(new_base, before_call)) {
        snap <- curation_review_snapshot(r, pos)
        row <- data.frame(Action = "Mark reviewed / keep call", Position = pos, Before = before_call, After = before_call,
                          Method = "Manual review", Evidence = pnd$evidence, Details = "Selected base equals current call; no sequence change", stringsAsFactors = FALSE)
        label <- paste0("Reviewed position ", pos)
        r2 <- curation_commit(r, snap, row, selected_processing_settings(), label)
        commit_curated_result(input$inspect_sample, r2, label)
        rv$pending_curation <- NULL
        removeModal(); showNotification(paste0(label, "."), type = "message")
        return()
      }
      snap <- curation_set_base_snapshot(r, pos, new_base)
      if (is.null(snap)) { showNotification("Could not apply the requested base edit.", type = "error"); return() }
      method <- if (new_base %in% c("R","Y","S","W","K","M")) "Manual IUPAC ambiguity" else "Manual base correction"
      row <- data.frame(Action = "Base edit", Position = pos, Before = before_call, After = new_base,
                        Method = method, Evidence = pnd$evidence, Details = "Confirmed by user", stringsAsFactors = FALSE)
      label <- paste0("Base edit at ", pos, ": ", before_call, "->", new_base)
    } else {
      side <- as.character(pnd$value)
      st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
      snap <- curation_trim_snapshot(r, pos, side)
      if (is.null(snap)) { showNotification("The requested trimming boundary is not valid.", type = "error"); return() }
      row <- data.frame(Action = if (side == "left") "Manual trim left" else "Manual trim right", Position = pos,
                        Before = paste0(st, "-", en), After = paste0(snap$trim_start, "-", snap$trim_end),
                        Method = "Manual chromatogram curation", Evidence = pnd$evidence,
                        Details = "Flagged position excluded from retained sequence", stringsAsFactors = FALSE)
      label <- paste0(if (side == "left") "Trim left through " else "Trim right from ", pos)
    }
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), label)
    commit_curated_result(input$inspect_sample, r2, label)
    rv$pending_curation <- NULL
    removeModal()
    showNotification(paste0(label, "."), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$auto_correct_settings, {
    p <- ambiguous_peak_params_from_settings(selected_processing_settings())
    showModal(modalDialog(
      title = "Auto-correction criteria",
      div(class = "compact-help-note",
          span(icon("info-circle"), " These values control only automatic correction proposals. Flags and manual curation remain available even when a position does not meet these criteria.")),
      fluidRow(
        column(6,
          numericInput("auto_cfg_alt_called", "Alternative / current signal >=", value = p$auto_min_alt_to_called, min = 1.01, max = 10, step = 0.05),
          numericInput("auto_cfg_alt_third", "Alternative / third channel >=", value = p$auto_min_alt_to_third, min = 1.01, max = 10, step = 0.05)
        ),
        column(6,
          numericInput("auto_cfg_peak_offset", "Maximum peak offset (trace samples)", value = p$auto_max_peak_offset, min = 0, max = 10, step = 1),
          numericInput("auto_cfg_relative_signal", "Alternative signal / retained median >=", value = p$auto_min_relative_signal, min = 0.05, max = 3, step = 0.05)
        )
      ),
      footer = tagList(
        actionButton("reset_auto_correct_criteria", "Reset defaults", icon = icon("refresh")),
        modalButton("Cancel"),
        actionButton("apply_auto_correct_criteria", "Apply criteria", class = "btn-primary")
      ),
      size = "m", easyClose = TRUE
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$reset_auto_correct_criteria, {
    d <- ambiguous_peak_defaults()
    updateNumericInput(session, "auto_cfg_alt_called", value = d$auto_min_alt_to_called)
    updateNumericInput(session, "auto_cfg_alt_third", value = d$auto_min_alt_to_third)
    updateNumericInput(session, "auto_cfg_peak_offset", value = d$auto_max_peak_offset)
    updateNumericInput(session, "auto_cfg_relative_signal", value = d$auto_min_relative_signal)
  }, ignoreInit = TRUE)

  observeEvent(input$apply_auto_correct_criteria, {
    vals <- c(
      alt_called = suppressWarnings(as.numeric(input$auto_cfg_alt_called)),
      alt_third = suppressWarnings(as.numeric(input$auto_cfg_alt_third)),
      peak_offset = suppressWarnings(as.numeric(input$auto_cfg_peak_offset)),
      relative_signal = suppressWarnings(as.numeric(input$auto_cfg_relative_signal))
    )
    if (any(!is.finite(vals)) || vals[["alt_called"]] <= 1 || vals[["alt_third"]] <= 1 ||
        vals[["peak_offset"]] < 0 || vals[["peak_offset"]] > 10 ||
        vals[["relative_signal"]] <= 0 || vals[["relative_signal"]] > 3) {
      showNotification("Check the auto-correction criteria values.", type = "error")
      return()
    }
    settings <- rv$settings
    settings$auto_correct_min_alt_to_called <- vals[["alt_called"]]
    settings$auto_correct_min_alt_to_third <- vals[["alt_third"]]
    settings$auto_correct_max_peak_offset <- as.integer(round(vals[["peak_offset"]]))
    settings$auto_correct_min_relative_signal <- vals[["relative_signal"]]
    rv$settings <- settings
    for (nm in names(rv$results)) {
      if (!is.list(rv$results[[nm]]$processing_settings)) rv$results[[nm]]$processing_settings <- settings
      rv$results[[nm]]$processing_settings$auto_correct_min_alt_to_called <- vals[["alt_called"]]
      rv$results[[nm]]$processing_settings$auto_correct_min_alt_to_third <- vals[["alt_third"]]
      rv$results[[nm]]$processing_settings$auto_correct_max_peak_offset <- as.integer(round(vals[["peak_offset"]]))
      rv$results[[nm]]$processing_settings$auto_correct_min_relative_signal <- vals[["relative_signal"]]
    }
    rv$auto_correct_preview_df <- data.frame()
    removeModal()
    showNotification("Auto-correction criteria updated.", type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$auto_correct_preview, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    df <- high_confidence_autocorrections(r, selected_processing_settings())
    if (nrow(df) && length(r$curation$reviewed_positions)) df <- df[!df$Position %in% r$curation$reviewed_positions, , drop = FALSE]
    if (!nrow(df)) {
      showNotification("No positions meet the current auto-correction criteria. Use Criteria to adjust the thresholds.", type = "message")
      return()
    }
    df$Proposed_base <- df$Competing_channel
    rv$auto_correct_preview_df <- df
    showModal(modalDialog(
      title = "Preview high-confidence base corrections",
      p(paste0(nrow(df), " position(s) meet all conservative auto-correction criteria. Nothing changes until you confirm.")),
      tableOutput("auto_correct_preview_table"),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_auto_correct", paste0("Apply ", nrow(df), " correction(s)"), class = "btn-danger")),
      size = "l", easyClose = FALSE
    ))
  }, ignoreInit = TRUE)

  output$auto_correct_preview_table <- renderTable({
    df <- rv$auto_correct_preview_df
    if (!is.data.frame(df) || !nrow(df)) return(NULL)
    out <- df[, c("Position","Call","Proposed_base","Alternative_to_called_ratio","Alternative_to_third_ratio","Competitor_peak_offset","Competitor_signal"), drop = FALSE]
    names(out) <- c("Position","Current","Proposed","Alternative / current","Alternative / third","Peak offset","Signal")
    out
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  observeEvent(input$confirm_auto_correct, {
    req(input$inspect_sample)
    df <- rv$auto_correct_preview_df
    if (!is.data.frame(df) || !nrow(df)) return()
    r <- ensure_curation_state(selected_result())
    temp <- r
    rows <- list()
    for (i in seq_len(nrow(df))) {
      pos <- as.integer(df$Position[i]); new_base <- as.character(df$Proposed_base[i])
      calls <- curated_raw_calls(temp); before <- calls[pos]
      snap_i <- curation_set_base_snapshot(temp, pos, new_base)
      if (is.null(snap_i)) next
      temp <- curation_restore_snapshot(temp, snap_i, selected_processing_settings())
      ev <- paste0("alternative/current=", df$Alternative_to_called_ratio[i],
                   "; alternative/third=", df$Alternative_to_third_ratio[i],
                   "; peak offset=", df$Competitor_peak_offset[i],
                   "; alternative signal=", df$Competitor_signal[i])
      rows[[length(rows)+1L]] <- data.frame(Action="Base edit", Position=pos, Before=before, After=new_base,
                                            Method="Auto high-confidence correction", Evidence=ev,
                                            Details="Applied after preview and user confirmation", stringsAsFactors=FALSE)
    }
    if (!length(rows)) { removeModal(); showNotification("No valid corrections remained to apply.", type="warning"); return() }
    action_rows <- do.call(rbind, rows)
    final_snap <- curation_snapshot(temp)
    label <- paste0("Auto-corrected ", nrow(action_rows), " high-confidence position(s)")
    r2 <- curation_commit(r, final_snap, action_rows, selected_processing_settings(), label)
    commit_curated_result(input$inspect_sample, r2, label)
    rv$auto_correct_preview_df <- data.frame()
    removeModal()
    showNotification(label, type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$curation_undo, {
    req(input$inspect_sample)
    ans <- curation_undo(selected_result(), selected_processing_settings())
    if (!isTRUE(ans$changed)) { showNotification("Nothing to undo for this sample.", type = "message"); return() }
    commit_curated_result(input$inspect_sample, ans$result, paste0("Undo: ", ans$label))
    showNotification(paste0("Undid: ", ans$label), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$curation_redo, {
    req(input$inspect_sample)
    ans <- curation_redo(selected_result(), selected_processing_settings())
    if (!isTRUE(ans$changed)) { showNotification("Nothing to redo for this sample.", type = "message"); return() }
    commit_curated_result(input$inspect_sample, ans$result, paste0("Redo: ", ans$label))
    showNotification(paste0("Redid: ", ans$label), type = "message")
  }, ignoreInit = TRUE)

  output$curation_history_table <- renderDT({
    r <- ensure_curation_state(selected_result())
    lg <- r$curation$audit_log
    if (!is.data.frame(lg) || !nrow(lg)) return(datatable(data.frame(Message="No manual curation actions for this sample."), rownames=FALSE, options=list(dom="t")))
    datatable(lg, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE, autoWidth=TRUE, order=list(list(0,"desc"))))
  })

  observeEvent(input$curation_history, {
    req(input$inspect_sample)
    showModal(modalDialog(title = paste0("Curation history - ", input$inspect_sample), DTOutput("curation_history_table"), size="l", easyClose=TRUE, footer=modalButton("Close")))
  }, ignoreInit = TRUE)

  observeEvent(input$curation_reset, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    if (!length(r$curation$undo_stack) && !nrow(r$curation$base_edits) && identical(as.integer(r$curation$trim_start), as.integer(r$curation$auto_trim_start)) && identical(as.integer(r$curation$trim_end), as.integer(r$curation$auto_trim_end))) {
      showNotification("This sample is already at the automatic trimming/base-call state.", type="message"); return()
    }
    showModal(modalDialog(
      title="Reset sample to automatic processing?",
      p("This will restore the original automatic trim boundaries and remove all manual base edits/review marks from the active curated sequence. The reset itself is logged and can be undone."),
      footer=tagList(modalButton("Cancel"), actionButton("confirm_curation_reset", "Reset sample", class="btn-danger")), easyClose=FALSE
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_curation_reset, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    snap <- curation_reset_snapshot(r)
    row <- data.frame(Action="Reset to automatic processing", Position=NA_integer_, Before=paste0(r$summary$trim_start[1], "-", r$summary$trim_end[1]),
                      After=paste0(r$curation$auto_trim_start, "-", r$curation$auto_trim_end), Method="Manual reset", Evidence="",
                      Details="Manual trims, base edits and review marks cleared from active curation state", stringsAsFactors=FALSE)
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), "Reset to automatic processing")
    commit_curated_result(input$inspect_sample, r2, "Reset to automatic processing")
    removeModal(); showNotification("Sample reset to the automatic processing state.", type="message")
  }, ignoreInit = TRUE)

  output$chromatogram_plot <- plotly::renderPlotly({
    plot_result <- selected_result()
    plot_result$display_name <- qc_display_name(selected_sample_key())
    make_chromatogram_plotly(
      plot_result, selected_processing_settings(),
      flags = current_peak_flags(),
      show_flags = isTRUE(input$show_peak_flags)
    )
  })

  output$qc_plot <- renderPlot({
    draw_qc_metrics(selected_result(), selected_processing_settings())
  })

