  # ---------------- QC summary ----------------
  output$qc_summary_cards <- renderUI({
    req(rv$summary)
    vals <- c(
      Total=nrow(rv$summary),
      OK=sum(rv$summary$status=="OK",na.rm=TRUE),
      Warnings=sum(rv$summary$status=="SHORT_AFTER_TRIMMING",na.rm=TRUE),
      Failed=sum(rv$summary$status %in% c("FAILED_TRIMMING","ERROR"),na.rm=TRUE)
    )
    fluidRow(lapply(names(vals), function(nm) column(3, div(class="summary-card", div(class="summary-number",vals[[nm]]), div(class="summary-label",nm)))))
  })

  output$summary_table <- renderDT({
    req(rv$summary)
    df <- rv$summary[,c("sample_id","target","raw_length","trimmed_length","trim_start","trim_end","collapse_index","reason","median_peak_ratio_trimmed","status")]
    df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))
    names(df) <- c("Sample","Target","Raw length","Trimmed length","Start","End","Collapse","Reason","Median peak ratio","Status")
    datatable(df, rownames=FALSE, filter="top", options=list(pageLength=15,scrollX=TRUE,autoWidth=TRUE))
  })

  selected_sample_key <- reactive({
    sid <- input$inspect_sample
    req(!is.null(sid), length(sid) == 1L, nzchar(sid), sid %in% names(rv$results))
    sid
  })

  selected_result <- reactive({
    sid <- selected_sample_key()
    r <- rv$results[[sid]]
    req(!is.null(r))
    r
  })

  selected_processing_settings <- reactive({
    settings <- settings_for_result(selected_result())
    req(is.list(settings))
    settings
  })

  observeEvent(input$inspect_sample, {
    r <- selected_result()
    updateTextAreaInput(session, "trimmed_sequence_preview", value = r$seq)
  })

  output$sequence_metrics <- renderDT({
    datatable(make_sequence_preview(selected_result()), rownames = FALSE, options = list(dom = "t"))
  })

