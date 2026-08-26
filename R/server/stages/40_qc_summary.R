  # ---------------- QC summary ----------------
  qc_has_uploaded_ab1 <- function() {
    (!is.null(input$ab1_files) && is.data.frame(input$ab1_files) && nrow(input$ab1_files) > 0) ||
      (is.list(rv$results) && length(rv$results) > 0)
  }

  qc_workspace_empty_message <- function() {
    "Trim & QC has nothing to show, make sure there is at least one AB1 file uploaded."
  }

  output$qc_summary_cards <- renderUI({
    # Only warn when there is truly no AB1 evidence; do not spam after a normal upload.
    if (is.null(rv$summary) || !is.data.frame(rv$summary) || !nrow(rv$summary)) {
      if (!qc_has_uploaded_ab1()) {
        return(div(
          class = "status-warning",
          style = "padding:14px 16px; border:1px solid #f3d7a0; border-radius:10px; background:#fffaf0; margin-bottom:16px; line-height:1.45;",
          qc_workspace_empty_message()
        ))
      }
      return(NULL)
    }
    vals <- c(
      Total=nrow(rv$summary),
      OK=sum(rv$summary$status=="OK",na.rm=TRUE),
      Warnings=sum(rv$summary$status=="SHORT_AFTER_TRIMMING",na.rm=TRUE),
      Failed=sum(rv$summary$status %in% c("FAILED_TRIMMING","ERROR"),na.rm=TRUE)
    )
    fluidRow(lapply(names(vals), function(nm) column(3, div(class="summary-card", div(class="summary-number",vals[[nm]]), div(class="summary-label",nm)))))
  })

  output$summary_table <- renderDT({
    if (is.null(rv$summary) || !is.data.frame(rv$summary) || !nrow(rv$summary)) {
      if (!qc_has_uploaded_ab1()) {
        return(datatable(
          data.frame(Message = qc_workspace_empty_message(), stringsAsFactors = FALSE),
          rownames = FALSE, selection = "none", options = list(dom = "t")
        ))
      }
      req(FALSE)
    }
    df <- rv$summary[,c("sample_id","target","raw_length","trimmed_length","trim_start","trim_end","collapse_index","reason","median_peak_ratio_trimmed","status")]
    df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))
    names(df) <- c("Sample","Target","Raw length","Trimmed length","Start","End","Collapse","Reason","Median peak ratio","Status")
    datatable(df, rownames=FALSE, filter="top", options=list(pageLength=15,scrollX=TRUE,autoWidth=TRUE))
  })

  selected_sample_key <- reactive({
    # Quiet empty states: never validate() here (that paints red errors on plots).
    # The upload message is shown only on summary cards/tables.
    sid <- input$inspect_sample
    if (is.null(sid) || length(sid) != 1L) req(FALSE)
    sid <- as.character(sid)[1]
    if (!nzchar(sid) || !(sid %in% names(rv$results))) req(FALSE)
    sid
  })

  selected_result <- reactive({
    sid <- selected_sample_key()
    r <- rv$results[[sid]]
    # Stale selection must not invent another sample; require an exact key match.
    if (is.null(r)) req(FALSE)
    r
  })

  selected_processing_settings <- reactive({
    settings <- settings_for_result(selected_result())
    req(is.list(settings))
    settings
  })

  observeEvent(input$inspect_sample, {
    sid <- input$inspect_sample
    if (is.null(sid) || length(sid) != 1L) return()
    sid <- as.character(sid)[1]
    if (!nzchar(sid) || !(sid %in% names(rv$results))) return()
    r <- rv$results[[sid]]
    if (is.null(r)) return()
    updateTextAreaInput(session, "trimmed_sequence_preview", value = r$seq)
  }, ignoreInit = TRUE)

  output$sequence_metrics <- renderDT({
    datatable(make_sequence_preview(selected_result()), rownames = FALSE, options = list(dom = "t"))
  })
