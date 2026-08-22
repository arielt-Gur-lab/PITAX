  # ---------------- Stage 4 multi-locus profile ----------------
  observeEvent(input$to_multilocus, {
    updateTabsetPanel(session, "pipeline_step", selected = "multilocus")
  })

  observeEvent(input$back_taxonomy_from_multilocus, {
    updateTabsetPanel(session, "pipeline_step", selected = "taxonomy")
  })

  multilocus_status_label <- function(status) {
    labels <- c(
      CONCORDANT_SPECIES = "Concordant species",
      CONCORDANT_GENUS = "Concordant genus",
      GENUS_CONFLICT = "Genus conflict",
      SPECIES_CONFLICT = "Species conflict",
      PARTIAL_EVIDENCE = "Partial evidence",
      NO_TAXONOMY = "No taxonomy",
      SINGLE_LOCUS = "Single locus"
    )
    status <- stage4_scalar_text(status, "PARTIAL_EVIDENCE")
    if (status %in% names(labels)) unname(labels[[status]]) else status
  }

  multilocus_status_class <- function(status) {
    status <- stage4_scalar_text(status)
    if (identical(status, "CONCORDANT_SPECIES")) return("ml-concordant")
    if (identical(status, "CONCORDANT_GENUS")) return("ml-genus")
    if (status %in% c("GENUS_CONFLICT", "SPECIES_CONFLICT")) return("ml-conflict")
    "ml-review"
  }

  multilocus_metric_text <- function(x) {
    value <- suppressWarnings(as.numeric(x[1]))
    if (!length(value) || !is.finite(value)) return("-")
    paste0(format(round(value, 2), trim = TRUE, scientific = FALSE), "%")
  }

  multilocus_metric_width <- function(x) {
    value <- suppressWarnings(as.numeric(x[1]))
    if (!length(value) || !is.finite(value)) return("0%")
    paste0(max(0, min(100, value)), "%")
  }

  observeEvent(input$build_multilocus_profile, {
    projects <- list()
    source_names <- character()
    source_md5 <- character()

    if (isTRUE(input$multilocus_include_current)) {
      projects[[length(projects) + 1L]] <- make_project_bundle()
      source_names <- c(source_names, "Current session")
      source_md5 <- c(source_md5, "")
    }

    uploads <- input$multilocus_projects
    if (!is.null(uploads) && nrow(uploads)) {
      for (i in seq_len(nrow(uploads))) {
        loaded <- tryCatch(readRDS(uploads$datapath[i]), error = function(e) e)
        if (inherits(loaded, "error")) {
          showNotification(paste0("Could not read ", uploads$name[i], ": ", conditionMessage(loaded)), type = "error", duration = 12)
          return()
        }
        projects[[length(projects) + 1L]] <- loaded
        source_names <- c(source_names, as.character(uploads$name[i]))
        source_md5 <- c(source_md5, unname(tools::md5sum(uploads$datapath[i])))
      }
    }

    if (!length(projects)) {
      showNotification("Include the current project or add at least one saved PITAX project.", type = "warning")
      return()
    }

    built <- tryCatch(
      stage4_build_profile(projects, source_names = source_names, source_md5 = source_md5),
      error = function(e) e
    )
    if (inherits(built, "error")) {
      showNotification(paste0("Multi-locus profile was not built: ", conditionMessage(built)), type = "error", duration = 14)
      return()
    }
    rv$multilocus_profile <- built
    rv$project_status_text <- paste0("Unsaved Stage 4 profile: ", nrow(built$profiles), " isolate(s), ", nrow(built$evidence), " locus evidence row(s).")
    gate_error <- stage4_profile_gate_error(built)
    if (is.null(gate_error)) {
      showNotification("Multi-locus profile built. Review concordance and conflicts before reporting.", type = "message", duration = 8)
    } else {
      showNotification(gate_error, type = "warning", duration = 10)
    }
  })

  observe({
    profile <- stage4_ensure_profile(rv$multilocus_profile)
    isolates <- if (nrow(profile$profiles)) as.character(profile$profiles$Isolate) else character()
    isolates <- isolates[nzchar(trimws(isolates))]
    selected <- isolate(input$multilocus_isolate)
    if (!length(selected) || !selected %in% isolates) selected <- if (length(isolates)) isolates[1] else character()
    choices <- if (length(isolates)) stats::setNames(isolates, isolates) else character()
    updateSelectInput(session, "multilocus_isolate", choices = choices, selected = selected)
  })

  selected_multilocus_profile <- reactive({
    stage4_isolate_profile(rv$multilocus_profile, input$multilocus_isolate)
  })

  selected_multilocus_evidence <- reactive({
    stage4_isolate_evidence(rv$multilocus_profile, input$multilocus_isolate)
  })

  output$multilocus_overview_cards <- renderUI({
    overview <- stage4_profile_overview(rv$multilocus_profile)
    metric <- function(value, label) {
      div(class = "multilocus-overview-metric",
          div(class = "multilocus-overview-value", value),
          div(class = "multilocus-overview-label", label))
    }
    div(class = "multilocus-overview-metrics",
        metric(overview$Isolates[1], "Isolates"),
        metric(overview$Loci[1], "Distinct loci"),
        metric(overview$Concordant[1], "Concordant profiles"),
        metric(overview$Requires_Attention[1], "Require attention"))
  })

  output$multilocus_profile_hero <- renderUI({
    row <- selected_multilocus_profile()
    if (!nrow(row)) {
      return(div(class = "status-warning", "Build a profile, then choose an isolate to view its combined evidence."))
    }
    status <- stage4_scalar_text(row$Profile_Status)
    div(
      class = paste("multilocus-profile-hero", multilocus_status_class(status)),
      div(class = "multilocus-profile-title-row",
          div(class = "multilocus-profile-isolate", stage4_scalar_text(row$Isolate)),
          span(class = "multilocus-status-badge", multilocus_status_label(status))),
      div(class = "multilocus-profile-conclusion", stage4_scalar_text(row$Profile_Conclusion)),
      div(class = "multilocus-next-action",
          strong(paste0(stage4_scalar_text(row$Locus_Count, "0"), " loci | ", stage4_scalar_text(row$Taxonomy_Complete, "0/0"), " taxonomically interpreted")),
          tags$br(),
          strong("Next action: "), stage4_scalar_text(row$Next_Action))
    )
  })

  output$multilocus_locus_cards <- renderUI({
    evidence <- selected_multilocus_evidence()
    if (!nrow(evidence)) return(div(class = "status-warning", "No locus evidence is available for the selected isolate."))

    cards <- lapply(seq_len(nrow(evidence)), function(i) {
      row <- evidence[i, , drop = FALSE]
      analyzed <- identical(stage4_scalar_text(row$Taxonomy_Status), "Analyzed")
      call <- stage4_scalar_text(row$Recommended_Identification, "Not analyzed")
      rank <- stage4_scalar_text(row$Recommended_Level, "unresolved")
      confidence <- stage4_scalar_text(row$Confidence, "not assigned")
      best_match <- stage4_scalar_text(row$Best_Molecular_Match, "No molecular match recorded")
      accession <- stage4_scalar_text(row$Best_Match_Accession)
      reference <- stage4_scalar_text(row$Reference_Support, "Reference support not recorded")
      discrimination <- stage4_scalar_text(row$Locus_Discrimination, "Locus discrimination not recorded")
      rid <- stage4_scalar_text(row$RID, "No RID")

      bar <- function(label, value) {
        div(
          div(class = "multilocus-locus-bar-label", span(label), span(multilocus_metric_text(value))),
          div(class = "multilocus-locus-bar",
              div(class = "multilocus-locus-bar-fill", style = paste0("width:", multilocus_metric_width(value))))
        )
      }

      div(
        class = paste("multilocus-locus-card", if (analyzed) "locus-analyzed" else "locus-missing"),
        div(class = "multilocus-locus-title",
            span(stage4_scalar_text(row$Locus, "Unknown locus")),
            span(class = "multilocus-locus-state", if (analyzed) "Analyzed" else "Not analyzed")),
        div(class = "multilocus-locus-call", call),
        div(class = "multilocus-locus-meta", paste0("Rank: ", rank, " | Confidence: ", confidence)),
        div(class = "multilocus-locus-bars",
            bar("Identity", row$Best_Match_Identity),
            bar("Coverage", row$Best_Match_Coverage)),
        div(class = "multilocus-locus-meta", strong("Best match: "), best_match,
            if (nzchar(accession)) paste0(" | ", accession) else NULL),
        div(class = "multilocus-locus-meta", strong("Reference: "), reference),
        div(class = "multilocus-locus-meta", strong("Discrimination: "), discrimination),
        div(class = "multilocus-locus-meta", strong("Source: "), stage4_scalar_text(row$Source), " | ", rid)
      )
    })
    div(class = "multilocus-locus-grid", tagList(cards))
  })

  output$multilocus_evidence_plot <- plotly::renderPlotly({
    evidence <- selected_multilocus_evidence()
    if (!nrow(evidence)) {
      return(plotly::plot_ly() |>
        plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                       annotations = list(list(text = "Choose an isolate with locus evidence.", showarrow = FALSE))))
    }
    evidence$Identity <- suppressWarnings(as.numeric(evidence$Best_Match_Identity))
    evidence$Coverage <- suppressWarnings(as.numeric(evidence$Best_Match_Coverage))
    plotted <- evidence[is.finite(evidence$Identity) & is.finite(evidence$Coverage), , drop = FALSE]
    if (!nrow(plotted)) {
      return(plotly::plot_ly() |>
        plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                       annotations = list(list(text = "No completed BLAST and taxonomy metrics for this isolate.", showarrow = FALSE))))
    }
    plotted$Hover <- paste0(
      "<b>", plotted$Locus, "</b>",
      "<br>Call: ", ifelse(nzchar(plotted$Recommended_Identification), plotted$Recommended_Identification, "Unresolved"),
      "<br>Confidence: ", ifelse(nzchar(plotted$Confidence), plotted$Confidence, "Not assigned"),
      "<br>Identity: ", round(plotted$Identity, 2), "%",
      "<br>Coverage: ", round(plotted$Coverage, 2), "%",
      "<br>Accession: ", ifelse(nzchar(plotted$Best_Match_Accession), plotted$Best_Match_Accession, "Not recorded")
    )
    lower <- max(0, floor(min(c(plotted$Identity, plotted$Coverage), na.rm = TRUE) - 3))
    plotly::plot_ly(
      plotted, x = ~Coverage, y = ~Identity, color = ~Locus, text = ~Hover,
      type = "scatter", mode = "markers+text", textposition = "top center",
      marker = list(size = 12, line = list(color = "#ffffff", width = 1.5)),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        showlegend = FALSE,
        xaxis = list(title = "Query coverage (%)", range = c(lower, 100.5), gridcolor = "#e5e7eb"),
        yaxis = list(title = "Identity (%)", range = c(lower, 100.5), gridcolor = "#e5e7eb"),
        plot_bgcolor = "#fbfcfe", paper_bgcolor = "rgba(0,0,0,0)",
        margin = list(l = 58, r = 18, t = 28, b = 54)
      ) |>
      plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
  })

  multilocus_current_session_matches <- reactive({
    profile <- stage4_ensure_profile(rv$multilocus_profile)
    if (!nrow(profile$evidence) || !any(profile$evidence$Source == "Current session")) return(TRUE)
    stage4_current_project_matches(profile, make_project_bundle(), "Current session")
  })

  multilocus_export_ready <- reactive({
    is.null(stage4_profile_gate_error(rv$multilocus_profile)) && isTRUE(multilocus_current_session_matches())
  })

  output$multilocus_gate_status <- renderUI({
    profile <- stage4_ensure_profile(rv$multilocus_profile)
    error <- stage4_profile_gate_error(profile)
    if (!isTRUE(multilocus_current_session_matches())) {
      return(div(class = "status-error", "Warning: The current session changed after this profile was built. Rebuild Stage 4 before export."))
    }
    if (!is.null(error)) return(div(class = "status-warning", paste0("Warning: ", error)))
    conflict_n <- if (nrow(profile$profiles)) sum(profile$profiles$Profile_Status %in% c("GENUS_CONFLICT", "SPECIES_CONFLICT")) else 0L
    if (conflict_n > 0L) {
      return(div(class = "status-warning", paste0("OK: Structure is valid. ", conflict_n, " isolate profile(s) retain a biological conflict; no combined call is made for them.")))
    }
    div(class = "status-ok", "OK: At least two loci are present and every Isolate/Locus evidence row is unique.")
  })

  output$multilocus_sources_table <- renderDT({
    df <- stage4_ensure_profile(rv$multilocus_profile)$sources
    if (!nrow(df)) return(datatable(data.frame(Message = "Build a profile to populate source provenance."), rownames = FALSE, options = list(dom = "t")))
    view <- df
    view$Fingerprint <- ifelse(nzchar(view$Source_MD5), substr(view$Source_MD5, 1, 12), "current session")
    view$Source_MD5 <- NULL
    datatable(view, rownames = FALSE, class = "compact stripe", options = list(dom = "t", scrollX = TRUE, autoWidth = TRUE))
  })

  output$multilocus_profiles_table <- renderDT({
    df <- stage4_ensure_profile(rv$multilocus_profile)$profiles
    if (!nrow(df)) return(datatable(data.frame(Message = "No isolate profiles have been built."), rownames = FALSE, options = list(dom = "t")))
    datatable(
      df, rownames = FALSE, filter = "top", class = "compact stripe",
      options = list(pageLength = 20, scrollX = TRUE, autoWidth = TRUE,
                     columnDefs = list(list(targets = c(5, 10), width = "360px", className = "dt-wrap")))
    )
  })

  output$multilocus_evidence_table <- renderDT({
    df <- stage4_ensure_profile(rv$multilocus_profile)$evidence
    if (!nrow(df)) return(datatable(data.frame(Message = "No per-locus evidence has been imported."), rownames = FALSE, options = list(dom = "t")))
    show <- intersect(c(
      "Isolate", "Locus", "Final_Name", "Analysis_Status", "Sequence_Length", "Consensus_Revision",
      "Taxonomy_Status", "Recommended_Identification", "Recommended_Level", "Confidence",
      "Best_Molecular_Match", "Best_Match_Accession", "Best_Match_Identity", "Best_Match_Coverage",
      "Reference_Support", "Locus_Discrimination", "RID", "Source"
    ), names(df))
    datatable(df[, show, drop = FALSE], rownames = FALSE, filter = "top", class = "compact stripe",
              options = list(pageLength = 25, scrollX = TRUE, autoWidth = TRUE))
  })

  output$download_multilocus_profiles <- downloadHandler(
    filename = function() "PITAX_multi_locus_isolate_profiles.csv",
    content = function(file) {
      req(nrow(rv$multilocus_profile$profiles), isTRUE(multilocus_export_ready()))
      write.csv(rv$multilocus_profile$profiles, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_multilocus_evidence <- downloadHandler(
    filename = function() "PITAX_multi_locus_per_locus_evidence.csv",
    content = function(file) {
      req(nrow(rv$multilocus_profile$evidence), isTRUE(multilocus_export_ready()))
      write.csv(rv$multilocus_profile$evidence, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_multilocus_fasta <- downloadHandler(
    filename = function() "PITAX_multi_locus_sequences.fasta",
    content = function(file) {
      req(nrow(rv$multilocus_profile$evidence), isTRUE(multilocus_export_ready()))
      writeLines(stage4_make_fasta(rv$multilocus_profile), file)
    }
  )

  output$download_multilocus_checkpoint <- downloadHandler(
    filename = function() "PITAX_checkpoint_F_multi_locus.zip",
    content = function(file) {
      req(nrow(rv$multilocus_profile$evidence), isTRUE(multilocus_export_ready()))
      write_stage4_checkpoint_zip(file, rv$multilocus_profile)
    }
  )

  output$changelog_text <- renderText({
    path <- "CHANGELOG.md"
    if (!file.exists(path)) return("No changelog file found.")
    paste(readLines(path, warn=FALSE), collapse="\n")
  })

  output$download_taxonomy_summary <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_taxonomic_summary.csv"),
    content=function(file) write.csv(selected_tax_summary(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_taxonomy_hits <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_taxonomy_enriched_hits.csv"),
    content=function(file) write.csv(selected_tax_hits(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_taxonomy_checkpoint <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_checkpoint_E_taxonomy.zip"),
    content=function(file) {
      req(nrow(selected_tax_summary()))
      source_hits <- blast_hits_for_sample(input$tax_sample)
      make_taxonomy_checkpoint_zip(file, selected_tax_summary(), selected_tax_hits(), selected_tax_counts(), source_hits)
    }
  )

