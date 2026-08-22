  # ---------------- Taxonomic interpretation ----------------
  latest_blast_rid_for_sample <- function(original_name) {
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == original_name,,drop=FALSE]
    if (!nrow(jobs)) return("")
    jobs$rid[nrow(jobs)]
  }

  blast_hits_for_sample <- function(original_name) {
    rid <- latest_blast_rid_for_sample(original_name)
    if (!nzchar(rid) || !nrow(rv$blast_hits)) return(rv$blast_hits[0,,drop=FALSE])
    rv$blast_hits[rv$blast_hits$rid == rid,,drop=FALSE]
  }

  update_tax_sample_choices <- function(selected = NULL) {
    if (!nrow(rv$blast_hits)) {
      updateSelectInput(session, "tax_sample", choices = character())
      return(invisible(NULL))
    }
    pairs <- unique(rv$blast_hits[,c("original_name","final_name"),drop=FALSE])
    pairs <- pairs[nzchar(as.character(pairs$original_name)),,drop=FALSE]
    choices <- setNames(as.character(pairs$original_name), as.character(pairs$final_name))
    if (is.null(selected) || !selected %in% choices) selected <- if (length(choices)) choices[1] else character()
    updateSelectInput(session, "tax_sample", choices = choices, selected = selected)
  }

  observeEvent(input$pipeline_step, {
    if (identical(input$pipeline_step, "taxonomy") && nrow(rv$blast_hits)) {
      selected <- if (!is.null(input$tax_sample) && nzchar(input$tax_sample)) input$tax_sample else if (!is.null(input$blast_sample) && nzchar(input$blast_sample)) input$blast_sample else NULL
      update_tax_sample_choices(selected)
    }
  }, ignoreInit = TRUE)

  # Keep the taxonomy selector synchronized as each sequence result arrives,
  # including results retrieved through the batch BLAST workflow.
  observe({
    rv$blast_hits
    if (nrow(rv$blast_hits)) {
      selected <- isolate(input$tax_sample)
      update_tax_sample_choices(selected)
    }
  })


  observeEvent(input$to_taxonomy, {
    if (!nrow(rv$blast_hits)) {
      showNotification("Retrieve at least one BLAST result before taxonomic interpretation.", type="warning")
      return()
    }
    selected <- if (!is.null(input$blast_sample) && nzchar(input$blast_sample)) input$blast_sample else NULL
    update_tax_sample_choices(selected)
    updateTabsetPanel(session, "pipeline_step", selected="taxonomy")
  })

  analyze_taxonomy_sample <- function(sample_key, quiet = FALSE) {
    src <- normalize_blast_hits_unique_accession(blast_hits_for_sample(sample_key))
    if (!nrow(src)) return(list(ok=FALSE, message="No retrieved unique accession-level BLAST hits are available for this sequence."))

    rid <- unique(as.character(src$rid))[1]
    final_name <- unique(as.character(src$final_name))[1]
    top_n <- nrow(src)
    sample_settings <- settings_for_result(rv$results[[sample_key]])
    sample_target <- if (is.list(sample_settings)) stage2_scalar_text(sample_settings$target, "Other") else "Other"

    enriched <- NULL
    result <- NULL
    err <- NULL
    enriched <- tryCatch(
      enrich_hits_with_taxonomy(src),
      error=function(e) { err <<- conditionMessage(e); NULL }
    )
    if (!is.null(enriched)) {
      result <- tryCatch(
        build_taxonomic_consensus(enriched, target=sample_target, top_n=top_n),
        error=function(e) { err <<- conditionMessage(e); NULL }
      )
    }
    if (is.null(result) || !nrow(result$summary)) {
      return(list(ok=FALSE, message=if (!is.null(err)) paste("Taxonomic analysis failed:", err) else "Taxonomic analysis did not produce a summary."))
    }

    tax_error <- attr(enriched, "taxonomy_error")
    sm <- result$summary
    sm$final_name <- final_name
    sm$original_name <- sample_key
    sm$rid <- rid
    sm$target <- sample_target
    sm$max_hits_inspected <- top_n
    sm$analyzed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    sm$app_version <- APP_VERSION
    order_cols <- c("final_name","original_name","rid","target","recommended_identification","recommended_level","confidence")
    sm <- sm[,c(order_cols, setdiff(names(sm), order_cols)),drop=FALSE]

    th <- result$hits
    th$final_name <- final_name
    th$original_name <- sample_key
    th$rid <- rid

    tc <- result$counts
    tc$final_name <- final_name
    tc$original_name <- sample_key
    tc$rid <- rid

    if (nrow(rv$taxonomy_summary) && "rid" %in% names(rv$taxonomy_summary)) rv$taxonomy_summary <- rv$taxonomy_summary[rv$taxonomy_summary$rid != rid,,drop=FALSE]
    if (nrow(rv$taxonomy_hits) && "rid" %in% names(rv$taxonomy_hits)) rv$taxonomy_hits <- rv$taxonomy_hits[rv$taxonomy_hits$rid != rid,,drop=FALSE]
    if (nrow(rv$taxonomy_counts) && "rid" %in% names(rv$taxonomy_counts)) rv$taxonomy_counts <- rv$taxonomy_counts[rv$taxonomy_counts$rid != rid,,drop=FALSE]

    rv$taxonomy_summary <- rbind_fill(rv$taxonomy_summary, sm)
    rv$taxonomy_hits <- rbind_fill(rv$taxonomy_hits, th)
    rv$taxonomy_counts <- rbind_fill(rv$taxonomy_counts, tc)

    msg <- paste0(
      "Analysis complete for ", final_name, ": ",
      sm$recommended_identification[1], " · ", sm$recommended_level[1],
      " · confidence ", sm$confidence[1], "."
    )
    if (!is.null(tax_error) && nzchar(tax_error)) msg <- paste0(msg, " Some NCBI lineage lookups reported: ", tax_error)
    list(ok=TRUE, message=msg, final_name=final_name)
  }

  observeEvent(input$run_taxonomy, {
    req(input$tax_sample)
    rv$taxonomy_status_text <- paste0("Resolving NCBI taxonomy for ", input$tax_sample, "...")
    ans <- NULL
    withProgress(message="Taxonomic interpretation", value=0, {
      incProgress(0.2, detail="Resolving NCBI taxonomy and lineages")
      ans <- analyze_taxonomy_sample(input$tax_sample)
      incProgress(0.8, detail="Preparing consensus")
    })
    rv$taxonomy_status_text <- ans$message
    if (isTRUE(ans$ok)) showNotification("Taxonomic consensus completed.", type="message") else showNotification(ans$message, type="error")
  })

  observeEvent(input$run_taxonomy_all, {
    if (!nrow(rv$blast_hits)) {
      showNotification("Retrieve BLAST results before running batch taxonomy.", type="warning")
      return()
    }
    pairs <- unique(rv$blast_hits[,c("original_name","final_name"),drop=FALSE])
    samples <- as.character(pairs$original_name)
    samples <- samples[nzchar(samples)]
    if (!length(samples)) return()

    ok_n <- 0L; fail_n <- 0L; failures <- character()
    withProgress(message="Analyzing all retrieved sequences", value=0, {
      for (i in seq_along(samples)) {
        incProgress(1/length(samples), detail=paste0(i, "/", length(samples), " · ", pairs$final_name[match(samples[i], pairs$original_name)]))
        ans <- analyze_taxonomy_sample(samples[i], quiet=TRUE)
        if (isTRUE(ans$ok)) ok_n <- ok_n + 1L else { fail_n <- fail_n + 1L; failures <- c(failures, paste0(samples[i], ": ", ans$message)) }
        if (i < length(samples)) Sys.sleep(0.35)
      }
    })
    rv$taxonomy_batch_status_text <- paste0("Batch taxonomy complete: ", ok_n, " analyzed, ", fail_n, " failed.")
    rv$taxonomy_status_text <- rv$taxonomy_batch_status_text
    if (fail_n) showNotification(paste0(rv$taxonomy_batch_status_text, " Review failed samples individually."), type="warning") else showNotification(rv$taxonomy_batch_status_text, type="message")
    update_tax_sample_choices(isolate(input$tax_sample))
  })

  selected_tax_summary <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_summary)) return(rv$taxonomy_summary[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_summary[rv$taxonomy_summary$rid == rid,,drop=FALSE]
  })

  selected_tax_hits <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_hits)) return(rv$taxonomy_hits[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_hits[rv$taxonomy_hits$rid == rid,,drop=FALSE]
  })

  selected_tax_counts <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_counts)) return(rv$taxonomy_counts[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_counts[rv$taxonomy_counts$rid == rid,,drop=FALSE]
  })

  output$taxonomy_status <- renderUI({
    div(class="tax-note", strong("Status: "), rv$taxonomy_status_text)
  })

  output$taxonomy_hero <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)

    txt <- function(name, fallback = "—") {
      if (!name %in% names(sm)) return(fallback)
      x <- as.character(sm[[name]][1])
      if (is.na(x) || !nzchar(trimws(x))) fallback else x
    }
    num <- function(name) {
      if (!name %in% names(sm)) return(NA_real_)
      suppressWarnings(as.numeric(sm[[name]][1]))
    }
    fmt_pct <- function(x, digits=2) if (is.finite(x)) paste0(round(x, digits), "%") else "—"

    recommendation <- txt("recommended_identification", "Unresolved")
    level <- txt("recommended_level", "unresolved")
    confidence <- txt("confidence", "Low / review")
    conf_low <- tolower(confidence)
    conf_key <- if (grepl("high", conf_low)) "high" else if (grepl("moderate", conf_low)) "moderate" else "low"
    hero_class <- paste("taxonomy-result-hero", paste0("conf-", conf_key))
    hero_icon <- if (identical(conf_key, "high")) "check-circle" else if (identical(conf_key, "moderate")) "info-circle" else "exclamation-triangle"

    best_match <- txt("best_molecular_match", txt("species_candidate", ""))
    best_acc <- txt("best_match_accession", "")
    best_id <- num("best_match_identity_percent")
    if (!is.finite(best_id)) best_id <- num("candidate_identity_percent")
    best_cov <- num("best_match_query_coverage_percent")
    if (!is.finite(best_cov)) best_cov <- num("candidate_query_coverage_percent")

    alt <- txt("closest_alternative_species", txt("species_best_competitor", ""))
    alt_id <- num("closest_alternative_identity_percent")
    alt_cov <- num("closest_alternative_query_coverage_percent")
    locus_disc <- txt("locus_discrimination", "Not assessed")
    species_conclusion <- txt("species_level_conclusion", "")
    tax_support <- txt("taxonomic_support")
    seq_evidence <- txt("sequence_evidence")
    ref_support <- txt("reference_support")
    best_n <- num("best_species_database_accessions")
    alt_n <- num("closest_species_database_accessions")

    why_title <- if (tolower(level) == "genus" && nzchar(best_match) && best_match != "—") {
      "Why genus and not species?"
    } else if (tolower(level) == "species") {
      "Why this species?"
    } else if (tolower(level) == "unresolved") {
      "Why unresolved?"
    } else {
      "Why this recommendation?"
    }
    why_text <- txt("decision_reason", "Review the best molecular match and close alternatives below.")

    best_line <- if (nzchar(best_match) && best_match != "—") {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("flask")),
          div(strong("Best molecular match: "), tags$em(best_match),
              if (nzchar(best_acc) && best_acc != "—") paste0(" · ", best_acc) else "",
              paste0(" · Identity ", fmt_pct(best_id), " · Coverage ", fmt_pct(best_cov, 1))))
    } else NULL

    alt_line <- if (nzchar(alt) && alt != "—") {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("exchange-alt")),
          div(strong("Closest alternative: "), tags$em(alt),
              if (is.finite(alt_id) || is.finite(alt_cov)) paste0(" · Identity ", fmt_pct(alt_id), " · Coverage ", fmt_pct(alt_cov, 1)) else ""))
    } else {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("check")),
          div(strong("Closest alternative: "), "No different resolved species inside the close-match window"))
    }

    db_line <- if (is.finite(best_n) && best_n > 0) {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("database")),
          div(strong("Database representation: "), paste0(round(best_n), " accession(s) for the best species"),
              if (is.finite(alt_n) && alt_n > 0 && nzchar(alt) && alt != "—") paste0(" · ", round(alt_n), " for the closest alternative") else "",
              span(style="color:#64748b;", " · context only")))
    } else NULL

    div(class = hero_class,
      div(class = "taxonomy-hero-icon-wrap",
        div(class = "taxonomy-hero-icon", icon(hero_icon))
      ),
      div(class = "taxonomy-hero-main",
        div(class = "taxonomy-hero-title",
          recommendation,
          span(class = "level", paste0(" (", level, ")")),
          span(class = "taxonomy-hero-confidence", paste0("— ", confidence))
        ),
        div(class = "taxonomy-hero-lines",
          best_line,
          alt_line,
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("chart-line")),
              div(strong("Evidence: "), "Taxonomic support ", strong(tax_support), " · Sequence evidence ", strong(seq_evidence))),
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("microscope")),
              div(strong("Locus discrimination: "), locus_disc,
                  if (nzchar(species_conclusion) && species_conclusion != "—") paste0(" · ", species_conclusion) else "")),
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("book")),
              div(strong("Reference support: "), ref_support)),
          db_line
        )
      ),
      div(class = "taxonomy-why",
        div(class = "taxonomy-why-title", span(why_title), info_tip("The recommendation is based on the best molecular match and whether other named taxa are nearly indistinguishable by Identity and query coverage.")),
        div(class = "taxonomy-why-text", why_text)
      )
    )
  })

  output$taxonomy_summary_cards <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)

    sval <- function(name, fallback = "—") {
      if (!name %in% names(sm)) return(fallback)
      x <- as.character(sm[[name]][1])
      if (is.na(x) || !nzchar(trimws(x))) fallback else x
    }
    nval <- function(name) {
      if (!name %in% names(sm)) return(NA_real_)
      suppressWarnings(as.numeric(sm[[name]][1]))
    }
    pct <- function(x, digits=1) if (is.finite(x)) paste0(round(x, digits), "%") else "—"

    ident <- nval("best_match_identity_percent")
    if (!is.finite(ident)) ident <- nval("candidate_identity_percent")
    cov <- nval("best_match_query_coverage_percent")
    if (!is.finite(cov)) cov <- nval("candidate_query_coverage_percent")
    close_n <- nval("close_species_count")
    locus_disc <- sval("locus_discrimination", "Not assessed")
    seq_ev <- sval("sequence_evidence")
    overall <- sval("confidence")
    level <- sval("recommended_level", "unresolved")
    overall_sub <- paste0("Recommendation level: ", level)

    div(class = "taxonomy-metric-grid",
      div(class = "tax-metric metric-purple",
        div(class = "tax-metric-icon", icon("dna")),
        div(class = "tax-metric-label", "Best-match Identity"),
        div(class = "tax-metric-value metric-colored", pct(ident, 2)),
        div(class = "tax-metric-sub", "Best Identity within the near-best coverage band")
      ),
      div(class = "tax-metric metric-blue",
        div(class = "tax-metric-icon", icon("arrows-alt-h")),
        div(class = "tax-metric-label", "Query coverage"),
        div(class = "tax-metric-value metric-colored", pct(cov, 1)),
        div(class = "tax-metric-sub", "Coverage of the selected best molecular match")
      ),
      div(class = "tax-metric metric-green",
        div(class = "tax-metric-icon", icon("clone")),
        div(class = "tax-metric-label", "Close species alternatives"),
        div(class = "tax-metric-value", if (is.finite(close_n)) as.character(round(close_n)) else "—"),
        div(class = "tax-metric-sub", "Nearly indistinguishable named species")
      ),
      div(class = "tax-metric metric-cyan",
        div(class = "tax-metric-icon", icon("microscope")),
        div(class = "tax-metric-label", "Locus discrimination"),
        div(class = "tax-metric-value metric-colored", locus_disc),
        div(class = "tax-metric-sub", paste0("Sequence evidence: ", seq_ev))
      ),
      div(class = "tax-metric metric-amber",
        div(class = "tax-metric-icon", icon("shield")),
        div(class = "tax-metric-label", "Confidence in recommendation"),
        div(class = "tax-metric-value metric-colored", overall),
        div(class = "tax-metric-sub", overall_sub)
      )
    )
  })

  output$taxonomy_summary_table <- renderDT({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(datatable(data.frame(Message="Run multi-hit taxonomic analysis for this sequence."), rownames=FALSE, options=list(dom="t")))

    keep <- intersect(c(
      "algorithm_version","recommended_identification","recommended_level","confidence",
      "best_molecular_match","species_level_conclusion","locus_discrimination"
    ), names(sm))
    df <- sm[, keep, drop=FALSE]
    friendly <- c(
      algorithm_version="Decision engine",
      recommended_identification="Recommended identification",
      recommended_level="Recommended level",
      confidence="Confidence",
      best_molecular_match="Best molecular match",
      species_level_conclusion="Species-level conclusion",
      locus_discrimination="Locus discrimination"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(
      df,
      rownames=FALSE,
      class="compact stripe tax-decision-table",
      options=list(dom="t", autoWidth=FALSE, ordering=FALSE,
                   columnDefs=list(list(targets="_all", className="dt-wrap")))
    )
  })

  output$taxonomy_locus_note <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)
    decision <- if ("decision_reason" %in% names(sm) && nzchar(as.character(sm$decision_reason[1]))) as.character(sm$decision_reason[1]) else "Review the best molecular match and close alternatives."
    locus_flag <- if ("locus_flag" %in% names(sm)) as.character(sm$locus_flag[1]) else ""
    tagList(
      div(class="tax-callout",
        div(class="tax-callout-icon", icon("info-circle")),
        div(strong("Decision explanation"), decision)
      ),
      if (nzchar(locus_flag)) div(class="tax-callout",
        div(class="tax-callout-icon", icon("flag")),
        div(strong("Locus flag"), locus_flag)
      ) else NULL
    )
  })

  output$taxonomy_counts_table <- renderDT({
    df <- selected_tax_counts()
    if (!nrow(df)) return(datatable(data.frame(Message="No resolved species evidence profile yet."), rownames=FALSE, options=list(dom="t")))
    keep <- intersect(c(
      "interpretation","taxon","best_identity_percent","best_query_coverage_percent",
      "delta_identity_pp","delta_coverage_pp","accession_count"
    ), names(df))
    df <- df[,keep,drop=FALSE]
    friendly <- c(
      interpretation="Interpretation",
      taxon="Species",
      best_identity_percent="Best Identity (%)",
      best_query_coverage_percent="Best coverage (%)",
      delta_identity_pp="ΔIdentity (pp)",
      delta_coverage_pp="ΔCoverage (pp)",
      accession_count="Accessions"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(
      df,
      rownames=FALSE,
      class="compact stripe species-profile-table",
      options=list(dom="t", autoWidth=FALSE, paging=FALSE, ordering=FALSE,
                   columnDefs=list(
                     list(targets=0, width="115px"),
                     list(targets=1, width="155px", className="dt-taxon"),
                     list(targets=c(2,3,4,5,6), width="74px", className="dt-compact-number")
                   ))
    )
  })

  output$taxonomy_score_plot <- plotly::renderPlotly({
    df <- selected_tax_hits()
    if (!nrow(df)) return(plotly::plot_ly() |> plotly::layout(title="Run taxonomic analysis to display the BLAST score landscape."))

    df$analysis_rank <- if ("analysis_rank" %in% names(df)) suppressWarnings(as.numeric(df$analysis_rank)) else seq_len(nrow(df))
    df$bit_score_num <- suppressWarnings(as.numeric(df$bit_score))
    sp <- clean_taxon_text(df$species)
    gn <- clean_taxon_text(df$genus)
    df$display_taxon <- ifelse(nzchar(sp), sp, ifelse(nzchar(gn), gn, clean_taxon_text(df$organism)))
    resolved_genera <- unique(gn[nzchar(gn)])
    resolved_species <- unique(sp[nzchar(sp)])
    if (length(resolved_genera) > 1) {
      df$color_taxon <- ifelse(nzchar(gn), gn, "Unresolved genus")
      legend_title <- "Genus"
    } else if (length(resolved_species) > 1) {
      df$color_taxon <- ifelse(nzchar(sp), sp, "Unresolved species")
      legend_title <- "Species"
    } else {
      df$color_taxon <- ifelse(nzchar(df$display_taxon), df$display_taxon, "Unresolved")
      legend_title <- "Taxon"
    }

    refq <- if ("reference_quality" %in% names(df)) as.character(df$reference_quality) else ""
    acc <- if ("accession" %in% names(df)) as.character(df$accession) else ""
    ident <- if ("identity_percent" %in% names(df)) suppressWarnings(as.numeric(df$identity_percent)) else NA_real_
    cov <- if ("query_coverage_percent" %in% names(df)) suppressWarnings(as.numeric(df$query_coverage_percent)) else NA_real_
    status <- if ("molecular_evidence_status" %in% names(df)) as.character(df$molecular_evidence_status) else "Database match"
    is_best <- if ("is_best_molecular_match" %in% names(df)) as.logical(df$is_best_molecular_match) else rep(FALSE, nrow(df))
    df$marker_size <- ifelse(is_best, 12, 7)

    df$hover <- paste0(
      "<b>", df$display_taxon, "</b>",
      "<br>Accession: ", acc,
      "<br>Identity: ", round(ident,2), "%",
      "<br>Query coverage: ", round(cov,1), "%",
      "<br>Bit score: ", round(df$bit_score_num,2),
      "<br>Evidence: ", status,
      "<br>Reference: ", refq
    )

    p <- plotly::plot_ly(df, x=~analysis_rank, y=~bit_score_num, type="scatter", mode="lines",
                         line=list(width=1.5, color="#cbd5e1"), hoverinfo="skip", showlegend=FALSE)
    p <- plotly::add_markers(
      p, data=df, x=~analysis_rank, y=~bit_score_num,
      color=~color_taxon, size=~marker_size, sizes=c(7,12),
      text=~hover, hoverinfo="text",
      marker=list(line=list(width=0.5, color="#ffffff")),
      showlegend=TRUE
    )
    p <- plotly::layout(
      p,
      xaxis=list(title="BLAST hit rank", dtick=1, gridcolor="#edf1f6", zeroline=FALSE, titlefont=list(color="#52647a"), tickfont=list(color="#64748b")),
      yaxis=list(title="Bit score", gridcolor="#edf1f6", zeroline=FALSE, titlefont=list(color="#52647a"), tickfont=list(color="#64748b")),
      legend=list(title=list(text=legend_title), orientation="h", x=0, y=1.14, font=list(color="#475569")),
      hovermode="closest",
      paper_bgcolor="rgba(0,0,0,0)",
      plot_bgcolor="rgba(0,0,0,0)",
      margin=list(l=62, r=20, t=35, b=55)
    )
    plotly::config(
      p,
      displaylogo=FALSE,
      modeBarButtonsToRemove=c("lasso2d", "select2d", "hoverCompareCartesian")
    )
  })

  output$taxonomy_hits_table <- renderDT({
    df <- selected_tax_hits()
    if (!nrow(df)) return(datatable(data.frame(Message="No taxonomy-enriched hits yet."), rownames=FALSE, options=list(dom="t")))
    keep <- intersect(c(
      "analysis_rank","organism","genus","species","family","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","delta_bit_from_best",
      "reference_quality","hsp_count"
    ), names(df))
    df <- df[,keep,drop=FALSE]
    friendly <- c(
      analysis_rank="Rank",
      organism="NCBI organism",
      genus="Genus",
      species="Species",
      family="Family",
      accession="Accession",
      taxid="TaxID",
      identity_percent="Identity (%)",
      query_coverage_percent="Coverage (%)",
      evalue="E-value",
      bit_score="Bit score",
      delta_bit_from_best="ΔBit",
      reference_quality="Reference",
      hsp_count="HSP"
    )
    names(df) <- unname(friendly[names(df)])

    idx <- function(nm) { w <- which(names(df) == nm); if (length(w)) w[1]-1L else integer() }
    defs <- list(
      list(targets=idx("Rank"), width="42px", className="dt-compact-number"),
      list(targets=idx("NCBI organism"), width="105px", className="dt-organism-compact"),
      list(targets=idx("Genus"), width="85px", className="dt-taxon-compact"),
      list(targets=idx("Species"), width="125px", className="dt-taxon-compact"),
      list(targets=idx("Family"), width="90px", className="dt-wrap-compact"),
      list(targets=idx("Accession"), width="88px", className="dt-nowrap"),
      list(targets=idx("TaxID"), width="62px", className="dt-compact-number"),
      list(targets=unlist(lapply(c("Identity (%)","Coverage (%)","E-value","Bit score","ΔBit"), idx)), width="68px", className="dt-compact-number"),
      list(targets=idx("Reference"), width="92px", className="dt-wrap-compact"),
      list(targets=idx("HSP"), width="40px", className="dt-compact-number")
    )
    defs <- Filter(function(x) length(x$targets) > 0, defs)

    datatable(
      df, rownames=FALSE, filter="top", class="compact stripe taxonomy-hits-compact",
      options=list(
        pageLength=25, scrollX=FALSE, autoWidth=FALSE,
        columnDefs=defs
      )
    )
  })


  all_qc_peak_flags <- reactive({
    if (is.null(rv$results) || !length(rv$results)) return(data.frame())
    collect_ambiguous_peak_flags(rv$results, scope = "trimmed", settings = rv$settings)
  })

  all_curation_log <- reactive({
    if (is.null(rv$results) || !length(rv$results)) return(data.frame())
    collect_curation_log(rv$results)
  })

  team_summary_df <- reactive({
    base <- if (!is.null(rv$summary) && nrow(rv$summary)) export_summary_df() else data.frame()
    if (!nrow(base)) return(data.frame())

    out <- data.frame(
      Sample = as.character(base$final_name),
      Original_sample = as.character(base$sample_id),
      Target = if ("target" %in% names(base)) as.character(base$target) else if (!is.null(rv$settings)) rv$settings$target else "",
      Raw_length = if ("raw_length" %in% names(base)) base$raw_length else NA,
      Trimmed_length = if ("trimmed_length" %in% names(base)) base$trimmed_length else NA,
      QC_status = if ("status" %in% names(base)) as.character(base$status) else "",
      Ambiguous_peak_flags = 0L,
      Strong_peak_flags = 0L,
      Manual_curation_actions = 0L,
      Manual_base_edits = 0L,
      Curation_revision = 0L,
      Identification = "Not analyzed",
      Identification_level = "",
      Overall_confidence = "",
      Taxonomic_support = "",
      Sequence_evidence = "",
      Reference_support = "",
      Best_molecular_match = "",
      Best_match_identity_percent = NA_real_,
      Best_match_query_coverage_percent = NA_real_,
      Closest_alternative_species = "",
      Closest_alternative_identity_percent = NA_real_,
      Closest_alternative_query_coverage_percent = NA_real_,
      Close_species_count = NA_integer_,
      Locus_discrimination = "",
      Species_level_conclusion = "",
      Best_species_database_accessions = NA_integer_,
      Species_candidate = "",
      Species_candidate_confidence = "",
      Candidate_identity_percent = NA_real_,
      Candidate_query_coverage_percent = NA_real_,
      Genus_best_competitor = "",
      Genus_delta_bit = NA_real_,
      Species_best_competitor = "",
      Species_delta_bit = NA_real_,
      BLAST_hits_used = NA_integer_,
      Top_accession = "",
      Top_NCBI_organism = "",
      RID = "",
      Comment = "",
      Application_version = APP_VERSION,
      stringsAsFactors = FALSE
    )

    qcf <- all_qc_peak_flags()
    if (nrow(qcf)) {
      for (i in seq_len(nrow(out))) {
        consensus_record <- rv$consensus_set$records[[out$Original_sample[i]]]
        source_ids <- if (is.list(consensus_record)) as.character(consensus_record$source_read_ids) else out$Original_sample[i]
        sf <- qcf[qcf$Sample %in% source_ids, , drop = FALSE]
        if ("Review_status" %in% names(sf)) sf <- sf[sf$Review_status == "Active", , drop = FALSE]
        out$Ambiguous_peak_flags[i] <- nrow(sf)
        out$Strong_peak_flags[i] <- if (nrow(sf)) sum(sf$Severity == "Strong", na.rm = TRUE) else 0L
      }
    }
    clog <- all_curation_log()
    for (i in seq_len(nrow(out))) {
      nm <- out$Original_sample[i]
      consensus_record <- rv$consensus_set$records[[nm]]
      source_ids <- if (is.list(consensus_record)) as.character(consensus_record$source_read_ids) else nm
      source_results <- rv$results[intersect(source_ids, names(rv$results))]
      if (length(source_results)) {
        active_changes <- base_edits <- revisions <- integer()
        for (r in source_results) {
        r <- ensure_curation_state(r)
        base_edit_count <- if (is.data.frame(r$curation$base_edits)) nrow(r$curation$base_edits) else 0L
        trim_left_changed <- is.finite(as.numeric(r$curation$trim_start)) && is.finite(as.numeric(r$curation$auto_trim_start)) &&
          as.integer(r$curation$trim_start) != as.integer(r$curation$auto_trim_start)
        trim_right_changed <- is.finite(as.numeric(r$curation$trim_end)) && is.finite(as.numeric(r$curation$auto_trim_end)) &&
          as.integer(r$curation$trim_end) != as.integer(r$curation$auto_trim_end)
          active_changes <- c(active_changes, as.integer(base_edit_count + trim_left_changed + trim_right_changed))
          base_edits <- c(base_edits, base_edit_count)
          revisions <- c(revisions, as.integer(r$curation$revision))
        }
        out$Manual_curation_actions[i] <- sum(active_changes)
        out$Manual_base_edits[i] <- sum(base_edits)
        out$Curation_revision[i] <- if (length(revisions)) max(revisions) else 0L
      }
    }

    if (nrow(rv$taxonomy_summary)) {
      for (i in seq_len(nrow(out))) {
        sm <- rv$taxonomy_summary[rv$taxonomy_summary$original_name == out$Original_sample[i],,drop=FALSE]
        if (!nrow(sm)) next
        sm <- sm[nrow(sm),,drop=FALSE]
        out$Identification[i] <- as.character(sm$recommended_identification[1])
        out$Identification_level[i] <- as.character(sm$recommended_level[1])
        out$Overall_confidence[i] <- as.character(sm$confidence[1])
        if ("taxonomic_support" %in% names(sm)) out$Taxonomic_support[i] <- as.character(sm$taxonomic_support[1])
        if ("sequence_evidence" %in% names(sm)) out$Sequence_evidence[i] <- as.character(sm$sequence_evidence[1])
        if ("reference_support" %in% names(sm)) out$Reference_support[i] <- as.character(sm$reference_support[1])
        if ("best_molecular_match" %in% names(sm)) out$Best_molecular_match[i] <- as.character(sm$best_molecular_match[1])
        if ("best_match_identity_percent" %in% names(sm)) out$Best_match_identity_percent[i] <- sm$best_match_identity_percent[1]
        if ("best_match_query_coverage_percent" %in% names(sm)) out$Best_match_query_coverage_percent[i] <- sm$best_match_query_coverage_percent[1]
        if ("closest_alternative_species" %in% names(sm)) out$Closest_alternative_species[i] <- as.character(sm$closest_alternative_species[1])
        if ("closest_alternative_identity_percent" %in% names(sm)) out$Closest_alternative_identity_percent[i] <- sm$closest_alternative_identity_percent[1]
        if ("closest_alternative_query_coverage_percent" %in% names(sm)) out$Closest_alternative_query_coverage_percent[i] <- sm$closest_alternative_query_coverage_percent[1]
        if ("close_species_count" %in% names(sm)) out$Close_species_count[i] <- sm$close_species_count[1]
        if ("locus_discrimination" %in% names(sm)) out$Locus_discrimination[i] <- as.character(sm$locus_discrimination[1])
        if ("species_level_conclusion" %in% names(sm)) out$Species_level_conclusion[i] <- as.character(sm$species_level_conclusion[1])
        if ("best_species_database_accessions" %in% names(sm)) out$Best_species_database_accessions[i] <- sm$best_species_database_accessions[1]
        if ("species_candidate" %in% names(sm)) out$Species_candidate[i] <- as.character(sm$species_candidate[1])
        if ("species_candidate_confidence" %in% names(sm)) out$Species_candidate_confidence[i] <- as.character(sm$species_candidate_confidence[1])
        if ("candidate_identity_percent" %in% names(sm)) out$Candidate_identity_percent[i] <- sm$candidate_identity_percent[1]
        if ("candidate_query_coverage_percent" %in% names(sm)) out$Candidate_query_coverage_percent[i] <- sm$candidate_query_coverage_percent[1]
        if ("genus_best_competitor" %in% names(sm)) out$Genus_best_competitor[i] <- as.character(sm$genus_best_competitor[1])
        if ("genus_delta_bit" %in% names(sm)) out$Genus_delta_bit[i] <- sm$genus_delta_bit[1]
        if ("species_best_competitor" %in% names(sm)) out$Species_best_competitor[i] <- as.character(sm$species_best_competitor[1])
        if ("species_delta_bit" %in% names(sm)) out$Species_delta_bit[i] <- sm$species_delta_bit[1]
        if ("hits_used" %in% names(sm)) out$BLAST_hits_used[i] <- sm$hits_used[1]
        out$RID[i] <- as.character(sm$rid[1])
        if ("decision_reason" %in% names(sm)) out$Comment[i] <- as.character(sm$decision_reason[1])

        bh <- rv$blast_hits[rv$blast_hits$rid == out$RID[i],,drop=FALSE]
        if (nrow(bh)) {
          if ("rank" %in% names(bh)) bh <- bh[order(suppressWarnings(as.numeric(bh$rank))),,drop=FALSE]
          if ("accession" %in% names(bh)) out$Top_accession[i] <- as.character(bh$accession[1])
          if ("organism" %in% names(bh)) out$Top_NCBI_organism[i] <- as.character(bh$organism[1])
        }
      }
    }
    out
  })

  output$team_summary_table <- renderDT({
    df <- team_summary_df()
    if (!nrow(df)) return(datatable(data.frame(Message="No processed sequences yet."), rownames=FALSE, options=list(dom="t")))
    show <- intersect(c(
      "Sample","Target","Trimmed_length","QC_status","Ambiguous_peak_flags","Manual_curation_actions",
      "Identification","Identification_level","Overall_confidence","Best_molecular_match",
      "Best_match_identity_percent","Best_match_query_coverage_percent","Closest_alternative_species",
      "Close_species_count","Locus_discrimination","Reference_support","Comment"
    ), names(df))
    comment_idx <- which(show == "Comment") - 1L
    taxon_idx <- which(show %in% c("Identification","Best_molecular_match","Closest_alternative_species")) - 1L
    defs <- list()
    if (length(comment_idx)) defs[[length(defs)+1]] <- list(targets=comment_idx, width="420px", className="dt-comment")
    if (length(taxon_idx)) defs[[length(defs)+1]] <- list(targets=taxon_idx, width="145px", className="dt-taxon-compact")
    team_view <- df[,show,drop=FALSE]
    if ("Manual_curation_actions" %in% names(team_view)) names(team_view)[names(team_view) == "Manual_curation_actions"] <- "Active_curation_changes"
    datatable(team_view, rownames=FALSE, filter="top", class="compact stripe",
              options=list(pageLength=25, scrollX=TRUE, autoWidth=TRUE, columnDefs=defs))
  })

  output$download_team_summary_csv <- downloadHandler(
    filename=function() paste0(project_export_stem(), "_team_identification_summary.csv"),
    content=function(file) write.csv(team_summary_df(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_team_summary_xlsx <- downloadHandler(
    filename=function() paste0(project_export_stem(), "_team_identification_summary.xlsx"),
    content=function(file) {
      req(nrow(team_summary_df()))
      wb <- openxlsx::createWorkbook(creator="Sanger Sequence Pipeline")
      header_style <- openxlsx::createStyle(fgFill="#1F4E78", fontColour="#FFFFFF", textDecoration="bold", halign="center", valign="center")
      wrap_style <- openxlsx::createStyle(wrapText=TRUE, valign="top")

      add_sheet <- function(name, df, freeze=TRUE) {
        if (is.null(df) || !ncol(df)) df <- data.frame(Message="No data available", stringsAsFactors=FALSE)
        openxlsx::addWorksheet(wb, name)
        openxlsx::writeData(wb, name, df, withFilter=nrow(df)>0, headerStyle=header_style)
        if (freeze) openxlsx::freezePane(wb, name, firstRow=TRUE)
        if (nrow(df)) openxlsx::addStyle(wb, name, wrap_style, rows=2:(nrow(df)+1), cols=seq_len(ncol(df)), gridExpand=TRUE, stack=TRUE)
        openxlsx::setColWidths(wb, name, cols=seq_len(ncol(df)), widths="auto")
        long_cols <- which(names(df) %in% c("Comment","decision_reason","record_title","hit_title","NCBI hit title"))
        if (length(long_cols)) openxlsx::setColWidths(wb, name, cols=long_cols, widths=45)
      }

      add_sheet("Summary", team_summary_df())
      add_sheet("BLAST Hits", rv$blast_hits)
      add_sheet("Taxonomy Details", rv$taxonomy_summary)
      add_sheet("Species Evidence", rv$taxonomy_counts)
      add_sheet("Consensus Summary", rv$consensus_set$summary)
      add_sheet("Source Read QC", read_export_summary_df())
      add_sheet("QC Flags", all_qc_peak_flags())
      add_sheet("Manual Curation", all_curation_log())
      if (!is.null(rv$rename)) add_sheet("Rename Map", rv$rename)
      if (is.data.frame(rv$read_assignments) && nrow(rv$read_assignments)) add_sheet("Read Assignments", rv$read_assignments)
      if (is.list(rv$architecture)) {
        add_sheet("Assays", rv$architecture$assays)
        add_sheet("Isolates", rv$architecture$isolates)
        add_sheet("Loci", rv$architecture$loci)
        add_sheet("Reads", rv$architecture$reads)
      }
      settings_df <- data.frame(
        Field=c("Application version","Exported at", if (!is.null(rv$settings)) names(rv$settings) else character()),
        Value=c(APP_VERSION, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), if (!is.null(rv$settings)) unlist(rv$settings, use.names=FALSE) else character()),
        stringsAsFactors=FALSE
      )
      add_sheet("Run Settings", settings_df, freeze=FALSE)
      openxlsx::saveWorkbook(wb, file, overwrite=TRUE)
    }
  )
