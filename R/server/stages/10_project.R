  # ---------------- Project Save / Load ----------------
  output$project_status <- renderUI({
    tagList(
      span(rv$project_status_text),
      uiOutput("share_session_banner")
    )
  })

  make_project_bundle <- function() {
    list(
      format = "SangerSequencePipelineProject",
      schema_version = PROJECT_SCHEMA_VERSION,
      app_version = APP_VERSION,
      saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      active_tab = if (!is.null(input$pipeline_step)) input$pipeline_step else "upload",
      ui_state = list(
        inspect_sample = input$inspect_sample,
        consensus_sample = input$consensus_sample,
        blast_sample = input$blast_sample,
        tax_sample = input$tax_sample,
        blast_database = input$blast_database,
        blast_hitlist = input$blast_hitlist
      ),
      state = list(
        results = rv$results,
        summary = rv$summary,
        rename = rv$rename,
        settings = rv$settings,
        assay_profiles = rv$assay_profiles,
        project_defaults = rv$project_defaults,
        project_mode = rv$project_mode,
        read_assignments = rv$read_assignments,
        architecture = rv$architecture,
        consensus_set = rv$consensus_set,
        multilocus_profile = rv$multilocus_profile,
        migration_log = rv$project_migration_log,
        blast_jobs = rv$blast_jobs,
        blast_raw = rv$blast_raw,
        blast_ids = rv$blast_ids,
        blast_hits = normalize_blast_hits_unique_accession(rv$blast_hits),
        blast_batch_status_text = rv$blast_batch_status_text,
        taxonomy_summary = rv$taxonomy_summary,
        taxonomy_hits = rv$taxonomy_hits,
        taxonomy_counts = rv$taxonomy_counts,
        taxonomy_status_text = rv$taxonomy_status_text,
        taxonomy_batch_status_text = rv$taxonomy_batch_status_text
      )
    )
  }

  output$save_project <- downloadHandler(
    filename = function() {
      paste0(project_export_stem(), "_", format(Sys.Date(), "%Y%m%d"), ".sangerproject")
    },
    content = function(file) {
      saveRDS(make_project_bundle(), file = file, compress = "xz")
      rv$project_status_text <- paste0("Project saved from v", APP_VERSION, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ".")
    }
  )

  rebuild_blast_ids <- function() {
    if (!is.data.frame(rv$blast_hits) || !nrow(rv$blast_hits)) {
      rv$blast_ids <- data.frame()
      return(invisible(NULL))
    }
    hits <- normalize_blast_hits_unique_accession(rv$blast_hits)
    rv$blast_hits <- hits
    if (!all(c("original_name", "rid") %in% names(hits))) {
      rv$blast_ids <- data.frame()
      return(invisible(NULL))
    }
    parts <- split(seq_len(nrow(hits)), paste(hits$original_name, hits$rid, sep = "\r"))
    tops <- lapply(parts, function(idx) {
      g <- hits[idx, , drop = FALSE]
      if ("rank" %in% names(g)) g <- g[order(suppressWarnings(as.numeric(g$rank)), na.last = TRUE), , drop = FALSE]
      g[1, , drop = FALSE]
    })
    top <- do.call(rbind, tops)
    keep <- intersect(c(
      "final_name","original_name","rid","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","match_support"
    ), names(top))
    rv$blast_ids <- top[, keep, drop = FALSE]
    rownames(rv$blast_ids) <- NULL
  }

  invalidate_downstream_for_sample <- function(original_name, reason = "Curated sequence changed") {
    if (is.null(original_name) || !nzchar(original_name)) return(invisible(NULL))
    stale_rids <- character()
    had_blast <- FALSE
    had_taxonomy <- FALSE
    if (is.data.frame(rv$blast_jobs) && nrow(rv$blast_jobs) && "original_name" %in% names(rv$blast_jobs)) {
      idx <- which(rv$blast_jobs$original_name == original_name)
      if (length(idx)) {
        had_blast <- TRUE
        stale_rids <- unique(as.character(rv$blast_jobs$rid[idx]))
        rv$blast_jobs$status[idx] <- "STALE"
        if ("auto_poll_enabled" %in% names(rv$blast_jobs)) {
          rv$blast_jobs$auto_poll_enabled[idx] <- FALSE
          rv$blast_jobs$next_poll_at[idx] <- ""
          rv$blast_jobs$manual_retrieval_required[idx] <- FALSE
        }
      }
    }
    if (is.data.frame(rv$blast_hits) && nrow(rv$blast_hits)) {
      if ("original_name" %in% names(rv$blast_hits)) {
        hit_idx <- rv$blast_hits$original_name == original_name
        had_blast <- had_blast || any(hit_idx)
        stale_rids <- unique(c(stale_rids, as.character(rv$blast_hits$rid[hit_idx])))
        rv$blast_hits <- rv$blast_hits[rv$blast_hits$original_name != original_name, , drop = FALSE]
      } else if (length(stale_rids) && "rid" %in% names(rv$blast_hits)) {
        had_blast <- had_blast || any(rv$blast_hits$rid %in% stale_rids)
        rv$blast_hits <- rv$blast_hits[!rv$blast_hits$rid %in% stale_rids, , drop = FALSE]
      }
    }
    rebuild_blast_ids()
    filter_stale_tax <- function(df) {
      if (!is.data.frame(df) || !nrow(df)) return(df)
      if ("original_name" %in% names(df)) {
        idx <- df$original_name == original_name
        had_taxonomy <<- had_taxonomy || any(idx)
        return(df[!idx, , drop = FALSE])
      }
      if (length(stale_rids) && "rid" %in% names(df)) {
        idx <- df$rid %in% stale_rids
        had_taxonomy <<- had_taxonomy || any(idx)
        return(df[!idx, , drop = FALSE])
      }
      df
    }
    rv$taxonomy_summary <- filter_stale_tax(rv$taxonomy_summary)
    rv$taxonomy_hits <- filter_stale_tax(rv$taxonomy_hits)
    rv$taxonomy_counts <- filter_stale_tax(rv$taxonomy_counts)
    cause <- stage2_scalar_text(reason, "Curated sequence changed")
    if (had_blast) {
      rv$blast_batch_status_text <- paste0("BLAST evidence for ", original_name, " was marked stale because the active sequence changed (", cause, "). Re-run BLAST on the current sequence.")
    }
    if (had_taxonomy) {
      rv$taxonomy_status_text <- paste0("The previous taxonomic interpretation for ", original_name, " was removed because the active sequence changed (", cause, "). Re-run BLAST, then Taxonomy.")
    }
    invisible(list(blast = had_blast, taxonomy = had_taxonomy))
  }

  commit_curated_result <- function(sample_name, new_result, label = "Manual curation") {
    old_result <- rv$results[[sample_name]]
    old_seq <- if (!is.null(old_result$seq)) as.character(old_result$seq) else ""
    new_seq <- if (!is.null(new_result$seq)) as.character(new_result$seq) else ""
    rv$results[[sample_name]] <- new_result
    sync_summary_from_results()
    if (identical(input$inspect_sample, sample_name)) {
      updateTextAreaInput(session, "trimmed_sequence_preview", value = new_seq)
    }
    if (!identical(old_seq, new_seq)) {
      old_revision <- if (is.list(old_result$curation)) suppressWarnings(as.integer(old_result$curation$revision[1])) else NA_integer_
      new_revision <- if (is.list(new_result$curation)) suppressWarnings(as.integer(new_result$curation$revision[1])) else NA_integer_
      revision_note <- if (is.finite(old_revision) && is.finite(new_revision)) paste0(label, "; revision ", old_revision, " -> ", new_revision) else label
      # Clearing the whole consensus set must stale BLAST/taxonomy for every prior
      # analysis ID, not only the edited read's isolate. Partial invalidation left
      # READY jobs on other isolates while Stage 3 was empty.
      prior_consensus_ids <- if (is.list(rv$consensus_set) && length(rv$consensus_set$records)) {
        names(rv$consensus_set$records)
      } else {
        character()
      }
      rv$consensus_set <- stage3_empty_consensus_set()
      for (consensus_id in prior_consensus_ids) invalidate_downstream_for_sample(consensus_id, revision_note)
      invalidate_downstream_for_sample(sample_name, revision_note)
    }
    rv$project_status_text <- paste0("Unsaved curation change: ", sample_name, " | ", label, ".")
    invisible(!identical(old_seq, new_seq))
  }

  restore_settings_inputs <- function(settings) {
    if (is.null(settings) || !is.list(settings)) return(invisible(NULL))
    if (!is.null(settings$target)) updateSelectInput(session, "target", selected = settings$target)
    if (!is.null(settings$forward_primer)) updateTextInput(session, "forward_primer", value = settings$forward_primer)
    if (!is.null(settings$forward_primer_seq)) updateTextInput(session, "forward_primer_seq", value = settings$forward_primer_seq)
    if (!is.null(settings$reverse_primer)) updateTextInput(session, "reverse_primer", value = settings$reverse_primer)
    if (!is.null(settings$reverse_primer_seq)) updateTextInput(session, "reverse_primer_seq", value = settings$reverse_primer_seq)
    if (!is.null(settings$sequencing_primer)) updateRadioButtons(session, "sequencing_primer", selected = settings$sequencing_primer)
    if (!is.null(settings$enable_primer_mapping)) updateCheckboxInput(session, "enable_primer_mapping", value = isTRUE(settings$enable_primer_mapping))
    if (!is.null(settings$expected_amplicon_len)) updateNumericInput(session, "expected_amplicon_len", value = settings$expected_amplicon_len)
    if (!is.null(settings$absolute_max_base_index)) updateNumericInput(session, "absolute_max_base_index", value = settings$absolute_max_base_index)
    if (!is.null(settings$window)) updateNumericInput(session, "window", value = settings$window)
    if (!is.null(settings$min_peak_ratio)) updateNumericInput(session, "min_peak_ratio", value = settings$min_peak_ratio)
    if (!is.null(settings$min_relative_signal)) updateNumericInput(session, "min_relative_signal", value = settings$min_relative_signal)
    if (!is.null(settings$min_len_before_collapse)) updateNumericInput(session, "min_len_before_collapse", value = settings$min_len_before_collapse)
    if (!is.null(settings$bad_run_windows)) updateNumericInput(session, "bad_run_windows", value = settings$bad_run_windows)
    if (!is.null(settings$min_usable_len)) updateNumericInput(session, "min_usable_len", value = settings$min_usable_len)
  }

  # Restore a deserialized SangerSequencePipelineProject into this session.
  # Used by Load project and by Share Project (?share=TOKEN) so both stay identical.
  apply_loaded_project_object <- function(obj, source_label = "project", notify = TRUE) {
    if (!is.list(obj) || !identical(obj$format, "SangerSequencePipelineProject") || is.null(obj$state)) {
      msg <- "This file is not a valid Sanger Sequence Pipeline project."
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    source_schema <- if (is.null(obj$schema_version)) 1L else suppressWarnings(as.integer(obj$schema_version))
    if (length(source_schema) != 1L || is.na(source_schema) || source_schema < 1L) {
      msg <- "This project has an invalid schema version."
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    if (source_schema > PROJECT_SCHEMA_VERSION) {
      msg <- "This project was created by a newer project schema and cannot be loaded safely."
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    if (source_schema < 5L) {
      msg <- "This project predates schema 5 and cannot be migrated safely in Alpha 10. Open it with PITAX Alpha 9.1 and resave it first."
      if (notify) showNotification(msg, type = "error", duration = 12)
      return(list(ok = FALSE, message = msg))
    }

    st <- obj$state
    if (source_schema == 5L) {
      st <- tryCatch(
        assay_migrate_schema5_state(st),
        error = function(e) structure(list(error = conditionMessage(e)), class = "schema5_migration_error")
      )
      if (inherits(st, "schema5_migration_error")) {
        if (notify) showNotification(st$error, type = "error", duration = 14)
        return(list(ok = FALSE, message = st$error))
      }
    }
    loaded_profiles <- assay_coerce_profiles(st$assay_profiles)
    profile_error <- assay_validate_profiles(loaded_profiles)
    if (!is.null(profile_error)) {
      msg <- paste("Project assay profiles are invalid:", profile_error)
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    loaded_assignments <- stage2_coerce_assignments(st$read_assignments)
    assignment_error <- stage2_validate_assignments(loaded_assignments, assay_profiles = loaded_profiles)
    if (length(st$results) && !is.null(assignment_error)) {
      msg <- paste("Project read architecture is invalid:", assignment_error)
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    loaded_architecture <- if (nrow(loaded_assignments) && is.null(assignment_error)) tryCatch(
      stage2_build_architecture(
        loaded_assignments,
        project_id = if (is.list(st$architecture)) stage2_scalar_text(st$architecture$project_id, "project") else "project",
        assay_profiles = loaded_profiles
      ),
      error = function(e) NULL
    ) else NULL
    if (length(st$results) && nrow(loaded_assignments) && is.null(loaded_architecture)) {
      msg <- "Project read architecture could not be rebuilt safely."
      if (notify) showNotification(msg, type = "error", duration = 10)
      return(list(ok = FALSE, message = msg))
    }
    rv$results <- if (!is.null(st$results)) st$results else list()
    if (length(rv$results)) {
      for (nm in names(rv$results)) {
        rv$results[[nm]] <- ensure_curation_state(rv$results[[nm]])
        loaded_settings <- if (is.list(rv$results[[nm]]$processing_settings)) rv$results[[nm]]$processing_settings else st$settings
        rv$results[[nm]] <- curation_rebuild(rv$results[[nm]], loaded_settings)
      }
    }
    rv$summary <- st$summary
    rv$rename <- st$rename
    rv$settings <- st$settings
    rv$assay_profiles <- loaded_profiles
    rv$project_defaults <- assay_project_defaults_from_legacy_settings(st$project_defaults)
    rv$project_mode <- if (!is.null(st$project_mode) && st$project_mode %in% c("simple", "paired_consensus")) st$project_mode else "paired_consensus"
    rv$read_assignments <- loaded_assignments
    rv$assignment_signature <- assignment_state_signature(loaded_assignments)
    rv$architecture <- loaded_architecture
    rv$consensus_set <- stage3_ensure_consensus_set(if (is.list(st$consensus_set)) st$consensus_set else stage3_empty_consensus_set())
    rv$multilocus_profile <- stage4_ensure_profile(if (is.list(st$multilocus_profile)) st$multilocus_profile else stage4_empty_profile())
    rv$project_migration_log <- stage2_scalar_text(st$migration_log)
    rv$blast_jobs <- ensure_blast_jobs_schema(if (is.data.frame(st$blast_jobs)) st$blast_jobs else NULL)
    # Safe reload: do not resume aggressive automatic NCBI polling for pending RIDs.
    if (nrow(rv$blast_jobs)) {
      pending <- rv$blast_jobs$status %in% c("SUBMITTED", "WAITING")
      if (any(pending)) {
        rv$blast_jobs$auto_poll_enabled[pending] <- FALSE
        rv$blast_jobs$manual_retrieval_required[pending] <- TRUE
        rv$blast_jobs$next_poll_at[pending] <- ""
      }
    }
    rv$blast_raw <- if (is.list(st$blast_raw)) st$blast_raw else list()
    rv$blast_hits <- if (is.data.frame(st$blast_hits)) normalize_blast_hits_unique_accession(st$blast_hits) else data.frame()
    rv$blast_batch_status_text <- if (!is.null(st$blast_batch_status_text)) st$blast_batch_status_text else "Loaded project."
    rv$taxonomy_summary <- if (is.data.frame(st$taxonomy_summary)) st$taxonomy_summary else data.frame()
    rv$taxonomy_hits <- if (is.data.frame(st$taxonomy_hits)) st$taxonomy_hits else data.frame()
    rv$taxonomy_counts <- if (is.data.frame(st$taxonomy_counts)) st$taxonomy_counts else data.frame()
    rv$taxonomy_status_text <- if (!is.null(st$taxonomy_status_text)) st$taxonomy_status_text else "No taxonomic analysis has been run yet."
    rv$taxonomy_batch_status_text <- if (!is.null(st$taxonomy_batch_status_text)) st$taxonomy_batch_status_text else "No batch taxonomic analysis has been run yet."
    rv$ncbi_last_contact <- as.POSIXct(NA)
    sync_summary_from_results()
    rebuild_blast_ids()

    restore_settings_inputs(rv$settings)
    updateRadioButtons(session, "project_mode", selected = rv$project_mode)
    sync_project_mode_navigation(rv$project_mode)
    if (is.list(rv$consensus_set$settings)) {
      cs <- rv$consensus_set$settings
      if (!is.null(cs$min_overlap)) updateNumericInput(session, "consensus_min_overlap", value = cs$min_overlap)
      if (!is.null(cs$min_identity)) updateNumericInput(session, "consensus_min_identity", value = cs$min_identity)
      if (!is.null(cs$quality_delta)) updateNumericInput(session, "consensus_quality_delta", value = cs$quality_delta)
      if (!is.null(cs$strong_quality)) updateNumericInput(session, "consensus_strong_quality", value = cs$strong_quality)
    }

    result_names <- names(rv$results)
    saved_ui <- obj$ui_state
    inspect_selected <- if (!is.null(saved_ui$inspect_sample) && saved_ui$inspect_sample %in% result_names) saved_ui$inspect_sample else if (length(result_names)) result_names[1] else character()
    sync_qc_sample_choices(inspect_selected)
    consensus_ids <- names(rv$consensus_set$records)
    consensus_labels <- if (length(consensus_ids)) setNames(consensus_ids, vapply(rv$consensus_set$records, function(x) stage3_scalar_text(x$final_name), character(1))) else character()
    consensus_selected <- if (!is.null(saved_ui$consensus_sample) && saved_ui$consensus_sample %in% consensus_ids) saved_ui$consensus_sample else if (length(consensus_ids)) consensus_ids[1] else character()
    updateSelectInput(session, "consensus_sample", choices = consensus_labels, selected = consensus_selected)

    final_names <- character()
    if (length(result_names)) {
      result_labels <- result_names
      if (!is.null(rv$rename) && nrow(rv$rename)) {
        rename_idx <- match(result_names, rv$rename$Original_name)
        resolved <- !is.na(rename_idx) & nzchar(trimws(rv$rename$New_name[rename_idx]))
        result_labels[resolved] <- rv$rename$New_name[rename_idx[resolved]]
      }
      final_names <- setNames(result_names, result_labels)
    }
    blast_selected <- if (!is.null(saved_ui$blast_sample) && saved_ui$blast_sample %in% unname(final_names)) saved_ui$blast_sample else if (length(final_names)) unname(final_names)[1] else character()
    updateSelectInput(session, "blast_sample", choices = final_names, selected = blast_selected)
    if (!is.null(saved_ui$blast_database)) updateSelectInput(session, "blast_database", selected = saved_ui$blast_database)
    if (!is.null(saved_ui$blast_hitlist)) updateNumericInput(session, "blast_hitlist", value = saved_ui$blast_hitlist)

    if (nrow(rv$blast_hits)) {
      pairs <- unique(rv$blast_hits[, intersect(c("original_name","final_name"), names(rv$blast_hits)), drop = FALSE])
      if (all(c("original_name","final_name") %in% names(pairs))) {
        tax_choices <- setNames(as.character(pairs$original_name), as.character(pairs$final_name))
        tax_selected <- if (!is.null(saved_ui$tax_sample) && saved_ui$tax_sample %in% unname(tax_choices)) saved_ui$tax_sample else if (length(tax_choices)) unname(tax_choices)[1] else character()
        updateSelectInput(session, "tax_sample", choices = tax_choices, selected = tax_selected)
      }
    } else {
      updateSelectInput(session, "tax_sample", choices = character())
    }

    active <- if (!is.null(obj$active_tab) && obj$active_tab %in% c("upload","settings","qc","rename","consensus","export","blast","taxonomy","multilocus","help")) obj$active_tab else if (nrow(rv$multilocus_profile$profiles)) "multilocus" else if (nrow(rv$taxonomy_summary)) "taxonomy" else if (nrow(rv$blast_hits)) "blast" else if (length(rv$consensus_set$records)) "consensus" else if (length(rv$results)) "qc" else "upload"
    if (identical(rv$project_mode, "simple") && identical(active, "consensus")) active <- if (length(rv$consensus_set$records)) "blast" else "qc"
    unlocked <- c("upload")
    if (!is.null(input$ab1_files) || length(rv$results) || nrow(rv$read_assignments)) unlocked <- c(unlocked, "settings", "rename")
    if (length(rv$results)) unlocked <- c(unlocked, "qc")
    if (length(rv$consensus_set$records)) {
      unlocked <- c(unlocked, "blast", "export")
      if (identical(rv$project_mode, "paired_consensus")) unlocked <- c(unlocked, "consensus")
    }
    if (nrow(rv$blast_hits)) unlocked <- c(unlocked, "blast", "taxonomy", "export")
    if (nrow(rv$taxonomy_summary)) unlocked <- c(unlocked, "taxonomy", "multilocus")
    if (nrow(rv$multilocus_profile$profiles) || nrow(rv$multilocus_profile$evidence)) unlocked <- c(unlocked, "multilocus")
    unlocked <- c(unlocked, active)
    rv$workflow_unlocked <- unique(unlocked)
    completed <- character(0)
    if (nrow(loaded_assignments) || length(rv$results)) completed <- c(completed, "upload", "settings")
    if (length(rv$results)) completed <- c(completed, "rename", "qc")
    if (identical(rv$project_mode, "paired_consensus") && length(rv$consensus_set$records)) {
      completed <- c(completed, "consensus", "export")
    } else if (length(rv$consensus_set$records)) {
      completed <- c(completed, "export")
    }
    if (nrow(rv$blast_hits)) completed <- c(completed, "blast")
    if (nrow(rv$taxonomy_summary)) completed <- c(completed, "taxonomy")
    if (nrow(rv$multilocus_profile$profiles) || nrow(rv$multilocus_profile$evidence)) {
      completed <- c(completed, "multilocus")
    }
    rv$workflow_completed <- unique(completed)
    updateTabsetPanel(session, "pipeline_step", selected = active)

    rv$project_loaded_name <- as.character(source_label)[1]
    rv$project_status_text <- paste0(
      "Loaded ", source_label,
      " | saved with app v", ifelse(is.null(obj$app_version), "unknown", obj$app_version),
      if (!is.null(obj$saved_at)) paste0(" | saved ", obj$saved_at) else "",
      if (source_schema < PROJECT_SCHEMA_VERSION) paste0(" | migrated project schema ", source_schema, " -> ", PROJECT_SCHEMA_VERSION, " in memory") else "",
      "."
    )
    if (notify) {
      showNotification(
        if (source_schema < PROJECT_SCHEMA_VERSION) {
          paste0("Older project loaded and migrated to schema ", PROJECT_SCHEMA_VERSION, ". Save it to persist the migration.")
        } else {
          "Project loaded successfully."
        },
        type = "message", duration = 8
      )
    }
    list(ok = TRUE, message = rv$project_status_text, source_schema = source_schema)
  }

  observeEvent(input$load_project, {
    req(input$load_project$datapath)
    obj <- tryCatch(readRDS(input$load_project$datapath), error = function(e) structure(list(error = conditionMessage(e)), class = "project_load_error"))
    if (inherits(obj, "project_load_error")) {
      showNotification(paste("Could not load project:", obj$error), type = "error", duration = 10)
      return()
    }
    rv$share_banner_text <- ""
    apply_loaded_project_object(obj, source_label = input$load_project$name, notify = TRUE)
  })

  # ---------------- Share Project (immutable snapshot + tokenized URL) ----------------
  pitax_share_cleanup_expired()

  observeEvent(input$share_project, {
    if (!length(rv$results) && (!is.data.frame(rv$summary) || !nrow(rv$summary)) &&
        (!is.list(rv$consensus_set) || !length(rv$consensus_set$records))) {
      showNotification("Nothing to share yet. Process or load a project first.", type = "warning", duration = 8)
      return()
    }
    showModal(modalDialog(
      title = "Share Project",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("share_project_confirm", "Create share link", class = "btn-primary", icon = icon("link"))
      ),
      p("Anyone with this link can open an independent copy of this PITAX project. Changes are not synchronized."),
      radioButtons(
        "share_expiry_days",
        "Link expires after",
        choices = c("1 day" = "1", "3 days" = "3", "7 days" = "7"),
        selected = "3",
        inline = TRUE
      )
    ))
  })

  observeEvent(input$share_project_confirm, {
    removeModal()
    days <- pitax_share_parse_expiry_days(input$share_expiry_days)
    created <- tryCatch({
      pitax_share_create(
        make_project_bundle(),
        expiry_days = days,
        app_version = APP_VERSION,
        schema_version = PROJECT_SCHEMA_VERSION
      )
    }, error = function(e) structure(list(error = conditionMessage(e)), class = "share_create_error"))
    if (inherits(created, "share_create_error")) {
      showNotification(paste("Could not create share link:", created$error), type = "error", duration = 10)
      return()
    }
    url <- pitax_share_public_url(created$token, session = session)
    rv$last_share_url <- url
    showModal(modalDialog(
      title = "Share link created",
      easyClose = TRUE,
      footer = modalButton("Close"),
      p(strong("Expires: "), created$snapshot$expires_at),
      tags$pre(
        id = "share_project_url_text",
        style = "white-space:pre-wrap; word-break:break-all; background:#f8fafc; border:1px solid #e2e8f0; padding:10px; border-radius:8px;",
        url
      ),
      actionButton("share_project_copy", "Copy link", icon = icon("copy"), class = "btn-primary"),
      p(class = "compact-hint", style = "margin-top:12px;",
        "Anyone with this link can open an independent copy of this PITAX project. Changes are not synchronized.")
    ))
    rv$project_status_text <- paste0("Share link created (expires ", created$snapshot$expires_at, ").")
  })

  observeEvent(input$share_project_copy, {
    url <- as.character(rv$last_share_url)[1]
    if (!length(url) || is.na(url) || !nzchar(url)) return()
    session$sendCustomMessage("copyText", list(text = url))
    showNotification("Share link copied.", type = "message", duration = 4)
  })

  output$share_session_banner <- renderUI({
    txt <- as.character(rv$share_banner_text)[1]
    if (!length(txt) || is.na(txt) || !nzchar(txt)) return(NULL)
    div(
      class = "tax-note",
      style = "margin:10px 0 0; padding:10px 12px; border:1px solid #bfdbfe; border-radius:10px; background:#eff6ff;",
      HTML(gsub("\n", "<br/>", txt, fixed = TRUE))
    )
  })

  # Open ?share=TOKEN into a fresh independent session (immutable snapshot copy).
  # Must run inside a reactive consumer (observe), not session$onFlushed:
  # apply_loaded_project_object reads/writes rv$* and would abort outside a consumer.
  share_query_handled <- FALSE
  observe({
    search <- session$clientData$url_search
    if (isTRUE(share_query_handled)) return()
    # NULL means clientData not ready yet; "" means ready with no query.
    if (is.null(search)) return()
    share_query_handled <<- TRUE

    qs <- tryCatch(shiny::parseQueryString(search), error = function(e) list())
    token <- as.character(qs$share)[1]
    if (!length(token) || is.na(token) || !nzchar(token)) return()

    opened <- pitax_share_open(token, current_schema = PROJECT_SCHEMA_VERSION)
    if (!isTRUE(opened$ok)) {
      msg <- if (!is.null(opened$message)) opened$message else "Could not open the shared project."
      showNotification(msg, type = "error", duration = 14)
      rv$share_banner_text <- msg
      rv$project_status_text <- "Shared project link could not be loaded."
      session$sendCustomMessage("hideLoader", list())
      return()
    }
    ans <- apply_loaded_project_object(
      opened$project,
      source_label = paste0("shared snapshot ", substr(token, 1, 8), "..."),
      notify = FALSE
    )
    if (!isTRUE(ans$ok)) {
      rv$share_banner_text <- ans$message
      session$sendCustomMessage("hideLoader", list())
      return()
    }
    # Prefer the share banner; keep project status short.
    rv$project_status_text <- "Shared project session (independent copy)."
    rv$share_banner_text <- paste0(
      "Shared project loaded\n",
      "Created: ", opened$snapshot$created_at, "\n",
      "Changes in this session do not affect the original project."
    )
    showNotification(
      "Shared project loaded. Changes in this session do not affect the original project.",
      type = "message", duration = 10
    )
    session$sendCustomMessage("hideLoader", list())
  })
