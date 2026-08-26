  # ---------------- BLAST workspace ----------------
  resolve_blast_original_name <- function(value, records = export_records()) {
    if (is.null(value) || !length(value) || !length(records)) return(NULL)
    value <- as.character(value)[1]
    if (!nzchar(value)) return(NULL)
    if (value %in% names(records)) return(value)
    final_names <- vapply(records, function(x) as.character(x$final_name), character(1))
    idx <- which(final_names == value)
    if (length(idx)) return(names(records)[idx[1]])
    NULL
  }

  # Keep the BLAST sequence selector synchronized with the current processed/renamed
  # record set without resetting the user's selection every time the tab is opened.
  observe({
    if (!length(rv$results) || is.null(rv$rename)) {
      updateSelectInput(session, "blast_sample", choices = character(), selected = character())
      return()
    }
    rec <- export_records()
    if (!length(rec)) {
      updateSelectInput(session, "blast_sample", choices = character(), selected = character())
      return()
    }
    final_names <- vapply(rec, function(x) as.character(x$final_name), character(1))
    choices <- setNames(names(rec), final_names)
    current <- isolate(input$blast_sample)
    selected <- resolve_blast_original_name(current, rec)
    if (is.null(selected) && length(rec)) selected <- names(rec)[1]
    updateSelectInput(session, "blast_sample", choices=choices, selected=if (is.null(selected)) character() else selected)
  })

  observeEvent(input$to_blast, {
    gate_error <- stage3_consensus_gate_error(rv$consensus_set, rv$results)
    if (!is.null(gate_error)) {
      showNotification(gate_error, type = "error", duration = 10)
      return()
    }
    workflow_mark_unlocked("blast", also_export = TRUE)
    if (identical(rv$project_mode, "paired_consensus")) workflow_mark_unlocked("consensus")
    updateTabsetPanel(session, "pipeline_step", selected="blast")
  })

  blast_selected <- reactive({
    rec <- export_records()
    original_name <- resolve_blast_original_name(input$blast_sample, rec)
    req(!is.null(original_name))
    rec[[original_name]]
  })

  output$blast_sequence_preview_ui <- renderUI({
    r <- blast_selected()
    seq_text <- if (is.null(r$seq)) "" else as.character(r$seq)
    if (!nzchar(seq_text)) {
      return(tags$pre(class="sequence-code-viewer sequence-code-empty", "No processed sequence is available for this selection."))
    }
    tags$pre(class="sequence-code-viewer", wrap_sequence(seq_text, 80))
  })

  observeEvent(input$copy_blast_sequence,
               session$sendCustomMessage("copyText", list(text=blast_selected()$seq)))
  observeEvent(input$open_ncbi_blast,
               session$sendCustomMessage("openUrl", list(url="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn")))

  output$download_selected_blast <- downloadHandler(
    filename=function() paste0(clean_fasta_name(blast_selected()$final_name), ".fasta"),
    content=function(file) writeLines(
      paste0(">", clean_fasta_name(blast_selected()$final_name), "\n", wrap_sequence(blast_selected()$seq, 80)),
      file
    )
  )

  # NCBI BLAST is a shared service. Keep automated BLAST contacts at least
  # 10 seconds apart. Per-RID timing uses RTOE floor + auto next_poll_at
  # (see R/services/blast_polling.R), not a fixed one-minute rule.
  wait_for_ncbi_contact_slot <- function(min_seconds = 10) {
    last <- rv$ncbi_last_contact
    if (!is.null(last) && length(last) && !is.na(last)) {
      elapsed <- as.numeric(difftime(Sys.time(), last, units="secs"))
      if (is.finite(elapsed) && elapsed < min_seconds) {
        Sys.sleep(min_seconds - elapsed)
      }
    }
    rv$ncbi_last_contact <- Sys.time()
  }

  # Non-blocking slot check for the auto-poll observer (no Sys.sleep).
  ncbi_contact_slot_ready <- function(min_seconds = 10) {
    last <- rv$ncbi_last_contact
    if (is.null(last) || !length(last) || is.na(last)) return(TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), last, units = "secs"))
    is.finite(elapsed) && elapsed >= min_seconds
  }

  apply_blast_job_updates <- function(job_idx, updates) {
    if (!length(updates) || is.na(job_idx)) return(invisible(NULL))
    for (nm in names(updates)) {
      if (nm %in% names(rv$blast_jobs)) rv$blast_jobs[[nm]][job_idx] <- updates[[nm]]
    }
    invisible(NULL)
  }

  format_elapsed_seconds <- function(secs) {
    secs <- suppressWarnings(as.integer(round(secs)))
    if (!length(secs) || !is.finite(secs[1]) || secs[1] < 0) return("unknown")
    secs <- secs[1]
    if (secs < 60L) return(paste0(secs, "s"))
    mins <- secs %/% 60L
    rem <- secs %% 60L
    if (mins < 60L) return(sprintf("%dm %ds", mins, rem))
    hours <- mins %/% 60L
    mins <- mins %% 60L
    sprintf("%dh %dm", hours, mins)
  }

  parse_blast_timestamp <- function(value) {
    blast_parse_timestamp(value)
  }

  blast_elapsed_since <- function(timestamp_text) {
    t0 <- parse_blast_timestamp(timestamp_text)
    if (is.na(t0)) return(NA_real_)
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  blast_job_status_line <- function(job_row, extra = NULL) {
    job_row <- job_row[1, , drop = FALSE]
    status <- as.character(job_row$status[1])
    rid <- as.character(job_row$rid[1])
    elapsed <- blast_elapsed_since(job_row$submitted_at[1])
    rtoe <- as.character(job_row$rtoe[1])
    checked <- as.character(job_row$last_checked_at[1])
    db <- as.character(job_row$database[1])
    parts <- c(
      paste0("Status ", status),
      paste0("RID ", rid),
      paste0("elapsed ", format_elapsed_seconds(elapsed))
    )
    if (nzchar(rtoe) && !identical(rtoe, "NA")) parts <- c(parts, paste0("NCBI RTOE ~", rtoe, "s"))
    if (nzchar(db) && !identical(db, "NA")) parts <- c(parts, paste0("DB ", db))
    if (nzchar(checked)) parts <- c(parts, paste0("last checked ", checked))
    if (!is.null(extra) && nzchar(as.character(extra)[1])) parts <- c(parts, as.character(extra)[1])
    paste(parts, collapse = " | ")
  }

  parse_submit_response <- function(txt) {
    rid_match <- regmatches(txt, regexpr("RID = [A-Z0-9-]+", txt))
    rtoe_match <- regmatches(txt, regexpr("RTOE = [0-9]+", txt))
    rid <- if (length(rid_match) && nzchar(rid_match)) sub("RID = ", "", rid_match, fixed=TRUE) else ""
    rtoe <- if (length(rtoe_match) && nzchar(rtoe_match)) sub("RTOE = ", "", rtoe_match, fixed=TRUE) else NA_character_
    list(rid=rid, rtoe=rtoe)
  }

  latest_job_index_for_sample <- function(original_name) {
    idx <- which(rv$blast_jobs$original_name == original_name)
    if (!length(idx)) return(NA_integer_)
    tail(idx, 1)
  }

  matching_active_blast_job_index <- function(original_name, database, hitlist_size, consensus_revision) {
    if (!is.data.frame(rv$blast_jobs) || !nrow(rv$blast_jobs)) return(NA_integer_)
    idx <- which(
      rv$blast_jobs$original_name == original_name &
      rv$blast_jobs$status %in% c("SUBMITTED", "WAITING", "READY") &
      nzchar(as.character(rv$blast_jobs$database)) &
      as.character(rv$blast_jobs$database) == as.character(database) &
      suppressWarnings(as.integer(rv$blast_jobs$hitlist_size)) == as.integer(hitlist_size) &
      suppressWarnings(as.integer(rv$blast_jobs$consensus_revision)) == as.integer(consensus_revision)
    )
    if (!length(idx)) return(NA_integer_)
    tail(idx, 1)
  }

  submit_one_blast <- function(original_name) {
    records <- export_records()
    if (!original_name %in% names(records)) {
      return(list(ok=FALSE, message="Sequence is not available in the processed record set."))
    }
    r <- records[[original_name]]
    if (is.null(r) || !nzchar(r$seq)) {
      return(list(ok=FALSE, message="Processed sequence is empty."))
    }

    hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
    wait_for_ncbi_contact_slot(10)

    req_obj <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
      httr2::req_body_form(
        CMD="Put", PROGRAM="blastn", DATABASE=input$blast_database,
        QUERY=paste0(">", clean_fasta_name(r$final_name), "\n", r$seq),
        HITLIST_SIZE=hitlist
      ) |>
      httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")

    txt <- tryCatch(
      httr2::resp_body_string(httr2::req_perform(req_obj)),
      error=function(e) structure(list(error=conditionMessage(e)), class="blast_submit_error")
    )
    if (inherits(txt, "blast_submit_error")) {
      return(list(ok=FALSE, message=txt$error))
    }

    parsed <- parse_submit_response(txt)
    if (!nzchar(parsed$rid)) {
      return(list(ok=FALSE, message="NCBI did not return a BLAST RID."))
    }

    submitted_at <- Sys.time()
    submitted_txt <- blast_format_timestamp(submitted_at)
    rtoe_secs <- blast_parse_rtoe_seconds(parsed$rtoe)
    earliest_txt <- blast_schedule_next_poll_at(submitted_at, blast_rtoe_floor_seconds(rtoe_secs))
    next_auto_txt <- blast_schedule_next_poll_at(submitted_at, blast_first_auto_poll_delay_seconds(rtoe_secs))

    rv$blast_jobs <- rbind(rv$blast_jobs, data.frame(
      final_name=r$final_name,
      original_name=original_name,
      rid=parsed$rid,
      rtoe=parsed$rtoe,
      database=as.character(input$blast_database),
      hitlist_size=hitlist,
      consensus_revision=if (is.list(r$consensus$curation)) as.integer(r$consensus$curation$revision) else 0L,
      status="SUBMITTED",
      submitted_at=submitted_txt,
      last_checked_at="",
      auto_poll_enabled=TRUE,
      auto_poll_attempts=0L,
      next_poll_at=next_auto_txt,
      earliest_retrieve_at=earliest_txt,
      manual_retrieval_required=FALSE,
      stringsAsFactors=FALSE
    ))

    list(ok=TRUE, rid=parsed$rid, rtoe=parsed$rtoe, final_name=r$final_name,
         earliest_retrieve_at=earliest_txt, next_poll_at=next_auto_txt)
  }

  observeEvent(input$submit_ncbi_blast, {
    rec <- export_records()
    original_name <- resolve_blast_original_name(input$blast_sample, rec)
    req(!is.null(original_name))
    requested_hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
    current_revision <- if (is.list(rec[[original_name]]$consensus$curation)) as.integer(rec[[original_name]]$consensus$curation$revision) else 0L
    existing_idx <- matching_active_blast_job_index(original_name, input$blast_database, requested_hitlist, current_revision)
    if (!is.na(existing_idx)) {
      existing <- rv$blast_jobs[existing_idx, , drop=FALSE]
      showNotification(
        paste0("Matching BLAST job already exists (RID ", existing$rid[1], ", ",
               existing$database[1], ", ", existing$hitlist_size[1], " hits; status ",
               existing$status[1], "). No new request was submitted."),
        type="message", duration=8
      )
      return()
    }
    ans <- submit_one_blast(original_name)
    if (!isTRUE(ans$ok)) {
      showNotification(ans$message, type="error")
      return()
    }
    showNotification(
      paste0(
        "Submitted to NCBI. RID: ", ans$rid,
        " | first automatic check after RTOE (~", blast_parse_rtoe_seconds(ans$rtoe), "s)"
      ),
      type="message", duration=8
    )
  })

  observeEvent(input$submit_all_ncbi_blast, {
    records <- export_records()
    sample_names <- names(records)[vapply(records, function(x) !is.null(x$seq) && nzchar(x$seq), logical(1))]
    if (!length(sample_names)) {
      showNotification("No processed sequences are available for BLAST submission.", type="warning")
      return()
    }

    submitted <- 0L
    skipped <- 0L
    failed <- 0L
    failures <- character()
    batch_started <- Sys.time()
    n_jobs <- length(sample_names)
    rv$blast_batch_status_text <- paste0(
      "Submitting ", n_jobs, " sequence(s) to NCBI... | elapsed ", format_elapsed_seconds(0)
    )

    withProgress(message = paste0("Submitting ", n_jobs, " sequence(s) to NCBI BLAST"), value = 0, {
      for (i in seq_along(sample_names)) {
        nm <- sample_names[i]
        requested_hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
        current_revision <- if (is.list(records[[nm]]$consensus$curation)) as.integer(records[[nm]]$consensus$curation$revision) else 0L
        idx <- matching_active_blast_job_index(nm, input$blast_database, requested_hitlist, current_revision)
        if (!is.na(idx)) {
          skipped <- skipped + 1L
          incProgress(
            1/n_jobs,
            detail = paste0(i, "/", n_jobs, " | skip active job: ", records[[nm]]$final_name)
          )
          next
        }

        incProgress(0, detail = paste0(i, "/", n_jobs, " | submitting ", records[[nm]]$final_name))
        ans <- submit_one_blast(nm)
        if (isTRUE(ans$ok)) {
          submitted <- submitted + 1L
        } else {
          failed <- failed + 1L
          failures <- c(failures, paste0(records[[nm]]$final_name, ": ", ans$message))
        }
        incProgress(1/n_jobs, detail = paste0(i, "/", n_jobs, " | elapsed ", format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs"))))
      }
    })

    elapsed_txt <- format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs"))
    rv$blast_batch_status_text <- paste0(
      "Batch submission complete | submitted: ", submitted,
      " | skipped matching active jobs: ", skipped,
      " | failed: ", failed,
      " | elapsed ", elapsed_txt,
      if (length(failures)) paste0(" | ", paste(failures, collapse=" | ")) else ""
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  store_retrieved_hits <- function(hits, r, original_name, rid) {
    # v2.9 BLAST parser: one biological/database hit per accession. Multiple
    # HSPs (local alignment segments) must not receive independent weight in
    # taxonomy or sequence-evidence calculations.
    hits <- normalize_blast_hits_unique_accession(hits)
    hits$final_name <- r$final_name
    hits$original_name <- original_name
    hits$rid <- rid

    hit_cols <- c(
      "final_name","original_name","rid","rank","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","alignment_length","hsp_count","match_support"
    )
    for (nm in setdiff(hit_cols, names(hits))) hits[[nm]] <- NA
    hits <- hits[, hit_cols, drop=FALSE]

    if (nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits)) {
      rv$blast_hits <- rv$blast_hits[rv$blast_hits$rid != rid, , drop=FALSE]
    }
    rv$blast_hits <- rbind(rv$blast_hits, hits)

    top <- hits[1, , drop=FALSE]
    top_cols <- c(
      "final_name","original_name","rid","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","match_support"
    )

    # Preliminary identification is one current top-hit row per sequence.
    if (nrow(rv$blast_ids) && "original_name" %in% names(rv$blast_ids)) {
      rv$blast_ids <- rv$blast_ids[rv$blast_ids$original_name != original_name, , drop=FALSE]
    }
    rv$blast_ids <- rbind(rv$blast_ids, top[, top_cols, drop=FALSE])
  }

  retrieve_blast_job <- function(job_idx, mode = c("manual", "auto")) {
    mode <- match.arg(mode)
    if (is.na(job_idx) || job_idx < 1 || job_idx > nrow(rv$blast_jobs)) {
      return(list(status="ERROR", contacted=FALSE, message="BLAST job index is invalid."))
    }
    if (identical(as.character(rv$blast_jobs$status[job_idx]), "STALE")) {
      apply_blast_job_updates(job_idx, blast_disable_auto_poll_fields())
      return(list(status="STALE", contacted=FALSE, message="This RID belongs to a sequence version that was changed during manual curation. Submit the current curated sequence as a new BLAST job."))
    }

    original_name <- rv$blast_jobs$original_name[job_idx]
    records <- export_records()
    if (!original_name %in% names(records)) {
      return(list(status="ERROR", contacted=FALSE, message="Processed sequence is no longer available."))
    }
    r <- records[[original_name]]
    job_revision <- suppressWarnings(as.integer(rv$blast_jobs$consensus_revision[job_idx]))
    current_revision <- if (is.list(r$consensus$curation)) {
      suppressWarnings(as.integer(r$consensus$curation$revision[1]))
    } else {
      0L
    }
    if (!is.finite(job_revision)) job_revision <- 0L
    if (!is.finite(current_revision)) current_revision <- 0L
    if (!identical(as.integer(job_revision), as.integer(current_revision))) {
      rv$blast_jobs$status[job_idx] <- "STALE"
      apply_blast_job_updates(job_idx, blast_disable_auto_poll_fields())
      return(list(
        status = "STALE",
        contacted = FALSE,
        message = paste0(
          "This RID belongs to consensus revision ", job_revision,
          ", but the active analysis sequence is revision ", current_revision,
          ". Submit the current sequence as a new BLAST job."
        )
      ))
    }
    rid <- rv$blast_jobs$rid[job_idx]
    job <- rv$blast_jobs[job_idx, , drop = FALSE]
    now <- Sys.time()

    # RTOE floor: neither manual nor auto may contact NCBI before earliest_retrieve_at.
    if (!blast_manual_retrieve_allowed(job, now)) {
      remain <- blast_seconds_until(job$earliest_retrieve_at[1], now)
      if (!is.finite(remain)) {
        remain <- blast_rtoe_floor_seconds(job$rtoe[1]) - blast_elapsed_since(job$submitted_at[1])
      }
      return(list(
        status = "TOO_SOON",
        contacted = FALSE,
        message = blast_job_status_line(
          job,
          paste0(
            "wait ", ceiling(max(0, remain)), "s more (NCBI RTOE floor) before Check now / automatic poll"
          )
        )
      ))
    }

    # Auto path also respects next_poll_at (manual may check earlier after RTOE floor).
    if (identical(mode, "auto") && !blast_auto_poll_due(job, now)) {
      remain <- blast_seconds_until(job$next_poll_at[1], now)
      return(list(
        status = "TOO_SOON",
        contacted = FALSE,
        message = blast_job_status_line(
          job,
          paste0("next automatic check in ~", ceiling(max(0, remain)), "s")
        )
      ))
    }

    if (identical(mode, "auto")) {
      if (!ncbi_contact_slot_ready(10)) {
        return(list(status = "TOO_SOON", contacted = FALSE, message = "NCBI contact slot busy; will retry shortly."))
      }
      rv$ncbi_last_contact <- Sys.time()
    } else {
      wait_for_ncbi_contact_slot(10)
    }
    rv$blast_jobs$last_checked_at[job_idx] <- blast_format_timestamp(Sys.time())

    req_obj <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
      httr2::req_url_query(CMD="Get", RID=rid, FORMAT_TYPE="XML2") |>
      httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")

    txt <- tryCatch(
      httr2::resp_body_string(httr2::req_perform(req_obj)),
      error=function(e) structure(list(error=conditionMessage(e)), class="blast_retrieve_error")
    )
    if (inherits(txt, "blast_retrieve_error")) {
      return(list(status="ERROR", contacted=TRUE, message=txt$error))
    }

    if (grepl("Status=WAITING", txt)) {
      rv$blast_jobs$status[job_idx] <- "WAITING"
      job <- rv$blast_jobs[job_idx, , drop = FALSE]
      return(list(
        status = "WAITING",
        contacted = TRUE,
        message = blast_job_status_line(job, "NCBI is still computing this BLAST job")
      ))
    }
    if (grepl("Status=FAILED|Status=UNKNOWN", txt)) {
      rv$blast_jobs$status[job_idx] <- "FAILED/UNKNOWN"
      apply_blast_job_updates(job_idx, blast_disable_auto_poll_fields())
      job <- rv$blast_jobs[job_idx, , drop = FALSE]
      fail_kind <- if (grepl("Status=UNKNOWN", txt)) "UNKNOWN" else "FAILED"
      return(list(
        status = fail_kind,
        contacted = TRUE,
        message = blast_job_status_line(job, paste0("NCBI reports ", tolower(fail_kind), " job"))
      ))
    }

    rv$blast_jobs$status[job_idx] <- "READY"
    rv$blast_raw[[rid]] <- txt

    requested_hits <- suppressWarnings(as.integer(rv$blast_jobs$hitlist_size[job_idx]))
    if (!is.finite(requested_hits) || requested_hits < 1) requested_hits <- as.integer(input$blast_hitlist)
    requested_hits <- min(100L, max(1L, requested_hits))

    hits <- parse_blast_xml2_hits(
      txt,
      query_length=nchar(r$seq),
      max_hits=requested_hits
    )

    if (!nrow(hits)) {
      # Defensive fallback to NCBI's headerless 12-field tabular CSV.
      wait_for_ncbi_contact_slot(10)
      csv_req <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
        httr2::req_url_query(CMD="Get", RID=rid, FORMAT_TYPE="CSV", ALIGNMENT_VIEW="Tabular") |>
        httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")
      csv_txt <- tryCatch(
        httr2::resp_body_string(httr2::req_perform(csv_req)),
        error=function(e) NULL
      )
      if (!is.null(csv_txt)) {
        hits <- parse_blast_csv_hits_fallback(csv_txt, query_length=nchar(r$seq))
        if (nrow(hits)) {
          hits <- head(hits, requested_hits)
          rv$blast_raw[[rid]] <- csv_txt
        }
      }
    }

    hits <- normalize_blast_hits_unique_accession(hits)
    if (nrow(hits) > requested_hits) hits <- hits[seq_len(requested_hits), , drop = FALSE]

    if (!nrow(hits)) {
      return(list(status="NO_HITS", contacted=TRUE, message="BLAST is ready, but no unique accession-level hits could be parsed."))
    }

    # Enrich every unique accession whose XML/tabular result lacks display metadata.
    needs_meta <- which(
      (!nzchar(ifelse(is.na(hits$organism), "", hits$organism))) |
      (!nzchar(ifelse(is.na(hits$record_title), "", hits$record_title)))
    )

    if (length(needs_meta)) {
      accessions_needed <- unique(hits$accession[needs_meta])
      meta_df <- tryCatch(
        fetch_ncbi_nucleotide_metadata_batch(accessions_needed),
        error=function(e) data.frame()
      )

      if (nrow(meta_df)) {
        strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))
        for (i in needs_meta) {
          acc <- hits$accession[i]
          idx <- which(meta_df$accession == acc)
          if (!length(idx) && "primary_accession" %in% names(meta_df)) idx <- which(meta_df$primary_accession == strip_version(acc))
          if (!length(idx)) idx <- which(strip_version(meta_df$accession) == strip_version(acc))
          if (length(idx)) {
            m <- meta_df[idx[1], , drop=FALSE]
            if ((!nzchar(hits$organism[i]) || is.na(hits$organism[i])) && nzchar(m$organism[1])) hits$organism[i] <- m$organism[1]
            if ((!nzchar(hits$record_title[i]) || is.na(hits$record_title[i])) && nzchar(m$title[1])) hits$record_title[i] <- m$title[1]
            if ((is.na(hits$taxid[i]) || !is.finite(hits$taxid[i])) && !is.na(m$taxid[1])) hits$taxid[i] <- m$taxid[1]
          }
        }
      }
    }

    store_retrieved_hits(hits, r, original_name, rid)
    list(status="READY", contacted=TRUE, hits=nrow(hits), message=paste0(nrow(hits), " hit(s) retrieved."))
  }

  observeEvent(input$retrieve_ncbi_blast, {
    selected_rows <- input$blast_jobs_table_rows_selected
    selected_rows <- suppressWarnings(as.integer(selected_rows))
    selected_rows <- selected_rows[is.finite(selected_rows) & selected_rows >= 1L & selected_rows <= nrow(rv$blast_jobs)]

    # Backward-compatible fallback: if no job-table rows are selected, retrieve
    # the newest job for the sequence currently shown in Query workspace.
    if (!length(selected_rows)) {
      rec <- export_records()
      original_name <- resolve_blast_original_name(input$blast_sample, rec)
      if (is.null(original_name)) {
        showNotification("Select one or more BLAST jobs, or choose a sequence in Query workspace.", type="warning")
        return()
      }
      job_idx <- latest_job_index_for_sample(original_name)
      if (is.na(job_idx)) {
        showNotification("No BLAST job exists for the selected sequence.", type="warning")
        return()
      }
      selected_rows <- job_idx
    }

    # Preserve table order and avoid contacting the same RID twice if the UI ever
    # reports duplicate row selections.
    selected_rows <- unique(selected_rows)
    ready <- 0L; waiting <- 0L; too_soon <- 0L; failed <- 0L; no_hits <- 0L; stale <- 0L
    batch_started <- Sys.time()
    n_jobs <- length(selected_rows)
    detail_notes <- character()

    withProgress(message = paste0("Retrieving ", n_jobs, " selected NCBI BLAST job(s)"), value = 0, {
      for (k in seq_along(selected_rows)) {
        idx <- selected_rows[k]
        rid <- rv$blast_jobs$rid[idx]
        final_name <- rv$blast_jobs$final_name[idx]
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(as.character(rv$blast_jobs$status[idx]), "READY") && already_stored) {
          ans <- list(status="READY", message="Already retrieved.")
          apply_blast_job_updates(idx, blast_disable_auto_poll_fields())
        } else {
          incProgress(0, detail = paste0(k, "/", n_jobs, " | checking ", final_name))
          ans <- retrieve_blast_job(idx, mode = "manual")
          # Manual Check now never restarts auto-poll; terminal outcomes stop auto.
          if (ans$status %in% c("READY", "FAILED", "UNKNOWN", "NO_HITS", "STALE")) {
            apply_blast_job_updates(idx, blast_disable_auto_poll_fields())
          }
        }

        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status %in% c("WAITING")) {
          waiting <- waiting + 1L
          if (isTRUE(rv$blast_jobs$manual_retrieval_required[idx])) {
            detail_notes <- c(
              detail_notes,
              paste0(
                "NCBI status: SEARCHING | Result is not ready yet. | Last checked: ",
                rv$blast_jobs$last_checked_at[idx], " | RID: ", rid
              )
            )
          } else {
            detail_notes <- c(detail_notes, ans$message)
          }
        }
        else if (ans$status == "TOO_SOON") {
          too_soon <- too_soon + 1L
          detail_notes <- c(detail_notes, ans$message)
        }
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(
          1/n_jobs,
          detail = paste0(k, "/", n_jobs, " | elapsed ", format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs")))
        )
      }
    })

    elapsed_txt <- format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs"))
    rv$blast_batch_status_text <- paste0(
      "Selected retrieval | jobs: ", n_jobs,
      " | ready: ", ready,
      " | still running: ", waiting,
      " | too soon to poll: ", too_soon,
      " | no parsed hits: ", no_hits,
      " | stale: ", stale,
      " | failed: ", failed,
      " | elapsed ", elapsed_txt,
      if (length(detail_notes)) paste0(" | ", paste(unique(detail_notes), collapse = " || ")) else ""
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  observeEvent(input$retrieve_all_ncbi_blast, {
    if (!nrow(rv$blast_jobs)) {
      showNotification("No BLAST jobs have been submitted yet.", type="warning")
      return()
    }

    # Only the newest RID for each sequence is relevant to the current workspace.
    latest_indices <- unname(vapply(unique(rv$blast_jobs$original_name), latest_job_index_for_sample, integer(1)))
    latest_indices <- latest_indices[is.finite(latest_indices)]
    if (!length(latest_indices)) return()

    ready <- 0L; waiting <- 0L; too_soon <- 0L; failed <- 0L; no_hits <- 0L; stale <- 0L
    batch_started <- Sys.time()
    n_jobs <- length(latest_indices)
    detail_notes <- character()
    rv$blast_batch_status_text <- paste0(
      "Checking ", n_jobs, " BLAST job(s)... | elapsed ", format_elapsed_seconds(0)
    )

    withProgress(message = paste0("Retrieving ", n_jobs, " NCBI BLAST job(s)"), value = 0, {
      for (k in seq_along(latest_indices)) {
        idx <- latest_indices[k]
        final_name <- rv$blast_jobs$final_name[idx]
        rid <- rv$blast_jobs$rid[idx]

        # A READY RID that is already present in the hit store does not need another server request.
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(rv$blast_jobs$status[idx], "READY") && already_stored) {
          apply_blast_job_updates(idx, blast_disable_auto_poll_fields())
          ready <- ready + 1L
          incProgress(
            1/n_jobs,
            detail = paste0(k, "/", n_jobs, " | ", final_name, " already retrieved")
          )
          next
        }

        incProgress(0, detail = paste0(k, "/", n_jobs, " | checking ", final_name))
        ans <- retrieve_blast_job(idx, mode = "manual")
        if (ans$status %in% c("READY", "FAILED", "UNKNOWN", "NO_HITS", "STALE")) {
          apply_blast_job_updates(idx, blast_disable_auto_poll_fields())
        }
        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status == "WAITING") {
          waiting <- waiting + 1L
          detail_notes <- c(detail_notes, ans$message)
        }
        else if (ans$status == "TOO_SOON") too_soon <- too_soon + 1L
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(
          1/n_jobs,
          detail = paste0(k, "/", n_jobs, " | elapsed ", format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs")))
        )
      }
    })

    elapsed_txt <- format_elapsed_seconds(difftime(Sys.time(), batch_started, units = "secs"))
    rv$blast_batch_status_text <- paste0(
      "Batch retrieval | ready: ", ready,
      " | still running: ", waiting,
      " | too soon to poll: ", too_soon,
      " | no parsed hits: ", no_hits,
      " | stale: ", stale,
      " | failed: ", failed,
      " | elapsed ", elapsed_txt,
      if (length(detail_notes)) paste0(" | ", paste(unique(detail_notes), collapse = " || ")) else ""
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  output$blast_batch_status <- renderUI({
    div(class="tax-note", strong("Batch status: "), rv$blast_batch_status_text)
  })

  # PITAX auto-poll observer (NCBI): separate from status-card UI invalidateLater.
  # UI timers must not call retrieve_blast_job; only this observer contacts NCBI.
  # Two-phase tick: (1) publish "in progress" so the status card can paint,
  # (2) perform the NCBI Get on the next invalidateLater without a busy loop.
  observe({
    jobs <- rv$blast_jobs
    if (!is.data.frame(jobs) || !nrow(jobs)) {
      rv$blast_auto_activity_text <- ""
      rv$blast_auto_pending_idx <- NA_integer_
      return()
    }
    auto_pending <- isTRUE(any(
      jobs$auto_poll_enabled %in% TRUE &
        jobs$status %in% c("SUBMITTED", "WAITING")
    ))
    if (!auto_pending) {
      rv$blast_auto_activity_text <- ""
      rv$blast_auto_pending_idx <- NA_integer_
      return()
    }
    invalidateLater(2000, session)

    pending_idx <- suppressWarnings(as.integer(rv$blast_auto_pending_idx)[1])
    if (is.finite(pending_idx) && pending_idx >= 1L && pending_idx <= nrow(rv$blast_jobs)) {
      idx <- pending_idx
      rv$blast_auto_pending_idx <- NA_integer_
      withProgress(
        message = paste0("Automatic NCBI BLAST check: ", rv$blast_jobs$final_name[idx]),
        detail = paste0("RID ", rv$blast_jobs$rid[idx]),
        value = 0.4,
        {
          ans <- retrieve_blast_job(idx, mode = "auto")
          updates <- blast_job_after_auto_poll(rv$blast_jobs[idx, , drop = FALSE], ans$status, now = Sys.time())
          apply_blast_job_updates(idx, updates)
        }
      )
      rv$blast_auto_activity_text <- ""
      removeNotification(id = "blast_auto_check")

      if (identical(ans$status, "READY")) {
        rv$blast_batch_status_text <- paste0(
          "Automatic check complete | READY | ", rv$blast_jobs$final_name[idx],
          " | RID ", rv$blast_jobs$rid[idx]
        )
        showNotification(
          paste0("Automatic BLAST retrieval ready for ", rv$blast_jobs$final_name[idx],
                 " (RID ", rv$blast_jobs$rid[idx], ")."),
          type = "message", duration = 8
        )
      } else if (identical(ans$status, "WAITING") && isTRUE(rv$blast_jobs$manual_retrieval_required[idx])) {
        rv$blast_batch_status_text <- paste0(
          "Automatic checks stopped after 3 attempts | still SEARCHING at NCBI | RID ",
          rv$blast_jobs$rid[idx], " | use Check now"
        )
        showNotification(
          paste0(
            "NCBI has not returned the BLAST result after 3 automatic checks. ",
            "The request may still be running at NCBI. ",
            "Automatic checking has stopped. Use Check now to retrieve manually. ",
            "(RID ", rv$blast_jobs$rid[idx], ")"
          ),
          type = "warning", duration = 12
        )
      } else if (identical(ans$status, "WAITING")) {
        attempts_now <- suppressWarnings(as.integer(rv$blast_jobs$auto_poll_attempts[idx]))
        rv$blast_batch_status_text <- paste0(
          "Automatic check ", attempts_now, " of 3 | NCBI still SEARCHING | RID ",
          rv$blast_jobs$rid[idx]
        )
        showNotification(
          paste0(
            "Automatic NCBI check ", attempts_now, " of 3: still searching (RID ",
            rv$blast_jobs$rid[idx], ")."
          ),
          type = "message", duration = 6
        )
      } else if (ans$status %in% c("FAILED", "UNKNOWN")) {
        rv$blast_batch_status_text <- paste0("Automatic check | ", ans$status, " | RID ", rv$blast_jobs$rid[idx])
        showNotification(ans$message, type = "error", duration = 10)
      } else if (identical(ans$status, "ERROR")) {
        rv$blast_batch_status_text <- paste0(
          "Automatic check network error (RID kept) | RID ", rv$blast_jobs$rid[idx]
        )
        showNotification(
          paste0("Network error during automatic BLAST check (RID kept): ", ans$message),
          type = "warning", duration = 10
        )
      } else if (identical(ans$status, "TOO_SOON")) {
        rv$blast_batch_status_text <- paste0("Automatic check deferred | ", ans$message)
      }
      return()
    }

    if (!ncbi_contact_slot_ready(10)) return()
    now <- Sys.time()
    due_idx <- which(vapply(seq_len(nrow(jobs)), function(i) {
      blast_auto_poll_due(jobs[i, , drop = FALSE], now)
    }, logical(1)))
    if (!length(due_idx)) return()

    idx <- due_idx[[1]]
    rid_busy <- as.character(rv$blast_jobs$rid[idx])
    name_busy <- as.character(rv$blast_jobs$final_name[idx])
    attempt_next <- suppressWarnings(as.integer(rv$blast_jobs$auto_poll_attempts[idx]))
    if (!is.finite(attempt_next)) attempt_next <- 0L
    msg <- paste0(
      "Automatic NCBI check in progress: ", name_busy,
      " | RID ", rid_busy,
      " | check ", attempt_next + 1L, " of 3"
    )
    rv$blast_auto_activity_text <- msg
    rv$blast_batch_status_text <- msg
    # Sticky toast so the user sees activity even while the session is busy.
    showNotification(msg, id = "blast_auto_check", type = "message", duration = NULL)
    rv$blast_auto_pending_idx <- idx
    # Give the browser time to paint the notification/status before the blocking Get.
    invalidateLater(450, session)
  })

  output$blast_job_status <- renderUI({
    # UI refresh timer only - does not call retrieve_blast_job / NCBI.
    # do not invalidate the jobs DataTable; that desyncs scrollX header/body.
    if (is.data.frame(rv$blast_jobs) && nrow(rv$blast_jobs) &&
        any(rv$blast_jobs$status %in% c("SUBMITTED", "WAITING"))) {
      invalidateLater(2000, session)
    }
    req(input$blast_sample)
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
    if (!nrow(jobs)) return(p(class="settings-note", "No NCBI submission yet for this sequence."))
    j <- jobs[nrow(jobs), , drop = FALSE]
    banner <- blast_status_banner_parts(j, now = Sys.time(), format_elapsed = format_elapsed_seconds)
    activity <- as.character(rv$blast_auto_activity_text)[1]
    # Show in-progress activity for this RID, or any active auto-check message.
    show_activity <- !is.null(activity) && nzchar(activity)
    activity_line <- if (show_activity) {
      tags$div(
        class = "status-warning",
        style = "font-weight:700; margin-top:8px; padding:8px 10px; border:1px solid #f3d7a0; border-radius:8px; background:#fff7ed;",
        activity
      )
    } else {
      NULL
    }
    div(
      class = banner$class,
      strong("Selected sequence"),
      tags$br(),
      HTML(paste(banner$lines, collapse = "<br/>")),
      activity_line
    )
  })

  blast_jobs_display_df <- function() {
    df <- rv$blast_jobs
    if (!is.data.frame(df) || !nrow(df)) {
      return(data.frame(Message = "No BLAST jobs submitted yet.", stringsAsFactors = FALSE))
    }
    keep <- c(
      "final_name", "original_name", "rid", "database", "hitlist_size",
      "consensus_revision", "rtoe", "status", "submitted_at", "last_checked_at"
    )
    df <- df[, intersect(keep, names(df)), drop = FALSE]
    df$elapsed <- vapply(
      df$submitted_at,
      function(ts) format_elapsed_seconds(blast_elapsed_since(ts)),
      character(1)
    )
    # Stable column order - never rely on setdiff reordering.
    col_order <- c(
      "final_name", "original_name", "rid", "database", "hitlist_size",
      "consensus_revision", "rtoe", "status", "submitted_at", "last_checked_at", "elapsed"
    )
    df <- df[, col_order, drop = FALSE]
    names(df) <- c(
      "Sample", "Original sample", "RID", "Database", "Hits requested",
      "Consensus revision", "Estimated wait (s)", "Status",
      "Submitted at", "Last checked at", "Elapsed"
    )
    df
  }

  output$blast_jobs_table <- renderDT({
    df <- blast_jobs_display_df()
    if (identical(names(df), "Message")) {
      return(datatable(
        df, rownames = FALSE, selection = "none",
        options = list(dom = "t", scrollX = TRUE, autoWidth = FALSE)
      ))
    }
    # autoWidth=FALSE + CSS sync below: scrollX+autoWidth was shifting body
    # columns one slot right of the header. Live elapsed stays on the status
    # card (invalidateLater there); do not recreate this table on a timer.
    datatable(
      df,
      rownames = FALSE,
      selection = list(mode = "multiple", target = "row"),
      class = "display nowrap",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        autoWidth = FALSE,
        deferRender = TRUE,
        initComplete = htmlwidgets::JS(
          "function(settings, json) {",
          "  var api = this.api();",
          "  setTimeout(function() { api.columns.adjust(); }, 0);",
          "}"
        ),
        drawCallback = htmlwidgets::JS(
          "function(settings) {",
          "  var api = this.api();",
          "  api.columns.adjust();",
          "  var root = $(api.table().container());",
          "  var body = root.find('.dataTables_scrollBody');",
          "  var head = root.find('.dataTables_scrollHead');",
          "  if (body.length && head.length) head.scrollLeft(body.scrollLeft());",
          "}"
        )
      )
    )
  })

  output$blast_identification_table <- renderDT({
    # Deliberately show all retrieved sequences here: one current top hit per
    # sequence. The selected-sequence filter belongs to the detailed hit table.
    df <- rv$blast_ids
    if (!nrow(df)) {
      return(datatable(
        data.frame(Message="No parsed identification results yet. Submit and retrieve one or more BLAST jobs."),
        rownames=FALSE, options=list(dom="t")
      ))
    }

    df <- df[order(tolower(as.character(df$final_name))), , drop=FALSE]
    keep <- intersect(c(
      "final_name","organism","accession","identity_percent","query_coverage_percent",
      "evalue","bit_score","match_support","record_title","taxid","rid"
    ), names(df))
    df <- df[, keep, drop=FALSE]
    friendly <- c(
      final_name="Sample", organism="Organism / taxon", accession="Accession",
      identity_percent="Identity (%)", query_coverage_percent="Query coverage (%)",
      evalue="E-value", bit_score="Bit score", match_support="Match support",
      record_title="NCBI hit title", taxid="NCBI TaxID", rid="RID"
    )
    names(df) <- unname(friendly[names(df)])

    datatable(
      df, rownames=FALSE, filter="top",
      options=list(
        pageLength=25, scrollX=TRUE, autoWidth=TRUE,
        columnDefs=list(
          list(targets=1, className="dt-organism", width="190px"),
          list(targets=8, className="dt-hit-title", width="440px"),
          list(targets=c(0,2,3,4,5,6,7,9,10), className="dt-nowrap")
        )
      )
    )
  })

  output$blast_hits_table <- renderDT({
    df <- rv$blast_hits
    if (nrow(df) && !is.null(input$blast_sample)) {
      jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
      if (nrow(jobs)) {
        latest_rid <- jobs$rid[nrow(jobs)]
        df <- df[df$rid == latest_rid, , drop=FALSE]
      } else {
        df <- df[0, , drop=FALSE]
      }
    }
    if (!nrow(df)) return(datatable(data.frame(Message="No retrieved hits yet for the selected sequence."), rownames=FALSE, options=list(dom="t")))

    keep <- intersect(c(
      "rank","organism","accession","identity_percent","query_coverage_percent",
      "evalue","bit_score","hsp_count","match_support","record_title","taxid"
    ), names(df))
    df <- df[, keep, drop=FALSE]
    friendly <- c(
      rank="Rank", organism="Organism / taxon", accession="Accession",
      identity_percent="Identity (%)", query_coverage_percent="Query coverage (%)",
      evalue="E-value", bit_score="Bit score", hsp_count="HSP count", match_support="Match support",
      record_title="NCBI hit title", taxid="NCBI TaxID"
    )
    names(df) <- unname(friendly[names(df)])

    datatable(
      df, rownames=FALSE, filter="top",
      options=list(
        pageLength=25, scrollX=TRUE, autoWidth=TRUE,
        columnDefs=list(
          list(targets=1, className="dt-organism", width="190px"),
          list(targets=9, className="dt-hit-title", width="440px"),
          list(targets=c(0,2,3,4,5,6,7,8,10), className="dt-nowrap")
        )
      )
    )
  })

  output$blast_raw_preview <- renderText({
    req(input$blast_sample)
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
    if (!nrow(jobs)) return("")
    rid <- jobs$rid[nrow(jobs)]
    txt <- rv$blast_raw[[rid]]
    if (is.null(txt)) return("")
    substr(txt, 1, min(5000, nchar(txt)))
  })

  output$download_blast_jobs <- downloadHandler(
    filename=function() "NCBI_BLAST_jobs_and_identifications.csv",
    content=function(file) {
      jobs <- rv$blast_jobs
      if (nrow(rv$blast_ids)) jobs <- merge(jobs, rv$blast_ids, by=c("final_name","original_name","rid"), all.x=TRUE)
      write.csv(jobs, file, row.names=FALSE, fileEncoding="UTF-8")
    }
  )

  output$download_blast_hits <- downloadHandler(
    filename=function() "NCBI_BLAST_all_hits.csv",
    content=function(file) write.csv(rv$blast_hits, file, row.names=FALSE, fileEncoding="UTF-8")
  )

