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

  # NCBI BLAST is a shared service. Keep all automated BLAST contacts at least
  # 10 seconds apart, in addition to the per-RID one-minute polling rule.
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

    rv$blast_jobs <- rbind(rv$blast_jobs, data.frame(
      final_name=r$final_name,
      original_name=original_name,
      rid=parsed$rid,
      rtoe=parsed$rtoe,
      database=as.character(input$blast_database),
      hitlist_size=hitlist,
      consensus_revision=if (is.list(r$consensus$curation)) as.integer(r$consensus$curation$revision) else 0L,
      status="SUBMITTED",
      submitted_at=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      last_checked_at="",
      stringsAsFactors=FALSE
    ))

    list(ok=TRUE, rid=parsed$rid, rtoe=parsed$rtoe, final_name=r$final_name)
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
    showNotification(paste("Submitted to NCBI. RID:", ans$rid), type="message")
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
    rv$blast_batch_status_text <- paste0("Submitting ", length(sample_names), " processed sequence(s) to NCBI...")

    withProgress(message="Submitting all sequences to NCBI BLAST", value=0, {
      for (i in seq_along(sample_names)) {
        nm <- sample_names[i]
        requested_hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
        current_revision <- if (is.list(records[[nm]]$consensus$curation)) as.integer(records[[nm]]$consensus$curation$revision) else 0L
        idx <- matching_active_blast_job_index(nm, input$blast_database, requested_hitlist, current_revision)
        if (!is.na(idx)) {
          skipped <- skipped + 1L
          incProgress(
            1/length(sample_names),
            detail=paste("Skipping matching active job:", records[[nm]]$final_name,
                         paste0("(", input$blast_database, ", ", requested_hitlist, " hits)"))
          )
          next
        }

        incProgress(0, detail=paste("Submitting", records[[nm]]$final_name))
        ans <- submit_one_blast(nm)
        if (isTRUE(ans$ok)) {
          submitted <- submitted + 1L
        } else {
          failed <- failed + 1L
          failures <- c(failures, paste0(records[[nm]]$final_name, ": ", ans$message))
        }
        incProgress(1/length(sample_names))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Batch submission complete | submitted: ", submitted,
      " | skipped matching active jobs: ", skipped,
      " | failed: ", failed,
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

  retrieve_blast_job <- function(job_idx) {
    if (is.na(job_idx) || job_idx < 1 || job_idx > nrow(rv$blast_jobs)) {
      return(list(status="ERROR", contacted=FALSE, message="BLAST job index is invalid."))
    }
    if (identical(as.character(rv$blast_jobs$status[job_idx]), "STALE")) {
      return(list(status="STALE", contacted=FALSE, message="This RID belongs to a sequence version that was changed during manual curation. Submit the current curated sequence as a new BLAST job."))
    }

    original_name <- rv$blast_jobs$original_name[job_idx]
    records <- export_records()
    if (!original_name %in% names(records)) {
      return(list(status="ERROR", contacted=FALSE, message="Processed sequence is no longer available."))
    }
    r <- records[[original_name]]
    rid <- rv$blast_jobs$rid[job_idx]

    # Do not poll a single RID more often than once per minute.
    reference_time <- rv$blast_jobs$last_checked_at[job_idx]
    if (!nzchar(reference_time)) reference_time <- rv$blast_jobs$submitted_at[job_idx]
    ref <- suppressWarnings(as.POSIXct(reference_time, format="%Y-%m-%d %H:%M:%S"))
    if (!is.na(ref)) {
      elapsed <- as.numeric(difftime(Sys.time(), ref, units="secs"))
      if (is.finite(elapsed) && elapsed < 60) {
        return(list(
          status="TOO_SOON", contacted=FALSE,
          message=paste0("Wait ", ceiling(60-elapsed), " more seconds before checking RID ", rid, ".")
        ))
      }
    }

    wait_for_ncbi_contact_slot(10)
    rv$blast_jobs$last_checked_at[job_idx] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

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
      return(list(status="WAITING", contacted=TRUE, message="NCBI BLAST is still running."))
    }
    if (grepl("Status=FAILED|Status=UNKNOWN", txt)) {
      rv$blast_jobs$status[job_idx] <- "FAILED/UNKNOWN"
      return(list(status="FAILED", contacted=TRUE, message="NCBI reports failed or unknown job."))
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

    withProgress(message="Retrieving selected NCBI BLAST jobs", value=0, {
      for (k in seq_along(selected_rows)) {
        idx <- selected_rows[k]
        rid <- rv$blast_jobs$rid[idx]
        final_name <- rv$blast_jobs$final_name[idx]
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(as.character(rv$blast_jobs$status[idx]), "READY") && already_stored) {
          ans <- list(status="READY", message="Already retrieved.")
        } else {
          incProgress(0, detail=paste("Checking", final_name))
          ans <- retrieve_blast_job(idx)
        }

        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status == "WAITING") waiting <- waiting + 1L
        else if (ans$status == "TOO_SOON") too_soon <- too_soon + 1L
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(1/length(selected_rows))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Selected retrieval | jobs: ", length(selected_rows),
      " | ready: ", ready,
      " | still running: ", waiting,
      " | too soon to poll: ", too_soon,
      " | no parsed hits: ", no_hits,
      " | stale: ", stale,
      " | failed: ", failed
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
    rv$blast_batch_status_text <- paste0("Checking ", length(latest_indices), " BLAST job(s)...")

    withProgress(message="Retrieving submitted NCBI BLAST jobs", value=0, {
      for (k in seq_along(latest_indices)) {
        idx <- latest_indices[k]
        final_name <- rv$blast_jobs$final_name[idx]
        rid <- rv$blast_jobs$rid[idx]

        # A READY RID that is already present in the hit store does not need another server request.
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(rv$blast_jobs$status[idx], "READY") && already_stored) {
          ready <- ready + 1L
          incProgress(1/length(latest_indices), detail=paste(final_name, "already retrieved"))
          next
        }

        incProgress(0, detail=paste("Checking", final_name))
        ans <- retrieve_blast_job(idx)
        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status == "WAITING") waiting <- waiting + 1L
        else if (ans$status == "TOO_SOON") too_soon <- too_soon + 1L
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(1/length(latest_indices))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Batch retrieval | ready: ", ready,
      " | still running: ", waiting,
      " | too soon to poll: ", too_soon,
      " | no parsed hits: ", no_hits,
      " | stale: ", stale,
      " | failed: ", failed
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  output$blast_batch_status <- renderUI({
    div(class="tax-note", strong("Batch status: "), rv$blast_batch_status_text)
  })

  output$blast_job_status <- renderUI({
    req(input$blast_sample)
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
    if (!nrow(jobs)) return(p(class="settings-note", "No NCBI submission yet for this sequence."))
    j <- jobs[nrow(jobs), ]
    div(
      class=if(j$status == "READY") "status-ok" else "status-warning",
      paste("Selected sequence | RID:", j$rid, "| Status:", j$status,
            "| Database:", if (nzchar(as.character(j$database))) j$database else "legacy/unknown",
            "| Requested hits:", j$hitlist_size,
            "| Estimated wait (RTOE):", j$rtoe, "seconds")
    )
  })

  output$blast_jobs_table <- renderDT({
    df <- rv$blast_jobs
    if (!nrow(df)) return(datatable(data.frame(Message="No BLAST jobs submitted yet."), rownames=FALSE, options=list(dom="t")))
    keep <- c("final_name","original_name","rid","database","hitlist_size","consensus_revision","rtoe","status","submitted_at","last_checked_at")
    df <- df[, intersect(keep, names(df)), drop=FALSE]
    friendly <- c(
      final_name="Sample", original_name="Original sample", rid="RID", database="Database", hitlist_size="Hits requested",
      consensus_revision="Consensus revision", rtoe="Estimated wait (s)", status="Status", submitted_at="Submitted at", last_checked_at="Last checked at"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(df, rownames=FALSE, selection=list(mode="multiple", target="row"), options=list(pageLength=25, scrollX=TRUE, autoWidth=TRUE))
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

