  rv <- reactiveValues(
    results = list(), summary = NULL, rename = NULL, settings = NULL,
    assay_profiles = assay_default_profiles(),
    project_defaults = assay_project_defaults_from_legacy_settings(),
    read_assignments = stage2_empty_assignments(), architecture = NULL,
    consensus_set = stage3_empty_consensus_set(),
    multilocus_profile = stage4_empty_profile(),
    project_mode = "simple",
    assignment_signature = "",
    project_migration_log = "",
    blast_jobs = data.frame(
      final_name=character(), original_name=character(), rid=character(), rtoe=character(),
      database=character(), hitlist_size=integer(), consensus_revision=integer(), status=character(),
      submitted_at=character(), last_checked_at=character(),
      auto_poll_enabled=logical(), auto_poll_attempts=integer(),
      next_poll_at=character(), earliest_retrieve_at=character(),
      manual_retrieval_required=logical(),
      stringsAsFactors=FALSE
    ),
    ncbi_last_contact = as.POSIXct(NA),
    blast_batch_status_text = "No batch operation has been run yet.",
    blast_auto_activity_text = "",
    blast_auto_pending_idx = NA_integer_,
    share_banner_text = "",
    last_share_url = "",
    blast_raw = list(), blast_ids = data.frame(), blast_hits = data.frame(),
    taxonomy_summary = data.frame(), taxonomy_hits = data.frame(), taxonomy_counts = data.frame(),
    taxonomy_status_text = "No taxonomic analysis has been run yet.",
    taxonomy_batch_status_text = "No batch taxonomic analysis has been run yet.",
    project_status_text = "Current session has not been saved as a project.",
    project_loaded_name = "",
    workflow_unlocked = c("upload"),
    workflow_completed = character(0),
    context_peak_flag = NULL,
    pending_curation = NULL,
    auto_correct_preview_df = data.frame()
  )

  ensure_blast_jobs_schema <- function(df) {
    template <- data.frame(
      final_name=character(), original_name=character(), rid=character(), rtoe=character(),
      database=character(), hitlist_size=integer(), consensus_revision=integer(), status=character(),
      submitted_at=character(), last_checked_at=character(),
      auto_poll_enabled=logical(), auto_poll_attempts=integer(),
      next_poll_at=character(), earliest_retrieve_at=character(),
      manual_retrieval_required=logical(),
      stringsAsFactors=FALSE
    )
    if (!is.data.frame(df) || !nrow(df)) return(template)
    if (!"database" %in% names(df)) df$database <- ""
    if (!"hitlist_size" %in% names(df)) df$hitlist_size <- NA_integer_
    if (!"auto_poll_enabled" %in% names(df)) df$auto_poll_enabled <- FALSE
    if (!"auto_poll_attempts" %in% names(df)) df$auto_poll_attempts <- 0L
    if (!"next_poll_at" %in% names(df)) df$next_poll_at <- ""
    if (!"earliest_retrieve_at" %in% names(df)) df$earliest_retrieve_at <- ""
    if (!"manual_retrieval_required" %in% names(df)) df$manual_retrieval_required <- FALSE
    for (nm in setdiff(names(template), names(df))) df[[nm]] <- template[[nm]][NA_integer_]
    df <- df[, names(template), drop=FALSE]
    df$database <- as.character(df$database)
    df$hitlist_size <- suppressWarnings(as.integer(df$hitlist_size))
    df$consensus_revision <- suppressWarnings(as.integer(df$consensus_revision))
    df$auto_poll_enabled <- as.logical(df$auto_poll_enabled)
    df$auto_poll_enabled[is.na(df$auto_poll_enabled)] <- FALSE
    df$auto_poll_attempts <- suppressWarnings(as.integer(df$auto_poll_attempts))
    df$auto_poll_attempts[is.na(df$auto_poll_attempts)] <- 0L
    df$next_poll_at <- as.character(df$next_poll_at)
    df$next_poll_at[is.na(df$next_poll_at)] <- ""
    df$earliest_retrieve_at <- as.character(df$earliest_retrieve_at)
    df$earliest_retrieve_at[is.na(df$earliest_retrieve_at)] <- ""
    df$manual_retrieval_required <- as.logical(df$manual_retrieval_required)
    df$manual_retrieval_required[is.na(df$manual_retrieval_required)] <- FALSE
    # Backfill earliest_retrieve_at for legacy rows from submitted_at + RTOE floor.
    missing_earliest <- !nzchar(df$earliest_retrieve_at) & nzchar(as.character(df$submitted_at))
    if (any(missing_earliest)) {
      for (i in which(missing_earliest)) {
        df$earliest_retrieve_at[i] <- blast_schedule_next_poll_at(
          df$submitted_at[i],
          blast_rtoe_floor_seconds(df$rtoe[i])
        )
      }
    }
    df
  }

  settings_for_result <- function(result, fallback = rv$settings) {
    if (is.list(result) && is.list(result$processing_settings)) result$processing_settings else fallback
  }

  current_upload_source_ids <- function() {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return(NULL)
    vapply(input$ab1_files$name, stage2_read_stem, character(1))
  }

  initialize_current_read_assignments <- function() {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return(stage2_empty_assignments())
    default_profile <- assay_coerce_profiles(rv$assay_profiles)
    default_assay_id <- if (nrow(default_profile)) default_profile$Assay_ID[1] else ""
    stage2_make_read_assignments(
      input$ab1_files$name,
      default_locus = input$target,
      default_direction = input$sequencing_primer,
      forward_primer = input$forward_primer,
      reverse_primer = input$reverse_primer,
      default_assay_id = default_assay_id
    )
  }

  sync_summary_from_results <- function() {
    if (is.null(rv$summary) || !is.data.frame(rv$summary) || !nrow(rv$summary) || is.null(rv$results) || !length(rv$results)) return(invisible(NULL))
    for (nm in names(rv$results)) {
      r <- ensure_curation_state(rv$results[[nm]])
      r <- curation_rebuild(r, settings_for_result(r))
      rv$results[[nm]] <- r
      sm <- r$summary
      idx <- which(as.character(rv$summary$sample_id) == nm)
      if (!length(idx) || !nrow(sm)) next
      for (col in names(sm)) {
        if (!col %in% names(rv$summary)) rv$summary[[col]] <- NA
        rv$summary[idx[1], col] <- sm[1, col]
      }
    }
    invisible(NULL)
  }
