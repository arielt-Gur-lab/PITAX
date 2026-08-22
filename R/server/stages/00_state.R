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
      database=character(), hitlist_size=integer(), consensus_revision=integer(), status=character(), submitted_at=character(), last_checked_at=character(), stringsAsFactors=FALSE
    ),
    ncbi_last_contact = as.POSIXct(NA),
    blast_batch_status_text = "No batch operation has been run yet.",
    blast_raw = list(), blast_ids = data.frame(), blast_hits = data.frame(),
    taxonomy_summary = data.frame(), taxonomy_hits = data.frame(), taxonomy_counts = data.frame(),
    taxonomy_status_text = "No taxonomic analysis has been run yet.",
    taxonomy_batch_status_text = "No batch taxonomic analysis has been run yet.",
    project_status_text = "Current session has not been saved as a project.",
    project_loaded_name = "",
    context_peak_flag = NULL,
    pending_curation = NULL,
    auto_correct_preview_df = data.frame()
  )

  ensure_blast_jobs_schema <- function(df) {
    template <- data.frame(
      final_name=character(), original_name=character(), rid=character(), rtoe=character(),
      database=character(), hitlist_size=integer(), consensus_revision=integer(), status=character(), submitted_at=character(),
      last_checked_at=character(), stringsAsFactors=FALSE
    )
    if (!is.data.frame(df) || !nrow(df)) return(template)
    if (!"database" %in% names(df)) df$database <- ""
    if (!"hitlist_size" %in% names(df)) df$hitlist_size <- NA_integer_
    for (nm in setdiff(names(template), names(df))) df[[nm]] <- template[[nm]][NA_integer_]
    df <- df[, names(template), drop=FALSE]
    df$database <- as.character(df$database)
    df$hitlist_size <- suppressWarnings(as.integer(df$hitlist_size))
    df$consensus_revision <- suppressWarnings(as.integer(df$consensus_revision))
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
