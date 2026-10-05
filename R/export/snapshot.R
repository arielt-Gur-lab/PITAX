# ============================================================
# Frozen export snapshot
# ============================================================

pitax_export_empty_upload_files <- function() {
  data.frame(
    name = character(), size = numeric(), datapath = character(),
    Source_ID = character(), Bytes_available = logical(),
    stringsAsFactors = FALSE
  )
}

pitax_export_capture_snapshot <- function(state, upload_files = NULL, pending_curation = NULL,
                                         auto_correct_preview = NULL, current_sample = NULL,
                                         source_stage = "export", project_mode = "simple",
                                         app_version = "unknown", schema_version = 6L) {
  upload_files <- if (is.null(upload_files) || !is.data.frame(upload_files) || !nrow(upload_files)) {
    pitax_export_empty_upload_files()
  } else {
    df <- upload_files
    if (!"Source_ID" %in% names(df)) {
      df$Source_ID <- vapply(as.character(df$name), function(nm) {
        tryCatch(stage2_read_stem(nm), error = function(e) tools::file_path_sans_ext(basename(nm)))
      }, character(1))
    }
    if (!"Bytes_available" %in% names(df)) {
      df$Bytes_available <- file.exists(as.character(df$datapath))
    }
    df
  }

  list(
    snapshot_id = paste0("snap_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", as.integer(runif(1, 1e5, 1e6 - 1))),
    captured_at = Sys.time(),
    timezone = format(Sys.time(), "%Z"),
    app_version = as.character(app_version)[1],
    schema_version = suppressWarnings(as.integer(schema_version)[1]),
    source_stage = as.character(source_stage)[1],
    project_mode = as.character(project_mode)[1],
    current_sample = if (is.null(current_sample)) NA_character_ else as.character(current_sample)[1],
    pending_curation = pending_curation,
    auto_correct_preview = auto_correct_preview,
    upload_files = upload_files,
    results = if (is.list(state$results)) state$results else list(),
    summary = if (is.data.frame(state$summary)) state$summary else data.frame(),
    rename = if (is.data.frame(state$rename)) state$rename else data.frame(),
    settings = if (is.list(state$settings)) state$settings else list(),
    assay_profiles = if (is.data.frame(state$assay_profiles)) state$assay_profiles else data.frame(),
    project_defaults = if (is.list(state$project_defaults)) state$project_defaults else list(),
    read_assignments = if (is.data.frame(state$read_assignments)) state$read_assignments else data.frame(),
    architecture = if (is.list(state$architecture)) state$architecture else list(),
    consensus_set = if (is.list(state$consensus_set)) state$consensus_set else list(),
    blast_jobs = if (is.data.frame(state$blast_jobs)) state$blast_jobs else data.frame(),
    blast_hits = if (is.data.frame(state$blast_hits)) state$blast_hits else data.frame(),
    blast_ids = if (is.data.frame(state$blast_ids)) state$blast_ids else data.frame(),
    blast_raw = if (is.list(state$blast_raw)) state$blast_raw else list(),
    taxonomy_summary = if (is.data.frame(state$taxonomy_summary)) state$taxonomy_summary else data.frame(),
    taxonomy_hits = if (is.data.frame(state$taxonomy_hits)) state$taxonomy_hits else data.frame(),
    taxonomy_counts = if (is.data.frame(state$taxonomy_counts)) state$taxonomy_counts else data.frame(),
    multilocus_profile = if (is.list(state$multilocus_profile)) state$multilocus_profile else list()
  )
}

pitax_export_snapshot_from_rv <- function(rv, input = NULL, source_stage = "export",
                                          app_version = NULL, schema_version = NULL) {
  upload_files <- NULL
  if (!is.null(input) && !is.null(input$ab1_files) && is.data.frame(input$ab1_files) && nrow(input$ab1_files)) {
    upload_files <- input$ab1_files
  }
  current_sample <- NULL
  if (!is.null(input)) {
    for (nm in c("inspect_sample", "blast_sample", "tax_sample", "consensus_selected_id")) {
      if (!is.null(input[[nm]]) && length(input[[nm]]) && nzchar(as.character(input[[nm]])[1])) {
        current_sample <- as.character(input[[nm]])[1]
        break
      }
    }
  }
  if (is.null(app_version)) {
    app_version <- tryCatch(get("APP_VERSION", inherits = TRUE), error = function(e) "unknown")
  }
  if (is.null(schema_version)) {
    schema_version <- tryCatch(get("PROJECT_SCHEMA_VERSION", inherits = TRUE), error = function(e) 6L)
  }
  state <- list(
    results = rv$results,
    summary = rv$summary,
    rename = rv$rename,
    settings = rv$settings,
    assay_profiles = rv$assay_profiles,
    project_defaults = rv$project_defaults,
    read_assignments = rv$read_assignments,
    architecture = rv$architecture,
    consensus_set = rv$consensus_set,
    blast_jobs = rv$blast_jobs,
    blast_hits = rv$blast_hits,
    blast_ids = rv$blast_ids,
    blast_raw = rv$blast_raw,
    taxonomy_summary = rv$taxonomy_summary,
    taxonomy_hits = rv$taxonomy_hits,
    taxonomy_counts = rv$taxonomy_counts,
    multilocus_profile = rv$multilocus_profile
  )
  pitax_export_capture_snapshot(
    state = state,
    upload_files = upload_files,
    pending_curation = rv$pending_curation,
    auto_correct_preview = rv$auto_correct_preview_df,
    current_sample = current_sample,
    source_stage = source_stage,
    project_mode = if (is.null(rv$project_mode)) "simple" else rv$project_mode,
    app_version = app_version,
    schema_version = schema_version
  )
}

pitax_export_named_list_to_df <- function(x, id_col = "Setting") {
  if (!is.list(x) || !length(x)) {
    return(data.frame(Setting = character(), Value = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    Setting = names(x),
    Value = vapply(x, function(v) paste(as.character(unlist(v)), collapse = "; "), character(1)),
    stringsAsFactors = FALSE
  )
}
