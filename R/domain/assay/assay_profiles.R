# ============================================================
# PITAX v3.0.0-alpha.10.1
# Alpha 10: controlled loci, assay profiles and schema-5 bridge
# ============================================================

assay_scalar_text <- function(x, fallback = "") {
  if (is.null(x) || !length(x) || is.na(x[1])) return(fallback)
  value <- trimws(as.character(x[1]))
  if (!nzchar(value)) fallback else value
}

pitax_locus_vocabulary <- function() {
  data.frame(
    Locus_ID = c("ITS", "LSU", "TEF1", "RPB2", "TUB2", "CYP51", "SDHB", "IGS"),
    Display_Name = c(
      "ITS", "LSU", "TEF1 / EF1-alpha", "RPB2", "Beta-tubulin / TUB2",
      "CYP51", "SDHB", "IGS"
    ),
    stringsAsFactors = FALSE
  )
}

pitax_locus_choices <- function() {
  vocabulary <- pitax_locus_vocabulary()
  stats::setNames(vocabulary$Locus_ID, vocabulary$Display_Name)
}

pitax_normalize_locus_id <- function(x, fallback = "") {
  value <- toupper(trimws(assay_scalar_text(x)))
  aliases <- c(
    "ITS" = "ITS",
    "LSU" = "LSU",
    "TEF1" = "TEF1",
    "TEF1 / EF1-ALPHA" = "TEF1",
    "TEF1 / EF1-\u03B1" = "TEF1",
    "EF1-ALPHA" = "TEF1",
    "EF1-\u03B1" = "TEF1",
    "RPB2" = "RPB2",
    "TUB2" = "TUB2",
    "BETA-TUBULIN" = "TUB2",
    "BETA-TUBULIN / TUB2" = "TUB2",
    "\u03B2-TUBULIN" = "TUB2",
    "CYP51" = "CYP51",
    "SDHB" = "SDHB",
    "IGS" = "IGS"
  )
  if (value %in% names(aliases)) unname(aliases[[value]]) else fallback
}

pitax_locus_display_name <- function(locus_id, fallback = "") {
  locus_id <- pitax_normalize_locus_id(locus_id)
  vocabulary <- pitax_locus_vocabulary()
  idx <- match(locus_id, vocabulary$Locus_ID)
  if (is.na(idx)) fallback else vocabulary$Display_Name[idx]
}

assay_clean_sequence <- function(x) {
  toupper(gsub("[^ACGTRYSWKMBDHVN]", "", assay_scalar_text(x), perl = TRUE))
}

assay_empty_profiles <- function() {
  data.frame(
    Assay_ID = character(),
    Assay_Name = character(),
    Locus_ID = character(),
    Locus_Display_Name = character(),
    Forward_Primer_Name = character(),
    Forward_Primer_Sequence = character(),
    Reverse_Primer_Name = character(),
    Reverse_Primer_Sequence = character(),
    Expected_Amplicon_Length = integer(),
    Maximum_Sequence_Position = integer(),
    stringsAsFactors = FALSE
  )
}

assay_profile_columns <- function() names(assay_empty_profiles())

assay_make_id <- function(locus_id, existing_ids = character()) {
  stem <- tolower(pitax_normalize_locus_id(locus_id, "assay"))
  stem <- gsub("[^a-z0-9]+", "-", stem, perl = TRUE)
  stem <- gsub("(^-+|-+$)", "", stem, perl = TRUE)
  if (!nzchar(stem)) stem <- "assay"
  candidate <- paste0("assay-", stem)
  suffix <- 2L
  while (candidate %in% existing_ids) {
    candidate <- paste0("assay-", stem, "-", suffix)
    suffix <- suffix + 1L
  }
  candidate
}

assay_coerce_profiles <- function(df) {
  if (!is.data.frame(df)) return(assay_empty_profiles())
  template <- assay_empty_profiles()
  for (nm in setdiff(names(template), names(df))) df[[nm]] <- template[[nm]][NA_integer_]
  df <- df[, names(template), drop = FALSE]
  text_columns <- setdiff(names(template), c("Expected_Amplicon_Length", "Maximum_Sequence_Position"))
  for (nm in text_columns) df[[nm]] <- as.character(df[[nm]])
  df$Locus_ID <- vapply(df$Locus_ID, pitax_normalize_locus_id, character(1))
  df$Locus_Display_Name <- vapply(seq_len(nrow(df)), function(i) {
    pitax_locus_display_name(df$Locus_ID[i], assay_scalar_text(df$Locus_Display_Name[i]))
  }, character(1))
  df$Forward_Primer_Sequence <- vapply(df$Forward_Primer_Sequence, assay_clean_sequence, character(1))
  df$Reverse_Primer_Sequence <- vapply(df$Reverse_Primer_Sequence, assay_clean_sequence, character(1))
  df$Expected_Amplicon_Length <- suppressWarnings(as.integer(df$Expected_Amplicon_Length))
  df$Maximum_Sequence_Position <- suppressWarnings(as.integer(df$Maximum_Sequence_Position))
  rownames(df) <- NULL
  df
}

assay_profile_from_legacy_settings <- function(settings = NULL, assay_id = NULL) {
  settings <- if (is.list(settings)) settings else list()
  locus_id <- pitax_normalize_locus_id(settings$target, "ITS")
  display_name <- pitax_locus_display_name(locus_id, locus_id)
  assay_id <- assay_scalar_text(assay_id, assay_make_id(locus_id))
  assay_name <- assay_scalar_text(settings$assay_name, display_name)
  data.frame(
    Assay_ID = assay_id,
    Assay_Name = assay_name,
    Locus_ID = locus_id,
    Locus_Display_Name = display_name,
    Forward_Primer_Name = assay_scalar_text(settings$forward_primer),
    Forward_Primer_Sequence = assay_clean_sequence(settings$forward_primer_seq),
    Reverse_Primer_Name = assay_scalar_text(settings$reverse_primer),
    Reverse_Primer_Sequence = assay_clean_sequence(settings$reverse_primer_seq),
    Expected_Amplicon_Length = suppressWarnings(as.integer(if (is.null(settings$expected_amplicon_len)) 650L else settings$expected_amplicon_len[1])),
    Maximum_Sequence_Position = suppressWarnings(as.integer(if (is.null(settings$absolute_max_base_index)) 680L else settings$absolute_max_base_index[1])),
    stringsAsFactors = FALSE
  )
}

assay_default_profiles <- function() assay_profile_from_legacy_settings(list(target = "ITS"))

assay_validate_profiles <- function(df) {
  df <- assay_coerce_profiles(df)
  if (!nrow(df)) return("At least one assay profile is required.")
  if (any(!nzchar(trimws(df$Assay_ID)))) return("Every assay requires an Assay_ID.")
  if (anyDuplicated(df$Assay_ID)) return("Assay_ID values must be unique.")
  if (any(!nzchar(trimws(df$Assay_Name)))) return("Every assay requires a display name.")
  valid_loci <- pitax_locus_vocabulary()$Locus_ID
  if (any(!df$Locus_ID %in% valid_loci)) return("Every assay must use a locus from the controlled PITAX vocabulary.")
  if (any(is.na(df$Expected_Amplicon_Length) | df$Expected_Amplicon_Length < 1L)) return("Expected amplicon length must be a positive integer for every assay.")
  if (any(is.na(df$Maximum_Sequence_Position) | df$Maximum_Sequence_Position < 50L)) return("Maximum sequence position must be at least 50 bp for every assay.")
  NULL
}

assay_project_defaults_from_legacy_settings <- function(settings = NULL) {
  settings <- if (is.list(settings)) settings else list()
  defaults <- list(
    window = 25L,
    min_peak_ratio = 3,
    min_relative_signal = 0.20,
    min_len_before_collapse = 350L,
    bad_run_windows = 12L,
    min_usable_len = 400L,
    enable_primer_mapping = FALSE,
    ambiguous_peak_strong_ratio = 1.25,
    ambiguous_peak_moderate_ratio = 1.75,
    ambiguous_peak_min_relative_signal = 0.20,
    auto_correct_min_alt_to_called = 1.80,
    auto_correct_min_alt_to_third = 2.00,
    auto_correct_max_peak_offset = 2L,
    auto_correct_min_relative_signal = 0.50
  )
  for (nm in names(defaults)) if (!is.null(settings[[nm]]) && length(settings[[nm]])) defaults[[nm]] <- settings[[nm]][1]
  defaults
}

assay_resolve_read_settings <- function(profile, project_defaults = NULL, direction = "Unknown") {
  profiles <- assay_coerce_profiles(profile)
  if (nrow(profiles) != 1L) stop("Exactly one assay profile is required to resolve read settings.", call. = FALSE)
  validation_error <- assay_validate_profiles(profiles)
  if (!is.null(validation_error)) stop(validation_error, call. = FALSE)
  settings <- assay_project_defaults_from_legacy_settings(project_defaults)
  settings$target <- profiles$Locus_ID[1]
  settings$assay_id <- profiles$Assay_ID[1]
  settings$assay_name <- profiles$Assay_Name[1]
  settings$forward_primer <- profiles$Forward_Primer_Name[1]
  settings$forward_primer_seq <- profiles$Forward_Primer_Sequence[1]
  settings$reverse_primer <- profiles$Reverse_Primer_Name[1]
  settings$reverse_primer_seq <- profiles$Reverse_Primer_Sequence[1]
  settings$expected_amplicon_len <- profiles$Expected_Amplicon_Length[1]
  settings$absolute_max_base_index <- profiles$Maximum_Sequence_Position[1]
  settings$sequencing_primer <- if (exists("stage2_normalize_direction", mode = "function")) stage2_normalize_direction(direction) else assay_scalar_text(direction, "Unknown")
  settings
}

assay_migrate_schema5_state <- function(state) {
  if (!is.list(state)) state <- list()
  legacy_settings <- if (is.list(state$settings)) state$settings else list()
  profiles <- assay_profile_from_legacy_settings(legacy_settings)
  defaults <- assay_project_defaults_from_legacy_settings(legacy_settings)
  assignments <- if (exists("stage2_coerce_assignments", mode = "function")) stage2_coerce_assignments(state$read_assignments) else state$read_assignments
  if (is.data.frame(assignments) && nrow(assignments)) {
    assignments$Assay_ID <- profiles$Assay_ID[1]
    assignments$Locus <- profiles$Locus_ID[1]
    if (exists("stage2_coerce_assignments", mode = "function")) assignments <- stage2_coerce_assignments(assignments)
  }
  state$assay_profiles <- profiles
  state$project_defaults <- defaults
  state$settings <- assay_resolve_read_settings(
    profiles,
    defaults,
    assay_scalar_text(legacy_settings$sequencing_primer, "Unknown")
  )
  state$read_assignments <- assignments
  if (is.data.frame(assignments) && nrow(assignments) && exists("stage2_build_architecture", mode = "function")) {
    state$architecture <- stage2_build_architecture(assignments, assay_profiles = profiles)
  }
  previous_log <- assay_scalar_text(state$migration_log)
  bridge_log <- paste0(
    "Migrated project schema 5 to schema 6: the legacy run-level settings became assay profile ",
    profiles$Assay_ID[1], ", and every existing read was linked to that assay."
  )
  state$migration_log <- paste(c(previous_log, bridge_log)[nzchar(c(previous_log, bridge_log))], collapse = " ")
  state
}
