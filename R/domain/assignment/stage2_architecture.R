# ============================================================
# PITAX v3.0.0-alpha.6
# Stage 2: Project -> Isolate -> Locus -> Read architecture
# ============================================================

stage2_scalar_text <- function(x, fallback = "") {
  if (is.null(x) || !length(x) || is.na(x[1])) return(fallback)
  value <- trimws(as.character(x[1]))
  if (!nzchar(value)) fallback else value
}

stage2_read_stem <- function(file_name) {
  sub("\\.[Aa][Bb]1$", "", basename(stage2_scalar_text(file_name)))
}

stage2_normalize_direction <- function(x, fallback = "Unknown") {
  value <- tolower(stage2_scalar_text(x, fallback))
  if (value %in% c("f", "for", "fwd", "forward")) return("Forward")
  if (value %in% c("r", "rev", "reverse")) return("Reverse")
  if (value %in% c("u", "unk", "unknown", "unknown / infer")) return("Unknown")
  fallback
}

stage2_direction_code <- function(x) {
  direction <- stage2_normalize_direction(x)
  if (direction == "Forward") return("F")
  if (direction == "Reverse") return("R")
  "U"
}

stage2_compose_read_name <- function(isolate, locus, direction) {
  clean_part <- function(x) {
    x <- gsub("\\s+", "-", stage2_scalar_text(x), perl = TRUE)
    x <- gsub("[^A-Za-z0-9_.-]", "-", x, perl = TRUE)
    x <- gsub("-+", "-", x, perl = TRUE)
    gsub("(^[-.]+|[-.]+$)", "", x, perl = TRUE)
  }
  isolate_part <- clean_part(isolate)
  locus_part <- clean_part(locus)
  direction_code <- stage2_direction_code(direction)
  if (!nzchar(isolate_part) || !nzchar(locus_part) || !direction_code %in% c("F", "R")) return("")
  paste(isolate_part, locus_part, direction_code, sep = "_")
}

stage2_initialize_read_assignment <- function(file_name, default_locus = "Unassigned", default_direction = "Unknown", forward_primer = "", reverse_primer = "", default_assay_id = "") {
  stem <- stage2_read_stem(file_name)
  direction <- stage2_normalize_direction(default_direction)
  locus <- stage2_scalar_text(default_locus, "Unassigned")
  primer <- if (direction == "Forward") stage2_scalar_text(forward_primer) else if (direction == "Reverse") stage2_scalar_text(reverse_primer) else ""

  data.frame(
    Read_ID = "",
    Source_ID = stem,
    File = basename(stage2_scalar_text(file_name)),
    Final_Name = "",
    Isolate = "",
    Assay_ID = stage2_scalar_text(default_assay_id),
    Locus = locus,
    Direction = direction,
    Primer = primer,
    Inference = "operator assignment required; locus/direction initialized from assay defaults",
    Notes = "",
    stringsAsFactors = FALSE
  )
}

stage2_make_read_assignments <- function(files, default_locus = "Unassigned", default_direction = "Unknown", forward_primer = "", reverse_primer = "", default_assay_id = "") {
  files <- as.character(files)
  if (!length(files)) return(stage2_empty_assignments())
  rows <- lapply(files, stage2_initialize_read_assignment,
                 default_locus = default_locus,
                 default_direction = default_direction,
                 forward_primer = forward_primer,
                 reverse_primer = reverse_primer,
                 default_assay_id = default_assay_id)
  out <- do.call(rbind, rows)
  out$Read_ID <- sprintf("read_%03d", seq_len(nrow(out)))
  rownames(out) <- NULL
  out
}

stage2_empty_assignments <- function() {
  data.frame(
    Read_ID = character(), Source_ID = character(), File = character(),
    Final_Name = character(), Isolate = character(), Assay_ID = character(), Locus = character(), Direction = character(),
    Primer = character(), Inference = character(), Notes = character(),
    stringsAsFactors = FALSE
  )
}

stage2_assignment_columns <- function() names(stage2_empty_assignments())

stage2_coerce_assignments <- function(df) {
  if (!is.data.frame(df)) return(stage2_empty_assignments())
  for (nm in setdiff(stage2_assignment_columns(), names(df))) df[[nm]] <- ""
  df <- df[, stage2_assignment_columns(), drop = FALSE]
  for (nm in names(df)) df[[nm]] <- as.character(df[[nm]])
  df$Direction <- vapply(df$Direction, stage2_normalize_direction, character(1))
  rownames(df) <- NULL
  df
}

stage2_sync_generated_names <- function(assignments, forward_primer = "", reverse_primer = "", assay_profiles = NULL) {
  assignments <- stage2_coerce_assignments(assignments)
  profiles <- if (exists("assay_coerce_profiles", mode = "function")) assay_coerce_profiles(assay_profiles) else data.frame()
  for (i in seq_len(nrow(assignments))) {
    if (nrow(profiles) && nzchar(assignments$Assay_ID[i])) {
      profile_idx <- match(assignments$Assay_ID[i], profiles$Assay_ID)
      if (!is.na(profile_idx)) {
        assignments$Locus[i] <- profiles$Locus_ID[profile_idx]
        if (assignments$Direction[i] == "Forward") assignments$Primer[i] <- profiles$Forward_Primer_Name[profile_idx]
        if (assignments$Direction[i] == "Reverse") assignments$Primer[i] <- profiles$Reverse_Primer_Name[profile_idx]
      }
    }
    assignments$Isolate[i] <- trimws(assignments$Isolate[i])
    assignments$Locus[i] <- trimws(assignments$Locus[i])
    assignments$Direction[i] <- stage2_normalize_direction(assignments$Direction[i])
    assignments$Final_Name[i] <- stage2_compose_read_name(assignments$Isolate[i], assignments$Locus[i], assignments$Direction[i])
    if (!nrow(profiles) && assignments$Direction[i] == "Forward" && nzchar(stage2_scalar_text(forward_primer))) assignments$Primer[i] <- stage2_scalar_text(forward_primer)
    if (!nrow(profiles) && assignments$Direction[i] == "Reverse" && nzchar(stage2_scalar_text(reverse_primer))) assignments$Primer[i] <- stage2_scalar_text(reverse_primer)
    assignments$Inference[i] <- if (nzchar(assignments$Final_Name[i])) "operator-assigned identity; final name generated by PITAX" else "operator assignment required"
  }
  assignments
}

stage2_identity_error <- function(assignments) {
  assignments <- stage2_coerce_assignments(assignments)
  if (!nrow(assignments)) return("No read assignments are available.")
  if (any(!nzchar(trimws(assignments$Isolate)))) return("Assign an isolate code to every read.")
  if (any(!nzchar(trimws(assignments$Locus)))) return("Assign a gene/locus to every read.")
  if (any(!assignments$Direction %in% c("Forward", "Reverse"))) return("Assign every read as Forward or Reverse.")
  generated <- vapply(seq_len(nrow(assignments)), function(i) stage2_compose_read_name(assignments$Isolate[i], assignments$Locus[i], assignments$Direction[i]), character(1))
  if (any(!nzchar(generated))) return("One or more final read names could not be generated.")
  if (anyDuplicated(generated)) return("Generated final read / FASTA names must be unique.")
  NULL
}

stage2_validate_assignments <- function(df, expected_source_ids = NULL, assay_profiles = NULL) {
  df <- stage2_coerce_assignments(df)
  if (!nrow(df)) return("No read assignments are available.")
  # Source/read-key integrity is the first safety boundary. Report collisions
  # before incomplete biological fields so duplicate AB1 basenames cannot be
  # hidden behind an unrelated pending-assignment message.
  if (anyDuplicated(df$Read_ID)) return("Read_ID values must be unique.")
  if (anyDuplicated(df$Source_ID)) return("Source_ID values must be unique; duplicate AB1 basenames are not safe.")
  required <- c("Read_ID", "Source_ID", "File", "Final_Name", "Isolate", "Locus", "Direction")
  for (nm in required) {
    if (any(!nzchar(trimws(df[[nm]])))) return(paste0(nm, " contains an empty value."))
  }
  if (any(!df$Direction %in% c("Forward", "Reverse", "Unknown"))) return("Direction must be Forward, Reverse, or Unknown.")
  if (!is.null(assay_profiles)) {
    profiles <- if (exists("assay_coerce_profiles", mode = "function")) assay_coerce_profiles(assay_profiles) else data.frame()
    profile_error <- if (exists("assay_validate_profiles", mode = "function")) assay_validate_profiles(profiles) else NULL
    if (!is.null(profile_error)) return(profile_error)
    if (any(!nzchar(trimws(df$Assay_ID)))) return("Assign an assay profile to every read.")
    if (any(!df$Assay_ID %in% profiles$Assay_ID)) return("One or more reads refer to an assay profile that does not exist.")
    expected_locus <- profiles$Locus_ID[match(df$Assay_ID, profiles$Assay_ID)]
    if (any(df$Locus != expected_locus)) return("Read locus must match the locus defined by its assay profile.")
  }
  if (!is.null(expected_source_ids)) {
    expected <- sort(unique(as.character(expected_source_ids)))
    observed <- sort(unique(df$Source_ID))
    if (!identical(expected, observed)) return("Read assignments do not match the currently uploaded AB1 files. Refresh the assignments.")
  }
  NULL
}

stage2_build_architecture <- function(assignments, project_id = "project", assay_profiles = NULL) {
  assignments <- stage2_coerce_assignments(assignments)
  error <- stage2_validate_assignments(assignments, assay_profiles = assay_profiles)
  if (!is.null(error)) stop(error, call. = FALSE)

  isolate_names <- unique(assignments$Isolate)
  isolates <- data.frame(
    Isolate_ID = sprintf("isolate_%03d", seq_along(isolate_names)),
    Isolate = isolate_names,
    stringsAsFactors = FALSE
  )
  assignments$Isolate_ID <- isolates$Isolate_ID[match(assignments$Isolate, isolates$Isolate)]

  locus_keys <- unique(paste(assignments$Isolate_ID, assignments$Locus, sep = "\r"))
  locus_parts <- strsplit(locus_keys, "\r", fixed = TRUE)
  loci <- data.frame(
    Locus_ID = sprintf("locus_%03d", seq_along(locus_keys)),
    Isolate_ID = vapply(locus_parts, `[`, character(1), 1),
    Locus = vapply(locus_parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
  assignments$Locus_ID <- loci$Locus_ID[match(paste(assignments$Isolate_ID, assignments$Locus, sep = "\r"), locus_keys)]

  reads <- assignments[, c("Read_ID", "Locus_ID", "Isolate_ID", "Assay_ID", "Source_ID", "File", "Final_Name", "Direction", "Primer", "Inference", "Notes"), drop = FALSE]
  list(
    schema = if (is.null(assay_profiles)) "pitax-project-architecture-v1" else "pitax-project-architecture-v2",
    project_id = stage2_scalar_text(project_id, "project"),
    assays = if (is.null(assay_profiles)) data.frame() else assay_coerce_profiles(assay_profiles),
    isolates = isolates,
    loci = loci,
    reads = reads
  )
}

stage2_architecture_summary <- function(architecture) {
  if (!is.list(architecture) || !is.data.frame(architecture$reads)) {
    return(data.frame(Isolates = 0L, Loci = 0L, Reads = 0L, Paired_loci = 0L, Single_read_loci = 0L, stringsAsFactors = FALSE))
  }
  reads <- architecture$reads
  by_locus <- split(reads, reads$Locus_ID)
  paired <- sum(vapply(by_locus, function(x) all(c("Forward", "Reverse") %in% x$Direction), logical(1)))
  singles <- sum(vapply(by_locus, nrow, integer(1)) == 1L)
  data.frame(
    Isolates = nrow(architecture$isolates),
    Loci = nrow(architecture$loci),
    Reads = nrow(reads),
    Paired_loci = paired,
    Single_read_loci = singles,
    stringsAsFactors = FALSE
  )
}

stage2_default_rename_map <- function(assignments) {
  assignments <- stage2_coerce_assignments(assignments)
  if (!nrow(assignments)) return(data.frame(Original_name = character(), New_name = character(), stringsAsFactors = FALSE))
  out <- ifelse(nzchar(trimws(assignments$Final_Name)), assignments$Final_Name, assignments$Source_ID)
  data.frame(Original_name = assignments$Source_ID, New_name = out, stringsAsFactors = FALSE)
}

stage2_migrate_v1_state <- function(state) {
  if (!is.list(state)) state <- list()
  result_ids <- if (is.list(state$results) && length(state$results)) names(state$results) else character()
  if (!length(result_ids) && is.data.frame(state$summary) && "sample_id" %in% names(state$summary)) result_ids <- as.character(state$summary$sample_id)
  if (!length(result_ids) && is.data.frame(state$rename) && "Original_name" %in% names(state$rename)) result_ids <- as.character(state$rename$Original_name)
  if (!length(result_ids)) {
    state$read_assignments <- stage2_empty_assignments()
    state$architecture <- NULL
    state$migration_log <- "Schema 1 contained no reads; initialized an empty Stage 2 architecture."
    return(state)
  }

  settings <- if (is.list(state$settings)) state$settings else list()
  locus <- stage2_scalar_text(settings$target, "Unassigned")
  direction <- stage2_normalize_direction(settings$sequencing_primer)
  primer <- if (direction == "Forward") stage2_scalar_text(settings$forward_primer) else if (direction == "Reverse") stage2_scalar_text(settings$reverse_primer) else ""
  isolate_names <- result_ids
  if (is.data.frame(state$rename) && all(c("Original_name", "New_name") %in% names(state$rename))) {
    idx <- match(result_ids, as.character(state$rename$Original_name))
    resolved <- as.character(state$rename$New_name[idx])
    good <- !is.na(resolved) & nzchar(trimws(resolved))
    isolate_names[good] <- resolved[good]
  }
  assignments <- data.frame(
    Read_ID = sprintf("read_%03d", seq_along(result_ids)),
    Source_ID = result_ids,
    File = paste0(result_ids, ".ab1"),
    Final_Name = isolate_names,
    Isolate = isolate_names,
    Locus = rep(locus, length(result_ids)),
    Direction = rep(direction, length(result_ids)),
    Primer = rep(primer, length(result_ids)),
    Inference = rep("migrated from project schema 1", length(result_ids)),
    Notes = rep("Original AB1 filename was not stored in schema 1.", length(result_ids)),
    stringsAsFactors = FALSE
  )
  state$read_assignments <- assignments
  state$architecture <- stage2_build_architecture(assignments)
  state$migration_log <- paste0("Migrated ", length(result_ids), " read(s) from project schema 1 without changing result, rename, QC, BLAST, or taxonomy records.")
  state
}

stage2_migrate_v2_state <- function(state) {
  if (!is.list(state)) state <- list()
  assignments <- stage2_coerce_assignments(state$read_assignments)
  if (nrow(assignments)) {
    for (i in seq_len(nrow(assignments))) {
      if (!nzchar(trimws(assignments$Final_Name[i]))) {
        j <- if (is.data.frame(state$rename)) match(assignments$Source_ID[i], as.character(state$rename$Original_name)) else NA_integer_
        assignments$Final_Name[i] <- if (!is.na(j)) as.character(state$rename$New_name[j]) else assignments$Source_ID[i]
      }
      assignments$Inference[i] <- paste0("migrated from project schema 2; review explicit isolate/locus/direction fields. Previous: ", assignments$Inference[i])
    }
  }
  state$read_assignments <- assignments
  state$architecture <- if (nrow(assignments) && is.null(stage2_validate_assignments(assignments))) stage2_build_architecture(assignments) else state$architecture
  state$migration_log <- "Migrated project schema 2 to schema 3. Existing evidence and output names were preserved; explicit isolate, locus and direction fields should be reviewed before any new trimming run."
  state
}
