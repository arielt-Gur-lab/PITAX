# ============================================================
# PITAX Stage 4: multi-locus isolate profile foundation
# ============================================================

stage4_scalar_text <- function(x, default = "") {
  if (is.null(x) || !length(x) || is.na(x[1])) return(default)
  out <- trimws(as.character(x[1]))
  if (!nzchar(out)) default else out
}

stage4_empty_sources <- function() {
  data.frame(
    Source = character(), Source_MD5 = character(), Project_App_Version = character(),
    Project_Saved_At = character(), Project_Schema = integer(), Locus = character(),
    Isolates = integer(), Analysis_Sequences = integer(), Taxonomy_Complete = integer(),
    stringsAsFactors = FALSE
  )
}

stage4_empty_evidence <- function() {
  data.frame(
    Source = character(), Source_MD5 = character(), Project_App_Version = character(),
    Project_Saved_At = character(), Project_Schema = integer(), Consensus_ID = character(),
    Final_Name = character(), Isolate = character(), Locus = character(),
    Analysis_Status = character(), Sequence_Length = integer(), Consensus_Revision = integer(),
    Sequence = character(), Taxonomy_Status = character(), Recommended_Identification = character(),
    Recommended_Level = character(), Confidence = character(), Supported_Genus = character(),
    Supported_Species = character(), Best_Molecular_Match = character(),
    Best_Match_Accession = character(), Best_Match_Identity = numeric(),
    Best_Match_Coverage = numeric(), Reference_Support = character(),
    Locus_Discrimination = character(), RID = character(), Taxonomy_Analyzed_At = character(),
    stringsAsFactors = FALSE
  )
}

stage4_empty_profile_summary <- function() {
  data.frame(
    Isolate = character(), Loci = character(), Locus_Count = integer(),
    Taxonomy_Complete = character(), Profile_Status = character(),
    Profile_Conclusion = character(), Supported_Genus = character(),
    Supported_Species = character(), Conflicting_Calls = character(),
    Reference_Context = character(), Next_Action = character(),
    stringsAsFactors = FALSE
  )
}

stage4_empty_profile <- function() {
  list(
    schema = "pitax-multilocus-profile-v1",
    algorithm = "pitax-evidence-preserving-profile-v1",
    built_at = "",
    sources = stage4_empty_sources(),
    evidence = stage4_empty_evidence(),
    profiles = stage4_empty_profile_summary(),
    warnings = character()
  )
}

stage4_df_value <- function(df, column, default = "") {
  if (!is.data.frame(df) || !nrow(df) || !column %in% names(df)) return(default)
  stage4_scalar_text(df[[column]], default)
}

stage4_df_number <- function(df, column) {
  if (!is.data.frame(df) || !nrow(df) || !column %in% names(df)) return(NA_real_)
  out <- suppressWarnings(as.numeric(df[[column]][1]))
  if (length(out) != 1L || !is.finite(out)) NA_real_ else out
}

stage4_taxon_parts <- function(taxonomy_row) {
  rank <- tolower(stage4_df_value(taxonomy_row, "recommended_level"))
  identification <- stage4_df_value(taxonomy_row, "recommended_identification")
  if (identification %in% c("Not analyzed", "Unresolved", "Insufficient evidence")) identification <- ""

  genus <- ""
  species <- ""
  if (rank == "species" && nzchar(identification)) {
    species <- identification
    genus <- stage4_df_value(taxonomy_row, "best_match_genus")
    if (!nzchar(genus)) genus <- strsplit(identification, "\\s+")[[1]][1]
  } else if (rank == "genus" && nzchar(identification)) {
    genus <- identification
  }
  list(rank = rank, identification = identification, genus = genus, species = species)
}

stage4_project_taxonomy_row <- function(taxonomy_summary, consensus_id, final_name) {
  if (!is.data.frame(taxonomy_summary) || !nrow(taxonomy_summary)) return(data.frame())
  idx <- integer()
  if ("original_name" %in% names(taxonomy_summary)) {
    idx <- which(as.character(taxonomy_summary$original_name) == consensus_id)
  }
  if (!length(idx) && "final_name" %in% names(taxonomy_summary)) {
    idx <- which(as.character(taxonomy_summary$final_name) == final_name)
  }
  if (!length(idx)) return(data.frame())
  taxonomy_summary[idx[length(idx)], , drop = FALSE]
}

stage4_validate_project <- function(project, source_name = "project") {
  if (!is.list(project) || !identical(project$format, "SangerSequencePipelineProject") || !is.list(project$state)) {
    stop(paste0(source_name, " is not a valid PITAX project."), call. = FALSE)
  }
  schema <- suppressWarnings(as.integer(project$schema_version))
  if (length(schema) != 1L || is.na(schema) || schema < 4L) {
    stop(paste0(source_name, " predates the Stage 3 analysis-sequence schema. Load and resave it in the current PITAX version first."), call. = FALSE)
  }
  consensus_set <- project$state$consensus_set
  results <- project$state$results
  gate_error <- stage3_consensus_gate_error(consensus_set, results)
  if (!is.null(gate_error)) stop(paste0(source_name, ": ", gate_error), call. = FALSE)

  records <- consensus_set$records
  loci <- unique(trimws(vapply(records, function(x) stage4_scalar_text(x$locus), character(1))))
  loci <- loci[nzchar(loci)]
  if (length(loci) != 1L) {
    stop(paste0(source_name, " must contain exactly one Gene/Locus; found ", length(loci), "."), call. = FALSE)
  }
  invisible(TRUE)
}

stage4_extract_project_evidence <- function(project, source_name, source_md5 = "") {
  stage4_validate_project(project, source_name)
  st <- project$state
  records <- st$consensus_set$records
  taxonomy_summary <- if (is.data.frame(st$taxonomy_summary)) st$taxonomy_summary else data.frame()
  rows <- vector("list", length(records))
  ids <- names(records)

  for (i in seq_along(records)) {
    record <- stage3_ensure_record_curation(records[[i]])
    consensus_id <- stage4_scalar_text(record$consensus_id, ids[i])
    final_name <- stage4_scalar_text(record$final_name, consensus_id)
    tax <- stage4_project_taxonomy_row(taxonomy_summary, consensus_id, final_name)
    parts <- stage4_taxon_parts(tax)
    seq_text <- stage4_scalar_text(record$sequence)
    rows[[i]] <- data.frame(
      Source = source_name,
      Source_MD5 = source_md5,
      Project_App_Version = stage4_scalar_text(project$app_version, "unknown"),
      Project_Saved_At = stage4_scalar_text(project$saved_at),
      Project_Schema = suppressWarnings(as.integer(project$schema_version)),
      Consensus_ID = consensus_id,
      Final_Name = final_name,
      Isolate = stage4_scalar_text(record$isolate),
      Locus = stage4_scalar_text(record$locus),
      Analysis_Status = stage4_scalar_text(record$status),
      Sequence_Length = nchar(seq_text),
      Consensus_Revision = if (is.list(record$curation)) suppressWarnings(as.integer(record$curation$revision)) else 0L,
      Sequence = seq_text,
      Taxonomy_Status = if (nrow(tax)) "Analyzed" else "Not analyzed",
      Recommended_Identification = parts$identification,
      Recommended_Level = parts$rank,
      Confidence = stage4_df_value(tax, "confidence"),
      Supported_Genus = parts$genus,
      Supported_Species = parts$species,
      Best_Molecular_Match = stage4_df_value(tax, "best_molecular_match"),
      Best_Match_Accession = stage4_df_value(tax, "best_match_accession"),
      Best_Match_Identity = stage4_df_number(tax, "best_match_identity_percent"),
      Best_Match_Coverage = stage4_df_number(tax, "best_match_query_coverage_percent"),
      Reference_Support = stage4_df_value(tax, "reference_support"),
      Locus_Discrimination = stage4_df_value(tax, "locus_discrimination"),
      RID = stage4_df_value(tax, "rid"),
      Taxonomy_Analyzed_At = stage4_df_value(tax, "analyzed_at"),
      stringsAsFactors = FALSE
    )
  }
  out <- if (length(rows)) do.call(rbind, rows) else stage4_empty_evidence()
  rownames(out) <- NULL
  out
}

stage4_profile_row <- function(evidence) {
  loci <- sort(unique(as.character(evidence$Locus)))
  analyzed <- evidence$Taxonomy_Status == "Analyzed"
  supported <- analyzed & nzchar(as.character(evidence$Supported_Genus))
  genera <- unique(as.character(evidence$Supported_Genus[supported]))
  genera <- genera[nzchar(genera)]
  species <- unique(as.character(evidence$Supported_Species[nzchar(as.character(evidence$Supported_Species))]))
  species <- species[nzchar(species)]
  calls <- unique(as.character(evidence$Recommended_Identification[nzchar(as.character(evidence$Recommended_Identification))]))
  calls <- calls[nzchar(calls)]
  reference_context <- unique(as.character(evidence$Reference_Support[nzchar(as.character(evidence$Reference_Support))]))
  reference_context <- reference_context[nzchar(reference_context)]

  status <- "PARTIAL_EVIDENCE"
  conclusion <- "Available loci do not yet support a combined taxonomic conclusion."
  next_action <- "Complete missing locus-level BLAST and taxonomic interpretation; retain every locus as separate evidence."

  if (length(loci) < 2L) {
    status <- "SINGLE_LOCUS"
    conclusion <- "Only one locus is present; this is not a multi-locus profile."
    next_action <- "Add a separately processed second locus for the same isolate."
  } else if (!any(analyzed)) {
    status <- "NO_TAXONOMY"
    conclusion <- "Sequences are present, but no locus has a completed taxonomic interpretation."
    next_action <- "Complete BLAST retrieval and taxonomic interpretation for the included loci."
  } else if (length(genera) > 1L) {
    status <- "GENUS_CONFLICT"
    conclusion <- paste0("Conflicting supported genera across loci: ", paste(sort(genera), collapse = ", "), ". No combined identification is reported.")
    next_action <- "Review isolate assignment, chromatograms, sequence curation and reference quality; repeat the affected locus if needed."
  } else if (length(species) > 1L) {
    status <- "SPECIES_CONFLICT"
    conclusion <- paste0("Loci support different species-level calls: ", paste(sort(species), collapse = ", "), ". The conflict is retained rather than voted away.")
    next_action <- "Review locus discrimination and reference provenance; use a taxon-informed marker or curated reference set before reporting species."
  } else if (length(species) == 1L && sum(as.character(evidence$Supported_Species) == species[1]) >= 2L) {
    status <- "CONCORDANT_SPECIES"
    conclusion <- paste0("Concordant species-level evidence across at least two loci: ", species[1], ".")
    next_action <- "Confirm that decisive matches use appropriate curated or type-linked references before final reporting."
  } else if (length(genera) == 1L && sum(as.character(evidence$Supported_Genus) == genera[1]) >= 2L) {
    status <- "CONCORDANT_GENUS"
    conclusion <- paste0("The loci are concordant at genus level: ", genera[1], "; species remains unresolved or locus-dependent.")
    next_action <- "Add or review a taxon-informed secondary marker; PITAX does not infer a species by flat locus voting."
  } else if (sum(analyzed) < nrow(evidence)) {
    conclusion <- "Some loci have taxonomic interpretation and others are still missing; no combined call is made."
  }

  data.frame(
    Isolate = as.character(evidence$Isolate[1]),
    Loci = paste(loci, collapse = ", "),
    Locus_Count = length(loci),
    Taxonomy_Complete = paste0(sum(analyzed), "/", nrow(evidence)),
    Profile_Status = status,
    Profile_Conclusion = conclusion,
    Supported_Genus = if (length(genera) == 1L) genera else "",
    Supported_Species = if (identical(status, "CONCORDANT_SPECIES")) species else "",
    Conflicting_Calls = if (status %in% c("GENUS_CONFLICT", "SPECIES_CONFLICT")) paste(sort(calls), collapse = " | ") else "",
    Reference_Context = paste(reference_context, collapse = " | "),
    Next_Action = next_action,
    stringsAsFactors = FALSE
  )
}

stage4_build_profile <- function(projects, source_names = names(projects), source_md5 = rep("", length(projects))) {
  if (!is.list(projects) || !length(projects)) stop("Add at least one completed PITAX project.", call. = FALSE)
  if (is.null(source_names) || length(source_names) != length(projects)) source_names <- paste0("project_", seq_along(projects))
  if (length(source_md5) != length(projects)) source_md5 <- rep("", length(projects))
  source_names <- make.unique(vapply(source_names, stage4_scalar_text, character(1), default = "project"))
  md5_present <- nzchar(source_md5)
  if (anyDuplicated(source_md5[md5_present])) stop("The same project file was added more than once.", call. = FALSE)

  parts <- Map(stage4_extract_project_evidence, projects, source_names, source_md5)
  evidence <- do.call(rbind, parts)
  rownames(evidence) <- NULL
  if (!nrow(evidence)) stop("The selected projects contain no analysis sequences.", call. = FALSE)
  if (any(!nzchar(evidence$Isolate)) || any(!nzchar(evidence$Locus))) {
    stop("Every imported analysis sequence must retain explicit Isolate and Locus fields.", call. = FALSE)
  }

  evidence_key <- paste(tolower(trimws(evidence$Isolate)), tolower(trimws(evidence$Locus)), sep = "\r")
  if (anyDuplicated(evidence_key)) {
    duplicate_rows <- which(duplicated(evidence_key) | duplicated(evidence_key, fromLast = TRUE))
    duplicate_labels <- unique(paste0(evidence$Isolate[duplicate_rows], " / ", evidence$Locus[duplicate_rows]))
    stop(paste0("Duplicate Isolate/Locus evidence is not allowed: ", paste(duplicate_labels, collapse = ", "), "."), call. = FALSE)
  }

  sources <- lapply(seq_along(parts), function(i) {
    ev <- parts[[i]]
    data.frame(
      Source = source_names[i], Source_MD5 = source_md5[i],
      Project_App_Version = stage4_scalar_text(projects[[i]]$app_version, "unknown"),
      Project_Saved_At = stage4_scalar_text(projects[[i]]$saved_at),
      Project_Schema = suppressWarnings(as.integer(projects[[i]]$schema_version)),
      Locus = paste(sort(unique(ev$Locus)), collapse = ", "),
      Isolates = length(unique(tolower(ev$Isolate))), Analysis_Sequences = nrow(ev),
      Taxonomy_Complete = sum(ev$Taxonomy_Status == "Analyzed"), stringsAsFactors = FALSE
    )
  })
  sources <- do.call(rbind, sources)
  rownames(sources) <- NULL

  isolate_key <- tolower(trimws(evidence$Isolate))
  groups <- split(seq_len(nrow(evidence)), isolate_key)
  profiles <- do.call(rbind, lapply(groups, function(idx) stage4_profile_row(evidence[idx, , drop = FALSE])))
  rownames(profiles) <- NULL

  warnings <- character()
  isolate_variants <- lapply(groups, function(idx) unique(as.character(evidence$Isolate[idx])))
  case_variants <- isolate_variants[vapply(isolate_variants, length, integer(1)) > 1L]
  if (length(case_variants)) warnings <- c(warnings, "Some isolate codes differed only by letter case and were combined; review the evidence table.")

  list(
    schema = "pitax-multilocus-profile-v1",
    algorithm = "pitax-evidence-preserving-profile-v1",
    built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    sources = sources, evidence = evidence, profiles = profiles, warnings = warnings
  )
}

stage4_ensure_profile <- function(profile) {
  if (!is.list(profile) || !identical(stage4_scalar_text(profile$schema), "pitax-multilocus-profile-v1")) return(stage4_empty_profile())
  if (!is.data.frame(profile$sources)) profile$sources <- stage4_empty_sources()
  if (!is.data.frame(profile$evidence)) profile$evidence <- stage4_empty_evidence()
  if (!is.data.frame(profile$profiles)) profile$profiles <- stage4_empty_profile_summary()
  if (is.null(profile$warnings)) profile$warnings <- character()
  profile
}

stage4_isolate_profile <- function(profile, isolate) {
  profile <- stage4_ensure_profile(profile)
  isolate <- stage4_scalar_text(isolate)
  if (!nzchar(isolate) || !nrow(profile$profiles) || !"Isolate" %in% names(profile$profiles)) {
    return(stage4_empty_profile_summary())
  }
  key <- tolower(trimws(isolate))
  profile$profiles[tolower(trimws(as.character(profile$profiles$Isolate))) == key, , drop = FALSE]
}

stage4_isolate_evidence <- function(profile, isolate) {
  profile <- stage4_ensure_profile(profile)
  isolate <- stage4_scalar_text(isolate)
  if (!nzchar(isolate) || !nrow(profile$evidence) || !"Isolate" %in% names(profile$evidence)) {
    return(stage4_empty_evidence())
  }
  key <- tolower(trimws(isolate))
  profile$evidence[tolower(trimws(as.character(profile$evidence$Isolate))) == key, , drop = FALSE]
}

stage4_profile_overview <- function(profile) {
  profile <- stage4_ensure_profile(profile)
  statuses <- if (nrow(profile$profiles)) as.character(profile$profiles$Profile_Status) else character()
  loci <- if (nrow(profile$evidence)) trimws(as.character(profile$evidence$Locus)) else character()
  loci <- loci[nzchar(loci)]
  data.frame(
    Isolates = nrow(profile$profiles),
    Loci = length(unique(tolower(loci))),
    Concordant = sum(statuses %in% c("CONCORDANT_SPECIES", "CONCORDANT_GENUS")),
    Requires_Attention = sum(statuses %in% c(
      "GENUS_CONFLICT", "SPECIES_CONFLICT", "PARTIAL_EVIDENCE", "NO_TAXONOMY", "SINGLE_LOCUS"
    )),
    stringsAsFactors = FALSE
  )
}

stage4_profile_gate_error <- function(profile) {
  profile <- stage4_ensure_profile(profile)
  if (!nrow(profile$evidence)) return("Build a profile from completed single-locus PITAX projects.")
  if (length(unique(tolower(trimws(profile$evidence$Locus)))) < 2L) return("A multi-locus profile requires at least two distinct loci.")
  NULL
}

stage4_current_project_matches <- function(profile, project, source_name = "Current session") {
  profile <- stage4_ensure_profile(profile)
  stored <- profile$evidence[profile$evidence$Source == source_name, , drop = FALSE]
  if (!nrow(stored)) return(TRUE)
  current <- tryCatch(stage4_extract_project_evidence(project, source_name, ""), error = function(e) NULL)
  if (is.null(current) || nrow(current) != nrow(stored)) return(FALSE)
  compare_cols <- c(
    "Isolate", "Locus", "Analysis_Status", "Consensus_Revision", "Sequence",
    "Taxonomy_Status", "Recommended_Identification", "Recommended_Level", "RID"
  )
  key <- function(df) paste(tolower(trimws(df$Isolate)), tolower(trimws(df$Locus)), sep = "\r")
  stored <- stored[order(key(stored)), compare_cols, drop = FALSE]
  current <- current[order(key(current)), compare_cols, drop = FALSE]
  rownames(stored) <- NULL
  rownames(current) <- NULL
  identical(stored, current)
}

stage4_make_fasta <- function(profile) {
  profile <- stage4_ensure_profile(profile)
  if (!nrow(profile$evidence)) return("")
  lines <- character()
  for (i in seq_len(nrow(profile$evidence))) {
    ev <- profile$evidence[i, , drop = FALSE]
    if (!nzchar(ev$Sequence)) next
    header <- paste0(
      ">", clean_fasta_name(paste(ev$Isolate, ev$Locus, sep = "_")),
      " isolate=", clean_fasta_name(ev$Isolate),
      " locus=", clean_fasta_name(ev$Locus),
      " consensus_revision=", ev$Consensus_Revision,
      " source=", clean_fasta_name(ev$Source)
    )
    lines <- c(lines, header, wrap_sequence(ev$Sequence, 80))
  }
  paste(lines, collapse = "\n")
}

write_stage4_checkpoint_zip <- function(file, profile) {
  profile <- stage4_ensure_profile(profile)
  temp_dir <- tempfile("pitax_multilocus_")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  write.csv(profile$sources, file.path(temp_dir, "source_projects.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(profile$profiles, file.path(temp_dir, "isolate_profiles.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(profile$evidence, file.path(temp_dir, "per_locus_evidence.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(stage4_make_fasta(profile), file.path(temp_dir, "multi_locus_sequences.fasta"))
  writeLines(c(
    "PITAX Stage 4 multi-locus checkpoint",
    paste0("schema: ", profile$schema), paste0("algorithm: ", profile$algorithm),
    paste0("built_at: ", profile$built_at),
    "",
    "Interpretation contract:",
    "- Each locus remains separate evidence.",
    "- Profile status is concordance/conflict logic, not a flat vote.",
    "- Conflicts are retained and block a combined taxonomic call.",
    "- Taxon-specific marker recommendations are not yet automated in this foundation."
  ), file.path(temp_dir, "README.txt"))
  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  zip::zipr(file, files, root = temp_dir)
}

stage4_migrate_v4_state <- function(state) {
  if (!is.list(state)) state <- list()
  state$multilocus_profile <- stage4_empty_profile()
  old_log <- stage4_scalar_text(state$migration_log)
  note <- "Migrated project schema 4 to schema 5. Existing Stage 3, BLAST and taxonomy evidence was preserved; the Stage 4 multi-locus profile starts empty."
  state$migration_log <- if (nzchar(old_log)) paste(old_log, note, sep = " ") else note
  state
}
