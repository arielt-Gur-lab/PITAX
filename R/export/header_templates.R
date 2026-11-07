# Structured Custom FASTA Header profiles for Sequence Library compatibility.
# Authoritative state is a versioned profile model; text headers are derived.
# Never change source/biological identity or invent missing scientific metadata.

if (!exists("%||%", mode = "function", inherits = TRUE)) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
}

PITAX_HEADER_SCHEMA_VERSION <- 2L
PITAX_HEADER_BROWSER_STORAGE_KEY <- "PITAX_FASTA_HEADER_PROFILES_V2"

# ---------------------------------------------------------------------------
# Field registry (canonical Sequence Library tags)
# ---------------------------------------------------------------------------

pitax_header_field_registry <- function() {
  # Downstream destinations mirror Sequence Library fasta_header.py.
  list(
    final_name = list(
      field_id = "final_name", canonical_tag = NA_character_, display = "Final name",
      category = "primary", group = "A",
      description = "Stable PITAX record identifier (first FASTA token).",
      destination = "Primary FASTA identifier", destination_class = "primary",
      value_source = "system", editable = FALSE, data_type = "text",
      required_default = TRUE, aliases = character()
    ),
    organism = list(
      field_id = "organism", canonical_tag = "organism", display = "Organism",
      category = "primary", group = "A",
      description = "Organism / taxon text.",
      destination = "organism_text (first-class)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("org", "Organism")
    ),
    isolate = list(
      field_id = "isolate", canonical_tag = "isolate", display = "Isolate",
      category = "primary", group = "A",
      description = "Biomaterial isolate.",
      destination = "Isolate find-or-create (first-class)", destination_class = "first_class",
      value_source = "system", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("Isolate")
    ),
    sample = list(
      field_id = "sample", canonical_tag = "sample", display = "Sample",
      category = "primary", group = "A",
      description = "Biomaterial sample. Not inferred from Source_ID.",
      destination = "Sample find-or-create (first-class)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character()
    ),
    locus = list(
      field_id = "locus", canonical_tag = "locus", display = "Locus",
      category = "primary", group = "A",
      description = "Laboratory locus / gene.",
      destination = "Locus link (first-class)", destination_class = "first_class",
      value_source = "system", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("gene", "Locus")
    ),
    assay = list(
      field_id = "assay", canonical_tag = "assay", display = "Assay",
      category = "primary", group = "A",
      description = "Assay name. Requires locus in the same header.",
      destination = "Assay row when locus present (first-class)", destination_class = "first_class",
      value_source = "system", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character(), depends_on = "locus"
    ),
    role = list(
      field_id = "role", canonical_tag = "role", display = "Role",
      category = "primary", group = "A",
      description = "Experimental role. Downstream accepts only known vocabulary values.",
      destination = "Experimental role (first-class, vocabulary-gated)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character(),
      warning = "Sequence Library rejects unknown role values. Compatibility is not guaranteed without a shared vocabulary."
    ),
    molecule = list(
      field_id = "molecule", canonical_tag = "molecule", display = "Molecule",
      category = "primary", group = "A",
      description = "Molecule type: DNA, RNA, protein, or unknown.",
      destination = "molecule (first-class)", destination_class = "first_class",
      value_source = "system", editable = TRUE, data_type = "enum",
      enum_values = c("DNA", "RNA", "protein", "unknown"),
      required_default = FALSE, aliases = c("mol", "molecule_type")
    ),
    accession = list(
      field_id = "accession", canonical_tag = "accession", display = "Accession",
      category = "primary", group = "A",
      description = "External accession. Never fabricated by PITAX.",
      destination = "Record accession (first-class)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("acc")
    ),
    direction = list(
      field_id = "direction", canonical_tag = "direction", display = "Direction",
      category = "primary", group = "A",
      description = "Read direction annotation.",
      destination = "Description annotation (first-class)", destination_class = "first_class",
      value_source = "system", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("dir", "Direction")
    ),
    notes = list(
      field_id = "notes", canonical_tag = "notes", display = "Notes",
      category = "primary", group = "A",
      description = "Free-text notes.",
      destination = "Notes (first-class)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("note")
    ),
    pitax_id = list(
      field_id = "pitax_id", canonical_tag = "pitax_id", display = "PITAX ID",
      category = "primary", group = "A",
      description = "Stable PITAX record identifier.",
      destination = "Description annotation (first-class)", destination_class = "first_class",
      value_source = "system", editable = FALSE, data_type = "text",
      required_default = FALSE, aliases = c("pitax_record_id", "Record_ID")
    ),
    name = list(
      field_id = "name", canonical_tag = "name", display = "Name override",
      category = "primary", group = "A",
      description = "Advanced: overrides the first FASTA token. Can collide with other records.",
      destination = "Overrides record name (first-class)", destination_class = "first_class",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character(), advanced = TRUE
    ),
    experiment = list(
      field_id = "experiment", canonical_tag = "experiment", display = "Experiment",
      category = "experimental", group = "B",
      description = "Experiment label. Stored downstream in Notes.",
      destination = "Notes as experiment=... (not a dedicated column)", destination_class = "notes_folded",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("exp", "Experiment")
    ),
    seq_date = list(
      field_id = "seq_date", canonical_tag = "seq_date", display = "Sequencing date",
      category = "experimental", group = "B",
      description = "User-supplied sequencing date. Never inferred from export time.",
      destination = "Notes as seq_date=... (not a dedicated column)", destination_class = "notes_folded",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("date", "sequencing_date", "Seq_Date")
    ),
    treatment = list(
      field_id = "treatment", canonical_tag = "treatment", display = "Treatment",
      category = "experimental", group = "B",
      description = "Treatment annotation. Stored downstream in Notes.",
      destination = "Notes as treatment=...", destination_class = "notes_folded",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character()
    ),
    fungicide = list(
      field_id = "fungicide", canonical_tag = "fungicide", display = "Fungicide",
      category = "experimental", group = "B",
      description = "Fungicide annotation. Stored downstream in Notes.",
      destination = "Notes as fungicide=...", destination_class = "notes_folded",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character()
    ),
    location = list(
      field_id = "location", canonical_tag = "location", display = "Location",
      category = "experimental", group = "B",
      description = "Location annotation. Stored downstream in Notes.",
      destination = "Notes as location=...", destination_class = "notes_folded",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = c("Location")
    )
  )
}

pitax_header_builtin_fields <- function() {
  # Legacy token names still accepted when migrating old templates.
  c("Final_Name", "Isolate", "Locus", "Direction", "Source_ID", "Record_ID")
}

pitax_header_alias_map <- function() {
  reg <- pitax_header_field_registry()
  out <- list()
  for (f in reg) {
    keys <- unique(tolower(c(f$field_id, f$canonical_tag, f$aliases)))
    keys <- keys[nzchar(keys) & !is.na(keys)]
    for (k in keys) out[[k]] <- f$field_id
  }
  # Legacy built-ins that are not export tags.
  out[["final_name"]] <- "final_name"
  out[["source_id"]] <- "source_id"
  out
}

pitax_header_resolve_field_id <- function(token) {
  token <- trimws(as.character(token)[1])
  if (is.na(token) || !nzchar(token)) return(NA_character_)
  map <- pitax_header_alias_map()
  hit <- map[[tolower(token)]]
  if (!is.null(hit)) return(hit)
  if (grepl("^[A-Za-z][A-Za-z0-9_]*$", token)) return(tolower(token))
  NA_character_
}

pitax_header_field_def <- function(field_id, config = NULL) {
  field_id <- as.character(field_id)[1]
  reg <- pitax_header_field_registry()
  if (!is.null(reg[[field_id]])) return(reg[[field_id]])
  custom <- NULL
  if (is.list(config) && is.list(config$custom_fields)) custom <- config$custom_fields[[field_id]]
  if (!is.null(custom)) return(custom)
  # Unknown custom field synthesized on demand.
  list(
    field_id = field_id, canonical_tag = field_id,
    display = field_id, category = "custom", group = "C",
    description = "Custom / Extras field.",
    destination = "Extras + Notes (custom)", destination_class = "extras",
    value_source = "user", editable = TRUE, data_type = "text",
    required_default = FALSE, aliases = character()
  )
}

pitax_header_library_table <- function(config = NULL) {
  reg <- pitax_header_field_registry()
  rows <- lapply(reg, function(f) {
    data.frame(
      field_id = f$field_id,
      display = f$display,
      canonical_tag = if (is.na(f$canonical_tag)) "(primary token)" else f$canonical_tag,
      group = f$group,
      category = f$category,
      destination = f$destination,
      destination_class = f$destination_class,
      value_source = f$value_source,
      description = f$description,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.list(config) && is.list(config$custom_fields) && length(config$custom_fields)) {
    extra <- lapply(config$custom_fields, function(f) {
      data.frame(
        field_id = f$field_id, display = f$display %||% f$field_id,
        canonical_tag = f$canonical_tag %||% f$field_id, group = "C",
        category = "custom", destination = f$destination %||% "Extras + Notes (custom)",
        destination_class = "extras", value_source = "user",
        description = f$description %||% "Custom / Extras field.",
        stringsAsFactors = FALSE
      )
    })
    out <- rbind(out, do.call(rbind, extra))
  }
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Profile / config constructors
# ---------------------------------------------------------------------------

pitax_header_new_profile <- function(name = "Sequence Library default",
                                     ordered_fields = c("isolate", "locus", "molecule", "direction"),
                                     tag_style = "bare",
                                     primary_identifier = "final_name",
                                     required = character(),
                                     id = NULL) {
  name <- trimws(as.character(name)[1])
  if (is.na(name) || !nzchar(name) || nchar(name) > 80L) {
    stop("Profile name must be 1-80 characters.", call. = FALSE)
  }
  tag_style <- match.arg(as.character(tag_style)[1], c("bare", "bracket"))
  primary_identifier <- as.character(primary_identifier)[1]
  if (!identical(primary_identifier, "final_name") && !identical(primary_identifier, "name")) {
    stop("primary_identifier must be final_name or name.", call. = FALSE)
  }
  ordered_fields <- unique(as.character(ordered_fields))
  ordered_fields <- ordered_fields[nzchar(ordered_fields) & ordered_fields != "final_name"]
  if (anyDuplicated(tolower(ordered_fields))) stop("Duplicate fields in profile.", call. = FALSE)
  # Assay requires locus in the profile structure.
  if ("assay" %in% ordered_fields && !"locus" %in% ordered_fields) {
    ordered_fields <- c(ordered_fields[seq_len(match("assay", ordered_fields) - 1L)],
                        "locus",
                        ordered_fields[match("assay", ordered_fields):length(ordered_fields)])
    ordered_fields <- unique(ordered_fields)
  }
  required <- unique(as.character(required))
  required <- intersect(required, ordered_fields)
  if (is.null(id) || !nzchar(as.character(id)[1])) {
    id <- paste0("prof_", gsub("[^A-Za-z0-9]", "", substr(digest_like(name), 1, 12)))
  }
  list(
    id = as.character(id)[1],
    name = name,
    schema_version = PITAX_HEADER_SCHEMA_VERSION,
    primary_identifier = primary_identifier,
    tag_style = tag_style,
    ordered_fields = ordered_fields,
    required = required
  )
}

digest_like <- function(x) {
  # Lightweight non-crypto token for ids (no digest package dependency).
  bytes <- charToRaw(paste0(x, Sys.time(), sample.int(1e9, 1)))
  paste(as.character(as.integer(bytes) %% 16), collapse = "")
}

pitax_header_empty_values <- function(field_ids = character()) {
  cols <- c("Record_Key", as.character(field_ids))
  as.data.frame(setNames(rep(list(character()), length(cols)), cols),
                stringsAsFactors = FALSE)
}

pitax_header_default_profile <- function() {
  pitax_header_new_profile(
    name = "Sequence Library default",
    ordered_fields = c("isolate", "locus", "molecule", "direction"),
    tag_style = "bare",
    primary_identifier = "final_name",
    id = "prof_default"
  )
}

pitax_header_default_config <- function() {
  draft <- pitax_header_default_profile()
  list(
    header_schema = PITAX_HEADER_SCHEMA_VERSION,
    profiles = list(draft),
    active_profile_id = draft$id,
    draft = draft,
    custom_fields = list(),
    values = pitax_header_empty_values(character()),
    dirty = FALSE,
    persistence = list(
      browser_saved = FALSE,
      project_included = TRUE,
      portable_exported = FALSE,
      last_browser_save_at = NA_character_,
      last_status = "Unsaved draft"
    ),
    legacy_backup = NULL,
    migration_warnings = character(),
    # Derived legacy-compatible views for RUN_INFO / older helpers.
    fields = character(),
    templates = list(),
    active_template = draft$name,
    draft_template = pitax_header_profile_to_template(draft)
  )
}

# ---------------------------------------------------------------------------
# Template derivation / serialization
# ---------------------------------------------------------------------------

pitax_header_quote_value <- function(value) {
  value <- as.character(value)[1]
  if (is.na(value)) value <- ""
  if (grepl("[\r\n>]", value)) {
    stop("Header values must not contain line breaks or FASTA '>' symbols.", call. = FALSE)
  }
  # Quote when spaces or special characters would break bare key=value parsing.
  needs_quote <- grepl("[[:space:]=\\[\\]\"']", value, perl = TRUE) || !nzchar(value)
  if (needs_quote) {
    if (grepl('"', value, fixed = TRUE) && grepl("'", value, fixed = TRUE)) {
      stop("Value cannot contain both single and double quotes safely.", call. = FALSE)
    }
    if (grepl('"', value, fixed = TRUE)) return(paste0("'", value, "'"))
    return(paste0('"', value, '"'))
  }
  value
}

pitax_header_format_tag <- function(tag, value, style = "bare") {
  tag <- as.character(tag)[1]
  value <- pitax_header_quote_value(value)
  if (identical(style, "bracket")) paste0("[", tag, "=", value, "]") else paste0(tag, "=", value)
}

pitax_header_profile_to_template <- function(profile, config = NULL) {
  profile <- pitax_header_coerce_profile(profile)
  parts <- character()
  if (identical(profile$primary_identifier, "name")) {
    parts <- c(parts, "{name}")
  } else {
    parts <- c(parts, "{final_name}")
  }
  for (fid in profile$ordered_fields) {
    def <- pitax_header_field_def(fid, config)
    tag <- def$canonical_tag
    if (is.null(tag) || is.na(tag) || !nzchar(tag)) next
    # Structural preview uses {field_id} tokens for migration / advanced view.
    token <- paste0("{", fid, "}")
    if (identical(profile$tag_style, "bracket")) {
      parts <- c(parts, paste0("[", tag, "=", token, "]"))
    } else {
      # Preview form: tag={field}. Actual export quotes dynamically.
      parts <- c(parts, paste0(tag, "=", token))
    }
  }
  paste(parts, collapse = " ")
}

pitax_header_coerce_profile <- function(profile) {
  if (!is.list(profile)) stop("Invalid header profile.", call. = FALSE)
  name <- trimws(as.character(profile$name %||% "Untitled")[1])
  if (!nzchar(name)) name <- "Untitled"
  id <- as.character(profile$id %||% paste0("prof_", gsub("[^A-Za-z0-9]", "", name)))[1]
  tag_style <- as.character(profile$tag_style %||% "bare")[1]
  if (!tag_style %in% c("bare", "bracket")) tag_style <- "bare"
  primary <- as.character(profile$primary_identifier %||% "final_name")[1]
  if (!primary %in% c("final_name", "name")) primary <- "final_name"
  ordered <- unique(as.character(profile$ordered_fields %||% character()))
  ordered <- ordered[nzchar(ordered) & !tolower(ordered) %in% c("final_name")]
  # Normalize known aliases to field ids.
  ordered <- vapply(ordered, function(x) {
    id2 <- pitax_header_resolve_field_id(x)
    if (is.na(id2)) tolower(x) else id2
  }, character(1), USE.NAMES = FALSE)
  ordered <- unique(ordered)
  if ("assay" %in% ordered && !"locus" %in% ordered) {
    stop("Profile includes assay but not locus. Add locus before assay.", call. = FALSE)
  }
  if (anyDuplicated(ordered)) stop("Duplicate canonical fields in profile.", call. = FALSE)
  required <- unique(as.character(profile$required %||% character()))
  required <- intersect(vapply(required, function(x) {
    id2 <- pitax_header_resolve_field_id(x)
    if (is.na(id2)) tolower(x) else id2
  }, character(1), USE.NAMES = FALSE), ordered)
  list(
    id = id,
    name = name,
    schema_version = as.integer(profile$schema_version %||% PITAX_HEADER_SCHEMA_VERSION)[1],
    primary_identifier = primary,
    tag_style = tag_style,
    ordered_fields = ordered,
    required = required
  )
}

# ---------------------------------------------------------------------------
# Profile operations
# ---------------------------------------------------------------------------

pitax_header_find_profile_index <- function(config, id_or_name) {
  key <- as.character(id_or_name)[1]
  ids <- vapply(config$profiles, function(p) p$id, character(1))
  names <- vapply(config$profiles, function(p) p$name, character(1))
  idx <- match(key, ids)
  if (is.na(idx)) idx <- match(key, names)
  idx
}

pitax_header_get_active_profile <- function(config) {
  config <- pitax_header_coerce_config(config)
  config$draft
}

pitax_header_profile_add_field <- function(config, field_id, after = NULL) {
  config <- pitax_header_coerce_config(config)
  field_id <- pitax_header_resolve_field_id(field_id)
  if (is.na(field_id) || identical(field_id, "final_name")) {
    stop("Select a valid metadata field.", call. = FALSE)
  }
  if (identical(field_id, "source_id")) {
    stop("Source_ID is not an export tag. Use sample or pitax_id.", call. = FALSE)
  }
  draft <- config$draft
  if (field_id %in% draft$ordered_fields) {
    return(config) # idempotent
  }
  # Known tags cannot be duplicated via aliases.
  tags <- vapply(draft$ordered_fields, function(fid) {
    as.character(pitax_header_field_def(fid, config)$canonical_tag %||% fid)
  }, character(1))
  new_tag <- as.character(pitax_header_field_def(field_id, config)$canonical_tag %||% field_id)
  if (tolower(new_tag) %in% tolower(tags)) {
    stop(paste0("Canonical tag already present: ", new_tag), call. = FALSE)
  }
  if (identical(field_id, "assay") && !"locus" %in% draft$ordered_fields) {
    draft$ordered_fields <- c(draft$ordered_fields, "locus")
  }
  if (is.null(after) || !nzchar(as.character(after)[1]) ||
      !as.character(after)[1] %in% draft$ordered_fields) {
    draft$ordered_fields <- c(draft$ordered_fields, field_id)
  } else {
    pos <- match(as.character(after)[1], draft$ordered_fields)
    draft$ordered_fields <- append(draft$ordered_fields, field_id, after = pos)
  }
  config$draft <- draft
  config$dirty <- TRUE
  pitax_header_refresh_derived(config)
}

pitax_header_profile_remove_field <- function(config, field_id) {
  config <- pitax_header_coerce_config(config)
  field_id <- pitax_header_resolve_field_id(field_id)
  if (is.na(field_id)) stop("Select an existing field.", call. = FALSE)
  draft <- config$draft
  if (!field_id %in% draft$ordered_fields) return(config)
  if (identical(field_id, "locus") && "assay" %in% draft$ordered_fields) {
    stop("Remove assay before removing locus.", call. = FALSE)
  }
  draft$ordered_fields <- setdiff(draft$ordered_fields, field_id)
  draft$required <- setdiff(draft$required, field_id)
  config$draft <- draft
  config$dirty <- TRUE
  pitax_header_refresh_derived(config)
}

pitax_header_profile_reorder_fields <- function(config, ordered_fields) {
  config <- pitax_header_coerce_config(config)
  ordered <- unique(vapply(as.character(ordered_fields), function(x) {
    id <- pitax_header_resolve_field_id(x)
    if (is.na(id)) stop(paste("Unknown field:", x), call. = FALSE)
    id
  }, character(1), USE.NAMES = FALSE))
  ordered <- ordered[ordered != "final_name"]
  current <- config$draft$ordered_fields
  if (!setequal(ordered, current)) {
    stop("Reorder must include exactly the current profile fields.", call. = FALSE)
  }
  if ("assay" %in% ordered && !"locus" %in% ordered) {
    stop("Assay requires locus in the profile.", call. = FALSE)
  }
  config$draft$ordered_fields <- ordered
  config$dirty <- TRUE
  pitax_header_refresh_derived(config)
}

pitax_header_save_profile <- function(config, name = NULL) {
  config <- pitax_header_coerce_config(config)
  if (!is.null(name) && nzchar(trimws(as.character(name)[1]))) {
    config$draft$name <- trimws(as.character(name)[1])
  }
  draft <- pitax_header_coerce_profile(config$draft)
  idx <- pitax_header_find_profile_index(config, draft$id)
  if (is.na(idx)) {
    # Prefer update by name when id is new but name exists.
    idx <- pitax_header_find_profile_index(config, draft$name)
  }
  if (is.na(idx)) {
    config$profiles <- c(config$profiles, list(draft))
  } else {
    # Keep existing id when saving over a named profile.
    draft$id <- config$profiles[[idx]]$id
    config$profiles[[idx]] <- draft
  }
  config$draft <- draft
  config$active_profile_id <- draft$id
  config$dirty <- FALSE
  config$persistence$project_included <- TRUE
  config$persistence$last_status <- "Saved in project state (save .sangerproject to persist on disk)"
  pitax_header_refresh_derived(config)
}

pitax_header_load_profile <- function(config, id_or_name) {
  config <- pitax_header_coerce_config(config)
  idx <- pitax_header_find_profile_index(config, id_or_name)
  if (is.na(idx)) stop("Unknown profile.", call. = FALSE)
  config$draft <- pitax_header_coerce_profile(config$profiles[[idx]])
  config$active_profile_id <- config$draft$id
  config$dirty <- FALSE
  pitax_header_refresh_derived(config)
}

pitax_header_duplicate_profile <- function(config, id_or_name = NULL, new_name = NULL) {
  config <- pitax_header_coerce_config(config)
  src <- if (is.null(id_or_name)) config$draft else {
    idx <- pitax_header_find_profile_index(config, id_or_name)
    if (is.na(idx)) stop("Unknown profile.", call. = FALSE)
    config$profiles[[idx]]
  }
  src <- pitax_header_coerce_profile(src)
  src$id <- paste0("prof_", digest_like(src$name))
  src$name <- if (!is.null(new_name) && nzchar(trimws(as.character(new_name)[1]))) {
    trimws(as.character(new_name)[1])
  } else paste0(src$name, " copy")
  if (!is.na(pitax_header_find_profile_index(config, src$name))) {
    stop(paste0("A profile named '", src$name, "' already exists."), call. = FALSE)
  }
  config$profiles <- c(config$profiles, list(src))
  config$draft <- src
  config$active_profile_id <- src$id
  config$dirty <- FALSE
  pitax_header_refresh_derived(config)
}

pitax_header_rename_profile <- function(config, id_or_name, new_name) {
  config <- pitax_header_coerce_config(config)
  new_name <- trimws(as.character(new_name)[1])
  if (is.na(new_name) || !nzchar(new_name) || nchar(new_name) > 80L) {
    stop("Profile name must be 1-80 characters.", call. = FALSE)
  }
  idx <- pitax_header_find_profile_index(config, id_or_name)
  if (is.na(idx)) stop("Unknown profile.", call. = FALSE)
  other <- pitax_header_find_profile_index(config, new_name)
  if (!is.na(other) && other != idx) stop("Another profile already uses that name.", call. = FALSE)
  config$profiles[[idx]]$name <- new_name
  if (identical(config$active_profile_id, config$profiles[[idx]]$id) ||
      identical(config$draft$id, config$profiles[[idx]]$id)) {
    config$draft$name <- new_name
  }
  config$dirty <- FALSE
  pitax_header_refresh_derived(config)
}

pitax_header_delete_profile <- function(config, id_or_name) {
  config <- pitax_header_coerce_config(config)
  if (length(config$profiles) <= 1L) stop("Cannot delete the last profile.", call. = FALSE)
  idx <- pitax_header_find_profile_index(config, id_or_name)
  if (is.na(idx)) stop("Unknown profile.", call. = FALSE)
  deleted_id <- config$profiles[[idx]]$id
  config$profiles <- config$profiles[-idx]
  if (identical(config$active_profile_id, deleted_id)) {
    config$draft <- pitax_header_coerce_profile(config$profiles[[1]])
    config$active_profile_id <- config$draft$id
  }
  config$dirty <- FALSE
  pitax_header_refresh_derived(config)
}

pitax_header_create_custom_field <- function(config, label, tag = NULL, activate = TRUE) {
  config <- pitax_header_coerce_config(config)
  label <- trimws(as.character(label)[1])
  if (is.na(label) || !nzchar(label)) stop("Enter a custom field name.", call. = FALSE)
  tag <- if (is.null(tag) || !nzchar(trimws(as.character(tag)[1]))) {
    gsub("[^A-Za-z0-9_]+", "_", label)
  } else trimws(as.character(tag)[1])
  tag <- tolower(tag)
  if (!grepl("^[a-z][a-z0-9_]*$", tag)) {
    stop("Custom tags must start with a letter and contain only letters, digits or underscores.", call. = FALSE)
  }
  # Reject collisions with canonical tags / aliases.
  resolved <- pitax_header_alias_map()[[tag]]
  if (!is.null(resolved) && resolved %in% names(pitax_header_field_registry())) {
    stop(paste0("Tag '", tag, "' is reserved for a Sequence Library field."), call. = FALSE)
  }
  if (!is.null(config$custom_fields[[tag]])) {
    field_id <- tag
    created <- FALSE
  } else {
    if (length(config$custom_fields) >= 40L) stop("Maximum 40 custom fields per project.", call. = FALSE)
    config$custom_fields[[tag]] <- list(
      field_id = tag, canonical_tag = tag, display = label,
      category = "custom", group = "C",
      description = "Custom / Extras field.",
      destination = "Extras + Notes (custom)", destination_class = "extras",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = character()
    )
    if (!tag %in% names(config$values)) config$values[[tag]] <- rep("", nrow(config$values))
    field_id <- tag
    created <- TRUE
  }
  if (isTRUE(activate)) config <- pitax_header_profile_add_field(config, field_id)
  list(config = pitax_header_refresh_derived(config), field = field_id, created = created)
}

pitax_header_delete_custom_field <- function(config, field_id) {
  config <- pitax_header_coerce_config(config)
  field_id <- pitax_header_resolve_field_id(field_id)
  if (is.na(field_id) || is.null(config$custom_fields[[field_id]])) {
    stop("Select an existing custom field.", call. = FALSE)
  }
  used_by <- character()
  for (p in config$profiles) {
    if (field_id %in% p$ordered_fields) used_by <- c(used_by, p$name)
  }
  if (field_id %in% config$draft$ordered_fields) used_by <- unique(c(used_by, config$draft$name))
  if (length(used_by)) {
    stop(paste0("Cannot delete custom field: still used by profile(s): ",
                paste(unique(used_by), collapse = ", "),
                ". Remove it from those profiles first."), call. = FALSE)
  }
  config$custom_fields[[field_id]] <- NULL
  if (field_id %in% names(config$values)) config$values[[field_id]] <- NULL
  pitax_header_refresh_derived(config)
}

# ---------------------------------------------------------------------------
# Derived views + coercion / migration
# ---------------------------------------------------------------------------

pitax_header_refresh_derived <- function(config) {
  profiles <- config$profiles
  names_vec <- vapply(profiles, function(p) p$name, character(1))
  templates <- setNames(lapply(profiles, function(p) {
    pitax_header_profile_to_template(p, config)
  }), names_vec)
  draft <- config$draft
  # Editable metadata columns for the active draft (exclude primary-only final_name).
  editable <- character()
  for (fid in draft$ordered_fields) {
    def <- pitax_header_field_def(fid, config)
    if (isTRUE(def$editable) || identical(def$value_source, "user") ||
        identical(def$value_source, "system")) {
      # System fields still appear so users can review/override where allowed.
      if (!identical(fid, "final_name")) editable <- c(editable, fid)
    }
  }
  # Ensure value columns exist for all known editable/custom fields used anywhere.
  all_fields <- unique(c(
    unlist(lapply(profiles, function(p) p$ordered_fields), use.names = FALSE),
    draft$ordered_fields,
    names(config$custom_fields)
  ))
  all_fields <- setdiff(all_fields, c("final_name"))
  if (!is.data.frame(config$values)) config$values <- pitax_header_empty_values(all_fields)
  if (!"Record_Key" %in% names(config$values)) {
    config$values <- cbind(Record_Key = character(nrow(config$values)), config$values,
                           stringsAsFactors = FALSE)
  }
  for (fid in all_fields) {
    if (!fid %in% names(config$values)) config$values[[fid]] <- rep("", nrow(config$values))
  }
  # Keep Record_Key + known fields; do not drop historical columns silently.
  keep <- unique(c("Record_Key", all_fields, setdiff(names(config$values), "Record_Key")))
  keep <- keep[keep %in% names(config$values)]
  config$values <- config$values[, keep, drop = FALSE]
  for (col in names(config$values)) {
    config$values[[col]] <- as.character(config$values[[col]])
    config$values[[col]][is.na(config$values[[col]])] <- ""
  }
  config$values <- config$values[!duplicated(config$values$Record_Key) &
                                   nzchar(config$values$Record_Key), , drop = FALSE]
  rownames(config$values) <- NULL

  config$fields <- editable
  config$templates <- templates
  config$active_template <- draft$name
  config$draft_template <- pitax_header_profile_to_template(draft, config)
  config$header_schema <- PITAX_HEADER_SCHEMA_VERSION
  config
}

pitax_header_is_legacy_config <- function(config) {
  if (!is.list(config)) return(FALSE)
  if (!is.null(config$header_schema) && as.integer(config$header_schema)[1] >= 2L &&
      is.list(config$draft) && !is.null(config$draft$ordered_fields)) {
    return(FALSE)
  }
  # Legacy shape: fields + templates + draft_template strings.
  is.character(config$fields) || is.character(config$draft_template) ||
    (is.list(config$templates) && length(config$templates) && is.character(config$templates[[1]]))
}

pitax_header_parse_legacy_template <- function(template) {
  template <- as.character(template)[1]
  if (is.na(template)) template <- ""
  template <- sub("^[[:space:]]*>[[:space:]]*", "", template)
  tokens <- regmatches(template, gregexpr("\\{[A-Za-z][A-Za-z0-9_]*\\}", template, perl = TRUE))[[1]]
  if (!length(tokens)) return(list(primary = "final_name", fields = character(), warnings = character()))
  names <- substring(tokens, 2, nchar(tokens) - 1L)
  warnings <- character()
  fields <- character()
  primary <- "final_name"
  for (tok in names) {
    fid <- pitax_header_resolve_field_id(tok)
    if (is.na(fid)) {
      warnings <- c(warnings, paste0("Unrecognized token {", tok, "} kept as custom field."))
      fid <- tolower(tok)
    }
    if (fid %in% c("final_name", "source_id")) {
      if (identical(fid, "final_name")) primary <- "final_name"
      next
    }
    if (identical(fid, "name") && identical(tok, names[1])) {
      primary <- "name"
      next
    }
    fields <- c(fields, fid)
  }
  list(primary = primary, fields = unique(fields), warnings = warnings,
       tag_style = if (grepl("\\[[A-Za-z]", template)) "bracket" else "bare")
}

pitax_header_migrate_legacy_config <- function(config) {
  warnings <- character()
  legacy_backup <- config
  default <- pitax_header_default_config()
  if (!is.list(config)) {
    default$migration_warnings <- "Missing header config; installed Sequence Library default."
    default$legacy_backup <- legacy_backup
    return(default)
  }
  old_fields <- as.character(config$fields %||% character())
  old_values <- config$values
  templates <- config$templates
  if (!is.list(templates) || !length(templates)) {
    templates <- list("Research metadata" = as.character(config$draft_template %||%
                                                           pitax_header_default_template())[1])
  }
  if (is.null(names(templates)) || any(!nzchar(names(templates)))) {
    names(templates) <- paste0("Profile ", seq_along(templates))
  }
  profiles <- list()
  for (nm in names(templates)) {
    parsed <- pitax_header_parse_legacy_template(templates[[nm]])
    warnings <- c(warnings, parsed$warnings)
    # Custom fields referenced only as tokens.
    for (fid in parsed$fields) {
      if (!fid %in% names(pitax_header_field_registry()) &&
          !fid %in% c("final_name", "source_id")) {
        # will register below
      }
    }
    profiles[[length(profiles) + 1L]] <- tryCatch(
      pitax_header_new_profile(
        name = nm,
        ordered_fields = parsed$fields,
        tag_style = parsed$tag_style,
        primary_identifier = parsed$primary,
        id = paste0("legacy_", gsub("[^A-Za-z0-9]+", "_", nm))
      ),
      error = function(e) {
        warnings <<- c(warnings, paste0("Profile '", nm, "' needs review: ", conditionMessage(e)))
        pitax_header_new_profile(name = nm, ordered_fields = character(),
                                 id = paste0("legacy_", gsub("[^A-Za-z0-9]+", "_", nm)))
      }
    )
  }
  if (!length(profiles)) profiles <- list(pitax_header_default_profile())

  # Migrate custom field definitions from legacy field list.
  custom_fields <- list()
  for (f in old_fields) {
    fid <- pitax_header_resolve_field_id(f)
    if (is.na(fid)) fid <- tolower(f)
    if (fid %in% names(pitax_header_field_registry())) next
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", f)) {
      warnings <- c(warnings, paste0("Dropped invalid legacy field name: ", f))
      next
    }
    custom_fields[[fid]] <- list(
      field_id = fid, canonical_tag = fid, display = f,
      category = "custom", group = "C",
      description = "Migrated custom field.",
      destination = "Extras + Notes (custom)", destination_class = "extras",
      value_source = "user", editable = TRUE, data_type = "text",
      required_default = FALSE, aliases = if (!identical(f, fid)) f else character()
    )
  }

  # Migrate values: map legacy column names to field ids; preserve Record_Key.
  values <- pitax_header_empty_values(character())
  if (is.data.frame(old_values) && "Record_Key" %in% names(old_values) && nrow(old_values)) {
    values <- data.frame(Record_Key = as.character(old_values$Record_Key), stringsAsFactors = FALSE)
    for (col in setdiff(names(old_values), "Record_Key")) {
      fid <- pitax_header_resolve_field_id(col)
      if (is.na(fid)) fid <- tolower(col)
      if (fid %in% c("final_name", "source_id")) next
      vals <- as.character(old_values[[col]])
      vals[is.na(vals)] <- ""
      if (fid %in% names(values)) {
        # Prefer non-empty when merging alias collisions.
        empty <- !nzchar(values[[fid]])
        values[[fid]][empty] <- vals[empty]
      } else {
        values[[fid]] <- vals
      }
    }
  }

  active_name <- as.character(config$active_template %||% names(templates)[1])[1]
  idx <- match(active_name, vapply(profiles, function(p) p$name, character(1)))
  if (is.na(idx)) idx <- 1L
  draft_template <- as.character(config$draft_template %||% templates[[active_name]] %||% "")[1]
  if (nzchar(draft_template)) {
    parsed_draft <- pitax_header_parse_legacy_template(draft_template)
    warnings <- c(warnings, parsed_draft$warnings)
    draft <- tryCatch(
      pitax_header_new_profile(
        name = profiles[[idx]]$name,
        ordered_fields = parsed_draft$fields,
        tag_style = parsed_draft$tag_style,
        primary_identifier = parsed_draft$primary,
        id = profiles[[idx]]$id
      ),
      error = function(e) {
        warnings <<- c(warnings, paste0("Draft migration issue: ", conditionMessage(e)))
        profiles[[idx]]
      }
    )
  } else {
    draft <- profiles[[idx]]
  }

  out <- list(
    header_schema = PITAX_HEADER_SCHEMA_VERSION,
    profiles = profiles,
    active_profile_id = draft$id,
    draft = draft,
    custom_fields = custom_fields,
    values = values,
    dirty = FALSE,
    persistence = list(
      browser_saved = FALSE,
      project_included = TRUE,
      portable_exported = FALSE,
      last_browser_save_at = NA_character_,
      last_status = "Migrated from legacy header templates"
    ),
    legacy_backup = legacy_backup,
    migration_warnings = unique(warnings),
    fields = character(),
    templates = list(),
    active_template = draft$name,
    draft_template = ""
  )
  pitax_header_refresh_derived(out)
}

pitax_header_coerce_config <- function(config) {
  if (is.null(config)) return(pitax_header_default_config())
  if (!is.list(config)) {
    out <- pitax_header_default_config()
    out$migration_warnings <- "Invalid header config type; defaults installed. Original value preserved in legacy_backup."
    out$legacy_backup <- config
    return(out)
  }
  if (pitax_header_is_legacy_config(config)) {
    return(pitax_header_migrate_legacy_config(config))
  }
  # Already schema 2(+): normalize carefully without discarding data.
  warnings <- as.character(config$migration_warnings %||% character())
  profiles <- config$profiles
  if (!is.list(profiles) || !length(profiles)) {
    warnings <- c(warnings, "No profiles found; installed default profile. Previous config kept in legacy_backup.")
    out <- pitax_header_default_config()
    out$legacy_backup <- config
    out$values <- if (is.data.frame(config$values)) config$values else out$values
    out$migration_warnings <- warnings
    return(out)
  }
  profiles <- lapply(profiles, function(p) {
    tryCatch(pitax_header_coerce_profile(p), error = function(e) {
      warnings <<- c(warnings, paste0("Profile '", p$name %||% "?", "' invalid: ", conditionMessage(e)))
      pitax_header_default_profile()
    })
  })
  # Deduplicate names by appending suffix rather than dropping.
  seen <- character()
  for (i in seq_along(profiles)) {
    nm <- profiles[[i]]$name
    if (nm %in% seen) {
      profiles[[i]]$name <- paste0(nm, " (", i, ")")
    }
    seen <- c(seen, profiles[[i]]$name)
  }
  draft <- config$draft
  draft <- tryCatch(pitax_header_coerce_profile(draft %||% profiles[[1]]), error = function(e) {
    warnings <<- c(warnings, paste0("Draft invalid: ", conditionMessage(e)))
    profiles[[1]]
  })
  custom_fields <- config$custom_fields
  if (!is.list(custom_fields)) custom_fields <- list()
  values <- config$values
  if (!is.data.frame(values) || !"Record_Key" %in% names(values)) {
    values <- pitax_header_empty_values(character())
  }
  out <- list(
    header_schema = PITAX_HEADER_SCHEMA_VERSION,
    profiles = profiles,
    active_profile_id = as.character(config$active_profile_id %||% draft$id)[1],
    draft = draft,
    custom_fields = custom_fields,
    values = values,
    dirty = isTRUE(config$dirty),
    persistence = if (is.list(config$persistence)) config$persistence else list(
      browser_saved = FALSE, project_included = TRUE, portable_exported = FALSE,
      last_browser_save_at = NA_character_, last_status = "Loaded"
    ),
    legacy_backup = config$legacy_backup,
    migration_warnings = unique(warnings),
    fields = character(),
    templates = list(),
    active_template = draft$name,
    draft_template = ""
  )
  pitax_header_refresh_derived(out)
}

# ---------------------------------------------------------------------------
# Portable JSON / browser payload (profiles only — never record values)
# ---------------------------------------------------------------------------

pitax_header_profiles_portable <- function(config) {
  config <- pitax_header_coerce_config(config)
  list(
    schema = PITAX_HEADER_SCHEMA_VERSION,
    storage_key = PITAX_HEADER_BROWSER_STORAGE_KEY,
    active_profile_id = config$active_profile_id,
    profiles = lapply(config$profiles, function(p) {
      list(
        id = p$id, name = p$name, schema_version = p$schema_version,
        primary_identifier = p$primary_identifier, tag_style = p$tag_style,
        ordered_fields = I(p$ordered_fields), required = I(p$required)
      )
    }),
    custom_fields = lapply(config$custom_fields, function(f) {
      list(field_id = f$field_id, canonical_tag = f$canonical_tag,
           display = f$display, category = f$category, group = f$group,
           description = f$description, destination = f$destination,
           destination_class = f$destination_class, value_source = f$value_source,
           editable = isTRUE(f$editable), data_type = f$data_type)
    })
  )
}

pitax_header_merge_portable_profiles <- function(config, portable,
                                                 precedence = c("project", "browser")) {
  precedence <- match.arg(precedence)
  config <- pitax_header_coerce_config(config)
  if (!is.list(portable) || is.null(portable$profiles)) return(config)
  incoming <- lapply(portable$profiles, function(p) {
    tryCatch(pitax_header_coerce_profile(p), error = function(e) NULL)
  })
  incoming <- Filter(Negate(is.null), incoming)
  if (is.list(portable$custom_fields)) {
    for (nm in names(portable$custom_fields)) {
      f <- portable$custom_fields[[nm]]
      if (is.null(config$custom_fields[[nm]])) config$custom_fields[[nm]] <- f
    }
  }
  for (p in incoming) {
    idx_id <- pitax_header_find_profile_index(config, p$id)
    idx_name <- pitax_header_find_profile_index(config, p$name)
    if (identical(precedence, "browser")) {
      # Fresh browser session: browser layouts replace matching ids and add new names.
      if (!is.na(idx_id)) {
        config$profiles[[idx_id]] <- p
        if (identical(config$active_profile_id, p$id) ||
            identical(config$draft$id, p$id)) {
          config$draft <- p
          config$active_profile_id <- p$id
        }
      } else if (is.na(idx_name)) {
        config$profiles <- c(config$profiles, list(p))
      } else {
        # Same name, different id: keep project/session copy; import renamed.
        p$name <- paste0(p$name, " (browser)")
        config$profiles <- c(config$profiles, list(p))
      }
    } else {
      # Project takes priority. Only add browser profiles that do not collide.
      if (!is.na(idx_id) || !is.na(idx_name)) next
      config$profiles <- c(config$profiles, list(p))
    }
  }
  pitax_header_refresh_derived(config)
}

pitax_header_import_portable_json <- function(config, json_text) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    # Minimal parser fallback: write temp and use jsonlite if available; else stop.
    stop("jsonlite is required to import portable header profiles.", call. = FALSE)
  }
  portable <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)
  pitax_header_merge_portable_profiles(config, portable, precedence = "project")
}

# ---------------------------------------------------------------------------
# Values / catalog / editor
# ---------------------------------------------------------------------------

pitax_header_record_key <- function(component_id, record_id) {
  component_id <- as.character(component_id)[1]
  namespace <- if (identical(component_id, "qc_trimmed_sequences")) "read" else
    if (component_id %in% c("consensus_analysis_sequences", "blast_query_sequences")) "analysis" else
    if (identical(component_id, "multilocus_sequences")) "multilocus" else
      stop(paste("Unknown sequence component:", component_id), call. = FALSE)
  paste0(namespace, "::", as.character(record_id)[1])
}

pitax_header_catalog <- function(snap) {
  rows <- list()
  push <- function(key, final_name, isolate = "", locus = "", direction = "",
                   source = "", record = "", assay = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      Record_Key = key, Final_Name = final_name, Isolate = isolate, Locus = locus,
      Direction = direction, Source_ID = source, Record_ID = record, Assay = assay,
      stringsAsFactors = FALSE
    )
  }
  asg <- snap$read_assignments
  if (is.data.frame(asg) && nrow(asg)) {
    for (i in seq_len(nrow(asg))) {
      sid <- as.character(asg$Source_ID[i])
      name <- if (exists("pitax_export_final_name", mode = "function")) {
        pitax_export_final_name(snap, sid)
      } else as.character(asg$Final_Name[i] %||% sid)
      assay <- if ("Assay_ID" %in% names(asg)) as.character(asg$Assay_ID[i]) else
        if ("Assay_Name" %in% names(asg)) as.character(asg$Assay_Name[i]) else ""
      push(paste0("read::", sid), name, as.character(asg$Isolate[i]),
           as.character(asg$Locus[i]), as.character(asg$Direction[i]), sid, sid, assay)
    }
  } else if (length(snap$results)) {
    for (sid in names(snap$results)) {
      nm <- if (exists("pitax_export_final_name", mode = "function")) {
        pitax_export_final_name(snap, sid)
      } else sid
      push(paste0("read::", sid), nm, source = sid, record = sid)
    }
  }
  cs <- snap$consensus_set
  if (is.list(cs) && is.list(cs$records)) {
    for (id in names(cs$records)) {
      rec <- cs$records[[id]]
      push(paste0("analysis::", id), as.character(rec$final_name %||% id),
           as.character(rec$isolate %||% ""), as.character(rec$locus %||% ""),
           "oriented", paste(as.character(rec$source_read_ids %||% ""), collapse = ","), id,
           as.character(rec$assay %||% rec$assay_name %||% ""))
    }
  }
  ev <- tryCatch(pitax_export_extract_multilocus_sequences(snap), error = function(e) list())
  if (length(ev)) for (id in names(ev)) {
    rec <- ev[[id]]
    push(paste0("multilocus::", id), as.character(rec$final_name %||% id),
         as.character(rec$isolate %||% ""), as.character(rec$locus %||% ""),
         "oriented", "", id, as.character(rec$assay %||% ""))
  }
  if (!length(rows)) return(data.frame(
    Record_Key = character(), Final_Name = character(), Isolate = character(),
    Locus = character(), Direction = character(), Source_ID = character(),
    Record_ID = character(), Assay = character(), stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$Record_Key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

pitax_header_system_value <- function(field_id, catalog_row, rec = NULL, record_id = NULL) {
  field_id <- as.character(field_id)[1]
  if (identical(field_id, "final_name")) {
    return(as.character(catalog_row$Final_Name %||% rec$final_name %||% record_id %||% "")[1])
  }
  if (identical(field_id, "isolate")) return(as.character(catalog_row$Isolate %||% rec$isolate %||% "")[1])
  if (identical(field_id, "locus")) return(as.character(catalog_row$Locus %||% rec$locus %||% "")[1])
  if (identical(field_id, "direction")) {
    raw <- as.character(catalog_row$Direction %||% rec$direction %||% "")[1]
    # Compact F/R when clearly Forward/Reverse; otherwise keep original.
    if (grepl("^forward", raw, ignore.case = TRUE)) return("F")
    if (grepl("^reverse", raw, ignore.case = TRUE)) return("R")
    return(raw)
  }
  if (identical(field_id, "assay")) return(as.character(catalog_row$Assay %||% rec$assay %||% "")[1])
  if (identical(field_id, "pitax_id")) return(as.character(catalog_row$Record_ID %||% record_id %||% "")[1])
  if (identical(field_id, "molecule")) return("DNA") # justified default for Sanger DNA exports
  ""
}

pitax_header_active_metadata_fields <- function(config) {
  config <- pitax_header_coerce_config(config)
  # User-visible editor columns for the active draft only.
  fields <- config$draft$ordered_fields
  # Always allow editing user/system-override fields present in profile.
  fields
}

pitax_header_active_edit_table <- function(snap, config) {
  config <- pitax_header_coerce_config(config)
  cat <- pitax_header_catalog(snap)
  active <- pitax_header_active_metadata_fields(config)
  for (fid in active) {
    idx <- match(cat$Record_Key, config$values$Record_Key)
    stored <- if (fid %in% names(config$values)) config$values[[fid]][idx] else rep(NA_character_, nrow(cat))
    stored[is.na(stored)] <- ""
    sys <- vapply(seq_len(nrow(cat)), function(i) {
      pitax_header_system_value(fid, cat[i, , drop = FALSE])
    }, character(1))
    # Prefer explicit stored values; fall back to system-derived.
    merged <- stored
    empty <- !nzchar(trimws(merged))
    merged[empty] <- sys[empty]
    # Column display uses field_id; UI maps labels.
    cat[[fid]] <- merged
  }
  columns <- c("Record_Key", "Final_Name", "Isolate", "Locus", "Direction", active)
  columns <- unique(columns[columns %in% names(cat)])
  cat[, columns, drop = FALSE]
}

pitax_header_set_value <- function(config, record_keys, field, value) {
  config <- pitax_header_coerce_config(config)
  field_id <- pitax_header_resolve_field_id(field)
  if (is.na(field_id)) stop("Unknown metadata field.", call. = FALSE)
  if (identical(field_id, "final_name") || identical(field_id, "pitax_id")) {
    stop("System identity fields cannot be overwritten in export metadata.", call. = FALSE)
  }
  def <- pitax_header_field_def(field_id, config)
  if (identical(def$data_type, "enum")) {
    allowed <- def$enum_values %||% character()
    value <- as.character(value)[1]
    if (!value %in% allowed && nzchar(value)) {
      stop(paste0("Invalid value for ", def$display, ". Allowed: ",
                  paste(allowed, collapse = ", ")), call. = FALSE)
    }
  } else {
    value <- as.character(value)[1]
  }
  if (is.na(value) || grepl("[\r\n>]", value)) {
    stop("Header values must not contain line breaks or FASTA '>' symbols.", call. = FALSE)
  }
  if (!field_id %in% names(config$values)) {
    config$values[[field_id]] <- rep("", nrow(config$values))
  }
  for (key in unique(as.character(record_keys))) {
    if (!nzchar(key)) next
    idx <- match(key, config$values$Record_Key)
    if (is.na(idx)) {
      row <- as.data.frame(as.list(setNames(rep("", length(names(config$values))),
                                            names(config$values))), stringsAsFactors = FALSE)
      row$Record_Key <- key
      config$values <- rbind(config$values, row)
      idx <- nrow(config$values)
    }
    config$values[[field_id]][idx] <- value
  }
  config$dirty <- TRUE
  pitax_header_refresh_derived(config)
}

pitax_header_apply_isolate_value <- function(config, snap, isolate, field, value, locus = NULL) {
  config <- pitax_header_coerce_config(config)
  cat <- pitax_header_catalog(snap)
  keys <- cat$Record_Key[cat$Isolate == as.character(isolate)[1]]
  if (!is.null(locus) && nzchar(as.character(locus)[1])) {
    keys <- cat$Record_Key[cat$Isolate == as.character(isolate)[1] &
                             cat$Locus == as.character(locus)[1]]
  }
  # When applying isolate-level metadata for a locus-specific field, require locus.
  field_id <- pitax_header_resolve_field_id(field)
  if (field_id %in% c("locus", "assay", "direction") &&
      (is.null(locus) || !nzchar(as.character(locus)[1]))) {
    stop("Locus-specific fields require an explicit locus filter so other loci are not overwritten.",
         call. = FALSE)
  }
  pitax_header_set_value(config, keys, field, value)
}

# ---------------------------------------------------------------------------
# Rendering / validation
# ---------------------------------------------------------------------------

pitax_header_resolve_record_values <- function(rec, component_id, record_id, config) {
  config <- pitax_header_coerce_config(config)
  key <- pitax_header_record_key(component_id, record_id)
  catalog_row <- data.frame(
    Record_Key = key,
    Final_Name = as.character(rec$final_name %||% record_id)[1],
    Isolate = as.character(rec$isolate %||% "")[1],
    Locus = as.character(rec$locus %||% "")[1],
    Direction = as.character(rec$direction %||% "")[1],
    Source_ID = if (startsWith(key, "read::")) as.character(record_id)[1] else "",
    Record_ID = as.character(record_id)[1],
    Assay = as.character(rec$assay %||% rec$assay_name %||% "")[1],
    stringsAsFactors = FALSE
  )
  idx <- match(key, config$values$Record_Key)
  out <- list()
  out$final_name <- pitax_header_system_value("final_name", catalog_row, rec, record_id)
  for (fid in unique(c(config$draft$ordered_fields, "name"))) {
    stored <- if (!is.na(idx) && fid %in% names(config$values)) config$values[[fid]][idx] else ""
    if (is.na(stored)) stored <- ""
    if (nzchar(trimws(stored))) {
      out[[fid]] <- stored
    } else {
      out[[fid]] <- pitax_header_system_value(fid, catalog_row, rec, record_id)
    }
  }
  out
}

pitax_header_validate_profile_structure <- function(config) {
  config <- pitax_header_coerce_config(config)
  draft <- config$draft
  issues <- character()
  if ("assay" %in% draft$ordered_fields && !"locus" %in% draft$ordered_fields) {
    issues <- c(issues, "Assay requires locus in the profile.")
  }
  if ("role" %in% draft$ordered_fields) {
    issues <- c(issues, paste(
      "Warning: role is vocabulary-gated downstream;",
      "PITAX cannot guarantee acceptance without a shared role vocabulary."
    ))
  }
  if (identical(draft$primary_identifier, "name")) {
    issues <- c(issues, "Warning: name override replaces the primary FASTA identifier.")
  }
  tags <- character()
  for (fid in draft$ordered_fields) {
    tag <- as.character(pitax_header_field_def(fid, config)$canonical_tag %||% fid)
    if (tolower(tag) %in% tolower(tags)) {
      issues <- c(issues, paste0("Duplicate canonical tag: ", tag))
    }
    tags <- c(tags, tag)
  }
  issues
}

pitax_header_render <- function(rec, component_id, record_id, config) {
  config <- pitax_header_coerce_config(config)
  draft <- config$draft
  vals <- pitax_header_resolve_record_values(rec, component_id, record_id, config)
  struct_issues <- pitax_header_validate_profile_structure(config)
  hard <- struct_issues[!startsWith(struct_issues, "Warning:")]
  if (length(hard)) stop(paste(hard, collapse = " "), call. = FALSE)

  if (identical(draft$primary_identifier, "name") && nzchar(vals$name %||% "")) {
    primary <- vals$name
  } else {
    primary <- vals$final_name
  }
  primary <- as.character(primary)[1]
  if (!nzchar(trimws(primary))) {
    stop(paste0("Missing primary identifier for ", record_id), call. = FALSE)
  }
  if (grepl("[[:space:]>]", primary)) {
    stop(paste0("Primary FASTA identifier for ", record_id,
                " must not contain spaces or '>'."), call. = FALSE)
  }

  parts <- primary
  for (fid in draft$ordered_fields) {
    def <- pitax_header_field_def(fid, config)
    tag <- def$canonical_tag
    if (is.null(tag) || is.na(tag) || !nzchar(tag)) next
    value <- as.character(vals[[fid]] %||% "")[1]
    if (is.na(value)) value <- ""
    value <- trimws(value)
    required <- fid %in% draft$required
    if (!nzchar(value)) {
      if (required) {
        stop(paste0("Missing required {", fid, "} for ",
                    as.character(rec$final_name %||% record_id)), call. = FALSE)
      }
      next # omit optional empty annotations entirely
    }
    if (identical(fid, "molecule")) {
      allowed <- c("DNA", "RNA", "protein", "unknown")
      if (!value %in% allowed) {
        stop(paste0("Invalid molecule for ", record_id, ": ", value), call. = FALSE)
      }
    }
    if (identical(fid, "assay") && !nzchar(as.character(vals$locus %||% ""))) {
      stop(paste0("Assay requires a populated locus for ",
                  as.character(rec$final_name %||% record_id)), call. = FALSE)
    }
    parts <- c(parts, pitax_header_format_tag(tag, value, draft$tag_style))
  }
  header <- paste(parts, collapse = " ")
  if (!nzchar(trimws(header))) stop("Empty FASTA header.", call. = FALSE)
  header
}

pitax_header_apply_to_records <- function(records, component_id, config) {
  if (!length(records)) return(records)
  ids <- character()
  for (record_id in names(records)) {
    rec <- records[[record_id]]
    if (isTRUE(rec$missing) || !nzchar(as.character(rec$sequence %||% ""))) next
    rendered <- pitax_header_render(rec, component_id, record_id, config)
    primary <- strsplit(trimws(rendered), "[[:space:]]+")[[1]][1]
    if (tolower(primary) %in% tolower(ids)) {
      stop(paste0("Duplicate first FASTA identifier: ", primary,
                  ". Keep unique Final names or avoid colliding name overrides."),
           call. = FALSE)
    }
    ids <- c(ids, primary)
    records[[record_id]]$custom_header <- rendered
  }
  records
}

# ---------------------------------------------------------------------------
# Legacy-compatible helpers used by older tests / UI remnants
# ---------------------------------------------------------------------------

pitax_header_default_template <- function() {
  pitax_header_profile_to_template(pitax_header_default_profile())
}

pitax_header_validate_field <- function(field, existing = character()) {
  field <- trimws(as.character(field)[1])
  if (is.na(field) || !grepl("^[A-Za-z][A-Za-z0-9_]*$", field)) {
    stop("Metadata field names must start with a letter and contain only letters, digits or underscores.",
         call. = FALSE)
  }
  fid <- pitax_header_resolve_field_id(field)
  if (!is.na(fid) && fid %in% names(pitax_header_field_registry()) &&
      !identical(pitax_header_field_def(fid)$destination_class, "extras")) {
    # Allow registering only truly custom names through this legacy helper.
    if (tolower(field) %in% tolower(c(pitax_header_builtin_fields(), "Record_Key", existing))) {
      stop(paste0("Metadata field is reserved or already exists: ", field), call. = FALSE)
    }
  }
  if (tolower(field) %in% tolower(c(pitax_header_builtin_fields(), "Record_Key", existing))) {
    stop(paste0("Metadata field is reserved or already exists: ", field), call. = FALSE)
  }
  field
}

pitax_header_register_field <- function(config, field, activate = TRUE) {
  pitax_header_create_custom_field(config, field, activate = activate)
}

pitax_header_activate_field <- function(config, field) {
  pitax_header_profile_add_field(config, field)
}

pitax_header_remove_from_draft <- function(config, field) {
  pitax_header_profile_remove_field(config, field)
}

pitax_header_registered_field_table <- function(config) {
  config <- pitax_header_coerce_config(config)
  active <- pitax_header_active_metadata_fields(config)
  # Union of custom fields + fields appearing in any profile + active draft.
  ids <- unique(c(
    names(config$custom_fields),
    unlist(lapply(config$profiles, function(p) p$ordered_fields), use.names = FALSE),
    config$draft$ordered_fields
  ))
  # Prefer showing user-facing custom + non-system library fields that were used.
  ids <- ids[nzchar(ids)]
  if (!length(ids)) {
    return(data.frame(Field = character(), Status = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    Field = vapply(ids, function(id) pitax_header_field_def(id, config)$display, character(1)),
    field_id = ids,
    Status = ifelse(ids %in% active, "In active header", "Available (not in header)"),
    stringsAsFactors = FALSE
  )
}

pitax_header_validate_template <- function(template, fields) {
  # Advanced/legacy path: still accept {Token} templates for migration previews.
  template <- as.character(template)[1]
  if (is.na(template) || !nzchar(trimws(template)) || grepl("[\r\n]", template)) {
    stop("FASTA header template must be one nonempty line.", call. = FALSE)
  }
  template <- sub("^[[:space:]]*>[[:space:]]*", "", template)
  if (grepl(">", template, fixed = TRUE)) {
    stop("FASTA header must not contain an extra '>' character.", call. = FALSE)
  }
  tokens <- regmatches(template, gregexpr("\\{[A-Za-z][A-Za-z0-9_]*\\}", template, perl = TRUE))[[1]]
  if (is.null(tokens)) tokens <- character()
  remainder <- gsub("\\{[A-Za-z][A-Za-z0-9_]*\\}", "", template, perl = TRUE)
  if (grepl("[{}]", remainder)) stop("Invalid token syntax; use {Field_Name}.", call. = FALSE)
  tokens <- unique(substring(tokens, 2, nchar(tokens) - 1))
  known <- unique(tolower(c(
    pitax_header_builtin_fields(), fields,
    names(pitax_header_field_registry()),
    names(pitax_header_alias_map()),
    vapply(pitax_header_field_registry(), function(f) f$display, character(1)),
    vapply(pitax_header_field_registry(), function(f) f$canonical_tag, character(1))
  )))
  known <- known[nzchar(known) & !is.na(known)]
  unknown <- tokens[!tolower(tokens) %in% known]
  if (length(unknown)) {
    stop(paste("Unknown FASTA header token(s):", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  template
}

pitax_header_referenced_fields <- function(template, fields) {
  validated <- pitax_header_validate_template(template, fields)
  tokens <- regmatches(validated, gregexpr("\\{[A-Za-z][A-Za-z0-9_]*\\}", validated, perl = TRUE))[[1]]
  if (!length(tokens)) return(character())
  unique(substring(tokens, 2, nchar(tokens) - 1L))
}

pitax_header_edit_table <- function(snap, config) {
  pitax_header_active_edit_table(snap, config)
}

pitax_header_builder_state_json <- function(config) {
  config <- pitax_header_coerce_config(config)
  lib <- pitax_header_library_table(config)
  list(
    draft = config$draft,
    profiles = lapply(config$profiles, function(p) list(id = p$id, name = p$name,
                                                        ordered_fields = I(p$ordered_fields),
                                                        tag_style = p$tag_style)),
    library = lib,
    active_fields = config$draft$ordered_fields,
    preview_template = config$draft_template,
    dirty = isTRUE(config$dirty),
    persistence = config$persistence,
    migration_warnings = config$migration_warnings,
    structure_issues = pitax_header_validate_profile_structure(config)
  )
}
