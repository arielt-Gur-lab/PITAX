# Custom FASTA headers: inert template expansion and project-scoped metadata.
# Never change source/biological identity or construct code from templates.

pitax_header_builtin_fields <- function() {
  c("Final_Name", "Isolate", "Locus", "Direction", "Source_ID", "Record_ID")
}

pitax_header_default_template <- function() {
  '{Final_Name} organism={Organism} experiment="{Experiment}" date={Seq_Date}'
}

pitax_header_default_config <- function() {
  list(
    fields = c("Organism", "Experiment", "Seq_Date"),
    values = data.frame(Record_Key = character(), Organism = character(),
                        Experiment = character(), Seq_Date = character(),
                        stringsAsFactors = FALSE),
    templates = list("Research metadata" = pitax_header_default_template()),
    active_template = "Research metadata",
    draft_template = pitax_header_default_template()
  )
}

pitax_header_validate_field <- function(field, existing = character()) {
  field <- trimws(as.character(field)[1])
  if (is.na(field) || !grepl("^[A-Za-z][A-Za-z0-9_]*$", field)) {
    stop("Metadata field names must start with a letter and contain only letters, digits or underscores.", call. = FALSE)
  }
  if (tolower(field) %in% tolower(c(pitax_header_builtin_fields(), existing))) {
    stop(paste0("Metadata field is reserved or already exists: ", field), call. = FALSE)
  }
  field
}

pitax_header_coerce_config <- function(config) {
  default <- pitax_header_default_config()
  if (!is.list(config)) return(default)
  fields <- as.character(config$fields %||% default$fields)
  if (length(fields) > 40L || anyNA(fields) || anyDuplicated(tolower(fields)) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", fields)) ||
      any(tolower(fields) %in% tolower(pitax_header_builtin_fields()))) {
    stop("Invalid custom FASTA metadata field definitions.", call. = FALSE)
  }
  values <- config$values
  if (!is.data.frame(values) || !"Record_Key" %in% names(values)) values <- default$values
  for (field in fields) if (!field %in% names(values)) values[[field]] <- rep("", nrow(values))
  values <- values[, c("Record_Key", fields), drop = FALSE]
  for (col in names(values)) {
    values[[col]] <- as.character(values[[col]])
    values[[col]][is.na(values[[col]])] <- ""
  }
  values <- values[!duplicated(values$Record_Key) & nzchar(values$Record_Key), , drop = FALSE]
  rownames(values) <- NULL
  templates <- config$templates
  if (!is.list(templates) || !length(templates) || is.null(names(templates)) ||
      any(!nzchar(names(templates))) || anyDuplicated(names(templates))) {
    templates <- default$templates
  }
  templates <- lapply(templates, function(x) as.character(x)[1])
  active <- as.character(config$active_template %||% default$active_template)[1]
  if (is.na(active) || !active %in% names(templates)) active <- names(templates)[1]
  draft <- as.character(config$draft_template %||% templates[[active]])[1]
  if (is.na(draft)) draft <- templates[[active]]
  list(fields = fields, values = values, templates = templates,
       active_template = active, draft_template = draft)
}

pitax_header_record_key <- function(component_id, record_id) {
  component_id <- as.character(component_id)[1]
  namespace <- if (identical(component_id, "qc_trimmed_sequences")) "read" else
    if (component_id %in% c("consensus_analysis_sequences", "blast_query_sequences")) "analysis" else
    if (identical(component_id, "multilocus_sequences")) "multilocus" else
      stop(paste("Unknown sequence component:", component_id), call. = FALSE)
  paste0(namespace, "::", as.character(record_id)[1])
}

pitax_header_validate_template <- function(template, fields) {
  template <- as.character(template)[1]
  if (is.na(template) || !nzchar(trimws(template)) || grepl("[\r\n]", template)) {
    stop("FASTA header template must be one nonempty line.", call. = FALSE)
  }
  template <- sub("^[[:space:]]*>[[:space:]]*", "", template)
  if (grepl(">", template, fixed = TRUE)) stop("FASTA header must not contain an extra '>' character.", call. = FALSE)
  tokens <- regmatches(template, gregexpr("\\{[A-Za-z][A-Za-z0-9_]*\\}", template, perl = TRUE))[[1]]
  if (is.null(tokens)) tokens <- character()
  remainder <- gsub("\\{[A-Za-z][A-Za-z0-9_]*\\}", "", template, perl = TRUE)
  if (grepl("[{}]", remainder)) stop("Invalid token syntax; use {Field_Name}.", call. = FALSE)
  tokens <- unique(substring(tokens, 2, nchar(tokens) - 1))
  unknown <- setdiff(tokens, c(pitax_header_builtin_fields(), fields))
  if (length(unknown)) stop(paste("Unknown FASTA header token(s):", paste(unknown, collapse = ", ")), call. = FALSE)
  template
}

pitax_header_catalog <- function(snap) {
  rows <- list()
  push <- function(key, final_name, isolate = "", locus = "", direction = "", source = "", record = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      Record_Key = key, Final_Name = final_name, Isolate = isolate, Locus = locus,
      Direction = direction, Source_ID = source, Record_ID = record, stringsAsFactors = FALSE
    )
  }
  asg <- snap$read_assignments
  if (is.data.frame(asg) && nrow(asg)) {
    for (i in seq_len(nrow(asg))) {
      sid <- as.character(asg$Source_ID[i])
      name <- pitax_export_final_name(snap, sid)
      push(paste0("read::", sid), name, as.character(asg$Isolate[i]),
           as.character(asg$Locus[i]), as.character(asg$Direction[i]), sid, sid)
    }
  } else if (length(snap$results)) {
    for (sid in names(snap$results)) push(paste0("read::", sid), pitax_export_final_name(snap, sid), source = sid, record = sid)
  }
  cs <- snap$consensus_set
  if (is.list(cs) && is.list(cs$records)) {
    for (id in names(cs$records)) {
      rec <- cs$records[[id]]
      push(paste0("analysis::", id), as.character(rec$final_name %||% id),
           as.character(rec$isolate %||% ""), as.character(rec$locus %||% ""),
           "oriented", paste(as.character(rec$source_read_ids %||% ""), collapse = ","), id)
    }
  }
  ev <- tryCatch(pitax_export_extract_multilocus_sequences(snap), error = function(e) list())
  if (length(ev)) for (id in names(ev)) {
    rec <- ev[[id]]
    push(paste0("multilocus::", id), as.character(rec$final_name %||% id),
         as.character(rec$isolate %||% ""), as.character(rec$locus %||% ""),
         "oriented", "", id)
  }
  if (!length(rows)) return(data.frame(
    Record_Key = character(), Final_Name = character(), Isolate = character(),
    Locus = character(), Direction = character(), Source_ID = character(),
    Record_ID = character(), stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$Record_Key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

pitax_header_edit_table <- function(snap, config) {
  config <- pitax_header_coerce_config(config)
  df <- pitax_header_catalog(snap)
  for (field in config$fields) {
    idx <- match(df$Record_Key, config$values$Record_Key)
    df[[field]] <- config$values[[field]][idx]
    df[[field]][is.na(df[[field]])] <- ""
  }
  df
}

pitax_header_set_value <- function(config, record_keys, field, value) {
  config <- pitax_header_coerce_config(config)
  if (!field %in% config$fields) stop("Only metadata fields are editable.", call. = FALSE)
  value <- as.character(value)[1]
  if (is.na(value) || grepl("[\r\n>]", value)) {
    stop("Header values must not contain line breaks or FASTA '>' symbols.", call. = FALSE)
  }
  for (key in unique(as.character(record_keys))) {
    if (!nzchar(key)) next
    idx <- match(key, config$values$Record_Key)
    if (is.na(idx)) {
      row <- as.data.frame(as.list(setNames(rep("", length(config$fields) + 1L),
                                            c("Record_Key", config$fields))), stringsAsFactors = FALSE)
      row$Record_Key <- key
      config$values <- rbind(config$values, row)
      idx <- nrow(config$values)
    }
    config$values[[field]][idx] <- value
  }
  config
}

pitax_header_render <- function(rec, component_id, record_id, config) {
  config <- pitax_header_coerce_config(config)
  template <- pitax_header_validate_template(config$draft_template, config$fields)
  key <- pitax_header_record_key(component_id, record_id)
  values <- c(
    Final_Name = as.character(rec$final_name %||% record_id)[1],
    Isolate = as.character(rec$isolate %||% "")[1],
    Locus = as.character(rec$locus %||% "")[1],
    Direction = as.character(rec$direction %||% "")[1],
    Source_ID = if (startsWith(key, "read::")) as.character(record_id)[1] else "",
    Record_ID = as.character(record_id)[1]
  )
  idx <- match(key, config$values$Record_Key)
  for (field in config$fields) {
    values[field] <- if (is.na(idx)) "" else config$values[[field]][idx]
  }
  matches <- gregexpr("\\{[A-Za-z][A-Za-z0-9_]*\\}", template, perl = TRUE)[[1]]
  if (identical(matches[1], -1L)) return(template)
  lengths <- attr(matches, "match.length")
  parts <- character()
  cursor <- 1L
  for (i in seq_along(matches)) {
    parts <- c(parts, substr(template, cursor, matches[i] - 1L))
    token <- substr(template, matches[i] + 1L, matches[i] + lengths[i] - 2L)
    value <- unname(values[token])
    if (!length(value) || is.na(value) || !nzchar(trimws(value))) {
      stop(paste0("Missing {", token, "} for ", as.character(rec$final_name %||% record_id)), call. = FALSE)
    }
    if (grepl("[\r\n>]", value)) stop(paste0("Invalid header value for {", token, "}"), call. = FALSE)
    parts <- c(parts, value)
    cursor <- matches[i] + lengths[i]
  }
  header <- paste0(c(parts, substr(template, cursor, nchar(template))), collapse = "")
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
                  ". Include {Final_Name} to distinguish reads/loci."), call. = FALSE)
    }
    ids <- c(ids, primary)
    records[[record_id]]$custom_header <- rendered
  }
  records
}
