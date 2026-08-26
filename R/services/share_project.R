# ============================================================
# PITAX share-project snapshots (immutable server-side copies)
# Tokens are possession-based; no accounts/ACLs in v3.2.0.
# ============================================================

PITAX_SHARE_TOKEN_BYTES <- 32L
PITAX_SHARE_TOKEN_HEX_LEN <- 64L
PITAX_SHARE_DEFAULT_DAYS <- 3L
PITAX_SHARE_ALLOWED_DAYS <- c(1L, 3L, 7L)

pitax_share_root <- function(app_root = getwd()) {
  normalizePath(file.path(app_root, "data", "shared"), winslash = "/", mustWork = FALSE)
}

pitax_share_ensure_dir <- function(app_root = getwd()) {
  root <- pitax_share_root(app_root)
  if (!dir.exists(root)) dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) stop("Could not create share storage directory.", call. = FALSE)
  root
}

pitax_share_generate_token <- function() {
  bytes <- as.raw(sample.int(256L, PITAX_SHARE_TOKEN_BYTES, replace = TRUE) - 1L)
  paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

pitax_share_token_valid <- function(token) {
  token <- as.character(token)[1]
  if (!length(token) || is.na(token) || !nzchar(token)) return(FALSE)
  grepl(paste0("^[a-f0-9]{", PITAX_SHARE_TOKEN_HEX_LEN, "}$"), token)
}

# Resolve a token to a path inside the share root only (blocks traversal).
pitax_share_snapshot_path <- function(token, app_root = getwd()) {
  if (!pitax_share_token_valid(token)) return(NA_character_)
  root <- pitax_share_ensure_dir(app_root)
  candidate <- normalizePath(file.path(root, paste0(token, ".rds")), winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  # Ensure candidate stays under root (Windows-safe prefix check).
  if (!startsWith(tolower(candidate), tolower(paste0(root_norm, "/")))) {
    return(NA_character_)
  }
  candidate
}

pitax_share_parse_expiry_days <- function(value, default = PITAX_SHARE_DEFAULT_DAYS) {
  days <- suppressWarnings(as.integer(value)[1])
  if (!is.finite(days) || !(days %in% PITAX_SHARE_ALLOWED_DAYS)) return(as.integer(default))
  days
}

pitax_share_format_time <- function(when = Sys.time()) {
  format(when, "%Y-%m-%d %H:%M:%S", tz = "UTC", usetz = TRUE)
}

pitax_share_parse_time <- function(value) {
  value <- as.character(value)[1]
  if (!nzchar(value)) return(as.POSIXct(NA, tz = "UTC"))
  # Accept both "UTC" suffix and plain timestamps written as UTC.
  value <- sub(" UTC$", "", value, ignore.case = TRUE)
  value <- trimws(value)
  suppressWarnings(as.POSIXct(value, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
}

pitax_share_build_snapshot <- function(project_bundle, expiry_days = PITAX_SHARE_DEFAULT_DAYS,
                                       token = pitax_share_generate_token(),
                                       now = Sys.time(),
                                       app_version = NULL,
                                       schema_version = NULL) {
  if (!is.list(project_bundle) || !identical(project_bundle$format, "SangerSequencePipelineProject")) {
    stop("Share snapshot requires a valid SangerSequencePipelineProject bundle.", call. = FALSE)
  }
  expiry_days <- pitax_share_parse_expiry_days(expiry_days)
  created <- as.POSIXct(now, tz = "UTC")
  expires <- created + (expiry_days * 24 * 60 * 60)
  list(
    format = "PITAXShareSnapshot",
    snapshot_id = paste0("snap_", token),
    token = token,
    created_at = pitax_share_format_time(created),
    expires_at = pitax_share_format_time(expires),
    pitax_version = if (is.null(app_version)) as.character(project_bundle$app_version) else as.character(app_version),
    schema_version = if (is.null(schema_version)) {
      suppressWarnings(as.integer(project_bundle$schema_version)[1])
    } else {
      suppressWarnings(as.integer(schema_version)[1])
    },
    project = project_bundle
  )
}

pitax_share_write_snapshot <- function(snapshot, app_root = getwd()) {
  if (!is.list(snapshot) || !identical(snapshot$format, "PITAXShareSnapshot")) {
    stop("Invalid share snapshot object.", call. = FALSE)
  }
  if (!pitax_share_token_valid(snapshot$token)) stop("Invalid share token.", call. = FALSE)
  path <- pitax_share_snapshot_path(snapshot$token, app_root = app_root)
  if (!nzchar(path) || is.na(path)) stop("Refusing to write share snapshot outside the share root.", call. = FALSE)
  saveRDS(snapshot, file = path, compress = "xz")
  invisible(path)
}

pitax_share_create <- function(project_bundle, expiry_days = PITAX_SHARE_DEFAULT_DAYS,
                               app_root = getwd(), app_version = NULL, schema_version = NULL) {
  pitax_share_cleanup_expired(app_root = app_root)
  token <- pitax_share_generate_token()
  # Extremely unlikely collision; regenerate once if needed.
  path <- pitax_share_snapshot_path(token, app_root = app_root)
  if (file.exists(path)) token <- pitax_share_generate_token()
  snapshot <- pitax_share_build_snapshot(
    project_bundle,
    expiry_days = expiry_days,
    token = token,
    app_version = app_version,
    schema_version = schema_version
  )
  written <- pitax_share_write_snapshot(snapshot, app_root = app_root)
  list(ok = TRUE, token = token, path = written, snapshot = snapshot)
}

pitax_share_is_expired <- function(snapshot, now = Sys.time()) {
  expires <- pitax_share_parse_time(snapshot$expires_at)
  if (is.na(expires)) return(TRUE)
  isTRUE(as.POSIXct(now, tz = "UTC") > expires)
}

# Returns list(ok=, status=, message=, snapshot=, project=)
pitax_share_open <- function(token, app_root = getwd(), now = Sys.time(),
                             current_schema = NULL) {
  if (!pitax_share_token_valid(token)) {
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(ok = FALSE, status = "invalid_token",
                message = "This share link is invalid.", snapshot = NULL, project = NULL))
  }
  path <- pitax_share_snapshot_path(token, app_root = app_root)
  if (!nzchar(path) || is.na(path) || !file.exists(path)) {
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(ok = FALSE, status = "not_found",
                message = "This shared PITAX project was not found.", snapshot = NULL, project = NULL))
  }
  snap <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(snap) || !identical(snap$format, "PITAXShareSnapshot") || is.null(snap$project)) {
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(ok = FALSE, status = "corrupt",
                message = "This shared PITAX project could not be read safely.", snapshot = NULL, project = NULL))
  }
  if (!identical(as.character(snap$token), as.character(token))) {
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(ok = FALSE, status = "token_mismatch",
                message = "This shared PITAX project could not be validated.", snapshot = NULL, project = NULL))
  }
  if (pitax_share_is_expired(snap, now = now)) {
    unlink(path, force = TRUE)
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(
      ok = FALSE, status = "expired",
      message = paste0(
        "This shared PITAX project has expired.\n",
        "It was available until ", as.character(snap$expires_at), "."
      ),
      snapshot = snap, project = NULL
    ))
  }
  project <- snap$project
  if (!is.list(project) || !identical(project$format, "SangerSequencePipelineProject")) {
    pitax_share_cleanup_expired(app_root = app_root, now = now)
    return(list(ok = FALSE, status = "bad_project",
                message = "This shared snapshot does not contain a valid PITAX project.", snapshot = snap, project = NULL))
  }
  source_schema <- suppressWarnings(as.integer(project$schema_version)[1])
  if (!is.null(current_schema)) {
    current_schema <- suppressWarnings(as.integer(current_schema)[1])
    if (is.finite(source_schema) && is.finite(current_schema) && source_schema > current_schema) {
      pitax_share_cleanup_expired(app_root = app_root, now = now)
      return(list(
        ok = FALSE, status = "schema_too_new",
        message = "This shared project was created by a newer project schema and cannot be loaded safely.",
        snapshot = snap, project = NULL
      ))
    }
  }
  pitax_share_cleanup_expired(app_root = app_root, now = now)
  list(ok = TRUE, status = "ok", message = "Shared project loaded.", snapshot = snap, project = project)
}

pitax_share_cleanup_expired <- function(app_root = getwd(), now = Sys.time()) {
  root <- pitax_share_root(app_root)
  if (!dir.exists(root)) return(invisible(0L))
  files <- list.files(root, pattern = "[.]rds$", full.names = TRUE)
  removed <- 0L
  for (f in files) {
    token <- sub("[.]rds$", "", basename(f))
    if (!pitax_share_token_valid(token)) {
      # Unexpected name in share root - remove for safety.
      unlink(f, force = TRUE)
      removed <- removed + 1L
      next
    }
    # Confirm path still under root.
    safe <- pitax_share_snapshot_path(token, app_root = app_root)
    if (is.na(safe) || !identical(normalizePath(f, winslash = "/", mustWork = FALSE), safe)) {
      next
    }
    snap <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.list(snap) || !identical(snap$format, "PITAXShareSnapshot") || pitax_share_is_expired(snap, now = now)) {
      unlink(f, force = TRUE)
      removed <- removed + 1L
    }
  }
  invisible(removed)
}

pitax_share_public_url <- function(token, session = NULL, base_url = NULL) {
  if (!pitax_share_token_valid(token)) return("")
  if (!is.null(base_url) && nzchar(as.character(base_url)[1])) {
    base <- sub("/+$", "", as.character(base_url)[1])
    return(paste0(base, "/?share=", token))
  }
  if (is.null(session)) return(paste0("?share=", token))
  proto <- as.character(session$clientData$url_protocol)[1]
  host <- as.character(session$clientData$url_hostname)[1]
  port <- as.character(session$clientData$url_port)[1]
  path <- as.character(session$clientData$url_pathname)[1]
  if (!nzchar(proto)) proto <- "http:"
  if (!nzchar(host)) host <- "localhost"
  if (!nzchar(path) || path == "/") path <- "/"
  # Shiny apps often live at a subpath; keep pathname and append query.
  port_part <- if (nzchar(port) && !port %in% c("80", "443")) paste0(":", port) else ""
  # pathname may already end with /; normalize to .../?share=
  path <- sub("/+$", "", path)
  if (!nzchar(path)) path <- ""
  paste0(proto, "//", host, port_part, path, "/?share=", token)
}
