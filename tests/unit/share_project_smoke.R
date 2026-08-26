# PITAX 3.2.0 - Share Project snapshot unit/integration smokes (no Shiny UI).

source(file.path("R", "services", "share_project.R"), local = TRUE)
source(file.path("tests", "helpers", "test_paths.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

tmp_root <- file.path(tempdir(), paste0("pitax_share_", as.integer(Sys.time()), "_", sample.int(1e6, 1)))
dir.create(file.path(tmp_root, "data", "shared"), recursive = TRUE, showWarnings = FALSE)

minimal_project <- function(marker = "A", schema = 6L) {
  list(
    format = "SangerSequencePipelineProject",
    schema_version = schema,
    app_version = "3.2.0",
    saved_at = "2026-08-26 10:00:00",
    active_tab = "qc",
    ui_state = list(inspect_sample = "S1"),
    state = list(
      results = list(S1 = list(sample_id = "S1", seq = marker, raw_seq = marker, summary = data.frame())),
      summary = data.frame(sample_id = "S1", status = "OK", stringsAsFactors = FALSE),
      rename = NULL,
      settings = list(target = "ITS"),
      assay_profiles = data.frame(),
      project_defaults = list(),
      project_mode = "simple",
      read_assignments = data.frame(),
      architecture = NULL,
      consensus_set = list(records = list(), summary = data.frame()),
      multilocus_profile = list(profiles = data.frame(), evidence = data.frame()),
      migration_log = "",
      blast_jobs = data.frame(),
      blast_raw = list(),
      blast_ids = data.frame(),
      blast_hits = data.frame(),
      blast_batch_status_text = "",
      taxonomy_summary = data.frame(),
      taxonomy_hits = data.frame(),
      taxonomy_counts = data.frame(),
      taxonomy_status_text = "",
      taxonomy_batch_status_text = ""
    )
  )
}

# --- Create snapshot ---
created <- pitax_share_create(minimal_project("ORIG"), expiry_days = 3L, app_root = tmp_root,
                              app_version = "3.2.0", schema_version = 6L)
assert_true(isTRUE(created$ok), "Share create must succeed.")
assert_true(pitax_share_token_valid(created$token), "Created token must be 64 hex chars.")
assert_true(file.exists(created$path), "Snapshot file must exist under data/shared.")
assert_true(grepl("data", created$path, fixed = TRUE) && grepl("shared", created$path, fixed = TRUE),
            "Snapshot must live under data/shared.")
expires <- pitax_share_parse_time(created$snapshot$expires_at)
created_at <- pitax_share_parse_time(created$snapshot$created_at)
assert_true(is.finite(as.numeric(expires - created_at)), "expires_at must parse.")
# ~3 days (allow clock skew)
delta_hours <- as.numeric(difftime(expires, created_at, units = "hours"))
assert_true(delta_hours > 70 && delta_hours < 74, paste0("Default expiry should be ~72h, got ", delta_hours))

# --- Open snapshot ---
opened <- pitax_share_open(created$token, app_root = tmp_root, current_schema = 6L)
assert_true(isTRUE(opened$ok), "Valid token must open.")
assert_true(identical(opened$project$state$results$S1$seq, "ORIG"), "Opened project must restore marker.")

# --- Independent sessions: mutate copies, not the file ---
session_b <- opened$project
session_c <- pitax_share_open(created$token, app_root = tmp_root, current_schema = 6L)$project
session_b$state$results$S1$seq <- "CHANGED_B"
assert_true(identical(session_c$state$results$S1$seq, "ORIG"), "Session C must keep original after B mutates.")
reopened <- pitax_share_open(created$token, app_root = tmp_root, current_schema = 6L)
assert_true(identical(reopened$project$state$results$S1$seq, "ORIG"),
            "Reopening the link must load the immutable snapshot, not session B.")

# --- Expiry ---
expired_snap <- pitax_share_build_snapshot(
  minimal_project("OLD"), expiry_days = 1L, token = pitax_share_generate_token(),
  now = Sys.time() - (2 * 24 * 60 * 60), app_version = "3.2.0", schema_version = 6L
)
pitax_share_write_snapshot(expired_snap, app_root = tmp_root)
expired_open <- pitax_share_open(expired_snap$token, app_root = tmp_root, current_schema = 6L)
assert_true(!isTRUE(expired_open$ok) && identical(expired_open$status, "expired"),
            "Expired token must fail with expired status.")
assert_true(grepl("expired", expired_open$message, ignore.case = TRUE),
            "Expired message must mention expiry.")
assert_true(!file.exists(pitax_share_snapshot_path(expired_snap$token, app_root = tmp_root)) ||
              !file.exists(file.path(tmp_root, "data", "shared", paste0(expired_snap$token, ".rds"))),
            "Expired open should remove the snapshot file.")

# --- Invalid / traversal ---
bad1 <- pitax_share_open("../etc/passwd", app_root = tmp_root, current_schema = 6L)
bad2 <- pitax_share_open("zzzz", app_root = tmp_root, current_schema = 6L)
bad3 <- pitax_share_open(paste(rep("a", 64), collapse = ""), app_root = tmp_root, current_schema = 6L)
assert_true(!isTRUE(bad1$ok) && identical(bad1$status, "invalid_token"), "Path-like token must be rejected.")
assert_true(!isTRUE(bad2$ok) && identical(bad2$status, "invalid_token"), "Short token must be rejected.")
assert_true(!isTRUE(bad3$ok), "Unknown valid-format token must fail safely.")
assert_true(is.na(pitax_share_snapshot_path("../x", app_root = tmp_root)),
            "Traversal must not resolve to a path.")

# --- Cleanup keeps valid snapshots ---
keep_token <- created$token
removed <- pitax_share_cleanup_expired(app_root = tmp_root)
assert_true(file.exists(pitax_share_snapshot_path(keep_token, app_root = tmp_root)),
            "Cleanup must not delete a still-valid snapshot.")

# --- Schema too new ---
future <- minimal_project("FUT", schema = 99L)
fut <- pitax_share_create(future, expiry_days = 1L, app_root = tmp_root, app_version = "9.9.9", schema_version = 99L)
fut_open <- pitax_share_open(fut$token, app_root = tmp_root, current_schema = 6L)
assert_true(!isTRUE(fut_open$ok) && identical(fut_open$status, "schema_too_new"),
            "Newer schema snapshots must be refused.")

# --- Contracts in app sources ---
proj_text <- pitax_read_text("R", "server", "stages", "10_project.R")
ui_text <- pitax_read_text("R", "ui", "app_ui.R")
boot_text <- pitax_read_text("R", "app", "bootstrap.R")
assert_true(grepl("apply_loaded_project_object", proj_text, fixed = TRUE),
            "Load/share must share apply_loaded_project_object.")
assert_true(grepl("pitax_share_open", proj_text, fixed = TRUE) && grepl("parseQueryString", proj_text, fixed = TRUE),
            "Session init must open ?share=TOKEN via parseQueryString.")
assert_true(grepl("share_query_handled", proj_text, fixed = TRUE),
            "Share bootstrap must run inside a one-shot observe (reactive consumer), not bare onFlushed.")
assert_true(!grepl("session\\$onFlushed\\(function\\(\\) \\{[\\s\\S]*pitax_share_open", proj_text, perl = TRUE),
            "Share load must not call pitax_share_open from session$onFlushed (breaks rv access).")
assert_true(grepl('actionButton("share_project"', ui_text, fixed = TRUE),
            "UI must expose Share Project.")
assert_true(grepl("share_project.R", boot_text, fixed = TRUE),
            "bootstrap must source share_project.R.")
assert_true(!grepl("Stage 3 sequence gate is green", pitax_read_text("R", "server", "stages", "70_consensus.R"), fixed = TRUE),
            "Internal consensus success toast must remain removed.")

unlink(tmp_root, recursive = TRUE, force = TRUE)
cat("Share project snapshot smokes passed.\n")
