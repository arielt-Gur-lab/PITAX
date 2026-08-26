# PITAX 3.1.0 - BLAST RTOE auto-poll unit smokes (no NCBI HTTP).

source(file.path("R", "services", "blast_polling.R"), local = TRUE)
source(file.path("tests", "helpers", "test_paths.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

# --- 1. RTOE scheduling ---
assert_true(identical(blast_parse_rtoe_seconds("30"), 30L), "RTOE parse must read integer seconds.")
assert_true(identical(blast_rtoe_floor_seconds(30), 30L), "RTOE floor for 30 is 30.")
assert_true(identical(blast_rtoe_floor_seconds(10), 30L), "RTOE floor must not go below 30s.")
assert_true(identical(blast_first_auto_poll_delay_seconds(30), 35L),
            "First auto poll delay is RTOE floor + 5s buffer.")
assert_true(identical(blast_next_wait_after_attempt(1L), 60L), "After attempt 1 wait 60s.")
assert_true(identical(blast_next_wait_after_attempt(2L), 90L), "After attempt 2 wait 90s.")
assert_true(is.na(blast_next_wait_after_attempt(3L)), "After attempt 3 there is no further auto schedule.")

t0 <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
earliest <- blast_schedule_next_poll_at(t0, blast_rtoe_floor_seconds(30))
first_auto <- blast_schedule_next_poll_at(t0, blast_first_auto_poll_delay_seconds(30))
assert_true(identical(earliest, format(t0 + 30, "%Y-%m-%d %H:%M:%S")),
            "earliest_retrieve_at = submit + RTOE floor.")
assert_true(identical(first_auto, format(t0 + 35, "%Y-%m-%d %H:%M:%S")),
            "next_poll_at first auto = submit + floor + buffer.")

# UI refresh interval (2s) must not equal NCBI schedule delays.
assert_true(blast_first_auto_poll_delay_seconds(30) > 2L, "First NCBI check is not UI-refresh cadence.")

# --- Manual RTOE floor ---
job0 <- data.frame(
  status = "SUBMITTED", rid = "RID1", rtoe = "30",
  submitted_at = format(t0, "%Y-%m-%d %H:%M:%S"),
  earliest_retrieve_at = earliest,
  next_poll_at = first_auto,
  auto_poll_enabled = TRUE, auto_poll_attempts = 0L,
  manual_retrieval_required = FALSE,
  stringsAsFactors = FALSE
)
before <- t0 + 20
after_floor <- t0 + 31
before_auto <- t0 + 32
assert_true(!blast_manual_retrieve_allowed(job0, before), "Manual Check now blocked before RTOE floor.")
assert_true(blast_manual_retrieve_allowed(job0, after_floor), "Manual Check now allowed after RTOE floor.")
assert_true(!blast_auto_poll_due(job0, before_auto), "Auto poll still waits for buffer after RTOE floor.")
assert_true(blast_auto_poll_due(job0, t0 + 36), "Auto poll due after next_poll_at.")

# --- 2. READY on first poll ---
u1 <- blast_job_after_auto_poll(job0, "READY")
assert_true(identical(u1$auto_poll_attempts, 1L), "READY on attempt 1 sets attempts=1.")
assert_true(isFALSE(u1$auto_poll_enabled), "READY stops auto poll.")
assert_true(identical(u1$next_poll_at, ""), "READY clears next_poll_at.")

# --- 3. WAITING then READY ---
u_w1 <- blast_job_after_auto_poll(job0, "WAITING", now = after_floor)
assert_true(identical(u_w1$auto_poll_attempts, 1L), "First WAITING increments to 1.")
assert_true(isTRUE(u_w1$auto_poll_enabled), "WAITING attempt 1 keeps auto enabled.")
job1 <- job0
job1$auto_poll_attempts <- u_w1$auto_poll_attempts
job1$next_poll_at <- u_w1$next_poll_at
job1$status <- "WAITING"
u_r2 <- blast_job_after_auto_poll(job1, "READY")
assert_true(identical(u_r2$auto_poll_attempts, 2L), "READY after WAITING records attempt 2.")
assert_true(isFALSE(u_r2$auto_poll_enabled), "READY after WAITING stops auto.")

# --- 4. Three WAITING ---
job_a <- job0
u_a <- blast_job_after_auto_poll(job_a, "WAITING")
job_a$auto_poll_attempts <- u_a$auto_poll_attempts
u_b <- blast_job_after_auto_poll(job_a, "WAITING")
job_a$auto_poll_attempts <- u_b$auto_poll_attempts
u_c <- blast_job_after_auto_poll(job_a, "WAITING")
assert_true(identical(u_c$auto_poll_attempts, 3L), "Third WAITING sets attempts=3.")
assert_true(isFALSE(u_c$auto_poll_enabled), "Third WAITING disables auto.")
assert_true(isTRUE(u_c$manual_retrieval_required), "Third WAITING requires manual retrieval.")
assert_true(identical(u_c$next_poll_at, ""), "Third WAITING clears next_poll_at (no 4th auto).")
job_stopped <- job_a
job_stopped$auto_poll_enabled <- FALSE
job_stopped$manual_retrieval_required <- TRUE
job_stopped$auto_poll_attempts <- 3L
job_stopped$next_poll_at <- ""
assert_true(!blast_auto_poll_due(job_stopped, Sys.time()), "Stopped job is never auto-due.")

# --- 5. Manual after stop does not re-enable auto ---
# Manual path applies blast_disable_auto_poll_fields, not after_auto_poll with enable.
dis <- blast_disable_auto_poll_fields()
assert_true(isFALSE(dis$auto_poll_enabled), "Manual terminal path keeps auto disabled.")
assert_true(isFALSE(dis$manual_retrieval_required), "Disable helper clears manual flag on terminal complete.")

# --- 6. Stale contract markers ---
blast_text <- pitax_read_text("R", "server", "stages", "90_blast.R")
assert_true(grepl("job_revision", blast_text, fixed = TRUE) && grepl("current_revision", blast_text, fixed = TRUE),
            "BLAST retrieve must compare stored and current consensus revisions.")
assert_true(grepl("blast_disable_auto_poll_fields", blast_text, fixed = TRUE),
            "STALE / terminal paths must disable auto-poll fields.")

# --- 7. Distinct failure states ---
u_fail <- blast_job_after_auto_poll(job0, "FAILED")
u_unk <- blast_job_after_auto_poll(job0, "UNKNOWN")
u_net <- blast_job_after_auto_poll(job0, "ERROR")
assert_true(isFALSE(u_fail$auto_poll_enabled) && isFALSE(u_unk$auto_poll_enabled),
            "FAILED/UNKNOWN stop auto without marking biological timeout.")
assert_true(isTRUE(u_net$auto_poll_enabled) && identical(u_net$auto_poll_attempts, 1L),
            "Network ERROR counts an attempt but keeps RID/auto until budget exhausted.")
assert_true(!identical(u_fail, u_net), "FAILED and NETWORK transitions differ.")

banner_manual <- blast_status_banner_parts(
  data.frame(
    status = "WAITING", rid = "RID1", rtoe = "30",
    submitted_at = format(t0, "%Y-%m-%d %H:%M:%S"),
    earliest_retrieve_at = earliest,
    next_poll_at = "",
    auto_poll_enabled = FALSE, auto_poll_attempts = 3L,
    manual_retrieval_required = TRUE,
    database = "core_nt", hitlist_size = 25L,
    last_checked_at = format(t0 + 300, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
)
assert_true(identical(banner_manual$kind, "manual_required"), "Banner marks manual_required after 3 checks.")
assert_true(any(grepl("Automatic checks completed:|Check now", banner_manual$lines)),
            "Stopped-after-3 copy must guide the user to Check now.")
assert_true(!any(grepl("^Status: FAILED", banner_manual$lines)),
            "WAITING-after-3 must not be labeled FAILED.")
assert_true(!any(grepl("First NCBI check in:|Next automatic check:", banner_manual$lines)),
            "Banner must not show automatic-check countdowns.")

banner_armed <- blast_status_banner_parts(job0)
assert_true(identical(banner_armed$kind, "auto_active"), "Waiting auto job uses auto_active banner.")
assert_true(any(grepl("Waiting for the first automatic NCBI check \\(after RTOE\\)", banner_armed$lines)),
            "First-wait banner must use the RTOE waiting sentence.")
assert_true(!any(grepl("armed|Automatic retrieval:", banner_armed$lines, ignore.case = TRUE)),
            "Banner must not use armed/retrieval jargon.")
assert_true(!any(grepl("First NCBI check in:|Next automatic check:", banner_armed$lines)),
            "Armed banner must not show countdown timers.")

# --- 8. UI / network separation contract markers ---
assert_true(grepl("PITAX auto-poll observer \\(NCBI\\)", blast_text),
            "Auto-poll observer must be explicitly marked.")
assert_true(grepl("UI refresh timer only", blast_text, fixed = TRUE) ||
              grepl("does not call retrieve_blast_job", blast_text, fixed = TRUE),
            "Status-card invalidateLater must be documented as UI-only.")
assert_true(grepl('mode = "auto"', blast_text, fixed = TRUE) && grepl('mode = "manual"', blast_text, fixed = TRUE),
            "retrieve_blast_job must distinguish auto vs manual modes.")
assert_true(grepl("blast_auto_check", blast_text, fixed = TRUE) &&
              grepl("Automatic NCBI check in progress", blast_text, fixed = TRUE),
            "Auto-poll must show an in-progress notification before the NCBI Get.")

load_text <- pitax_read_text("R", "server", "stages", "10_project.R")
assert_true(grepl("Safe reload", load_text, fixed = TRUE) &&
              grepl("auto_poll_enabled\\[pending\\] <- FALSE", load_text),
            "Project load must disable auto-poll for pending RIDs.")

cat("BLAST RTOE auto-poll smokes passed.\n")
