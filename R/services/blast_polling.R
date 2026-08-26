# ============================================================
# BLAST RTOE scheduling and auto-poll state transitions (pure R)
# No HTTP / Shiny side effects - unit-tested offline.
# ============================================================

blast_parse_rtoe_seconds <- function(rtoe, default = 30L) {
  secs <- suppressWarnings(as.integer(as.character(rtoe)[1]))
  if (!length(secs) || !is.finite(secs) || secs < 0L) return(as.integer(default))
  secs
}

# Minimum wait before any NCBI retrieve (manual Check now or automatic).
blast_rtoe_floor_seconds <- function(rtoe, min_seconds = 30L) {
  max(blast_parse_rtoe_seconds(rtoe), as.integer(min_seconds))
}

# First automatic poll: RTOE floor + small safety buffer (RTOE is not a poll interval).
blast_first_auto_poll_delay_seconds <- function(rtoe, min_seconds = 30L, buffer_seconds = 5L) {
  blast_rtoe_floor_seconds(rtoe, min_seconds) + as.integer(buffer_seconds)
}

# Delay until the next automatic attempt after a WAITING result.
# attempt is the count *after* increment (1 -> wait 60s for attempt 2; 2 -> 90s for attempt 3).
blast_next_wait_after_attempt <- function(attempt) {
  attempt <- suppressWarnings(as.integer(attempt)[1])
  if (!is.finite(attempt) || attempt < 1L) return(NA_integer_)
  if (attempt == 1L) return(60L)
  if (attempt == 2L) return(90L)
  NA_integer_
}

blast_format_timestamp <- function(when = Sys.time()) {
  format(when, "%Y-%m-%d %H:%M:%S")
}

blast_parse_timestamp <- function(value) {
  value <- as.character(value)[1]
  if (!nzchar(value) || identical(value, "NA")) return(as.POSIXct(NA))
  suppressWarnings(as.POSIXct(value, format = "%Y-%m-%d %H:%M:%S"))
}

blast_schedule_next_poll_at <- function(from_time, delay_seconds) {
  delay_seconds <- suppressWarnings(as.integer(delay_seconds)[1])
  if (!is.finite(delay_seconds) || delay_seconds < 0L) return("")
  if (inherits(from_time, "POSIXt")) {
    t0 <- from_time
  } else {
    t0 <- blast_parse_timestamp(from_time)
  }
  if (is.na(t0)) t0 <- Sys.time()
  blast_format_timestamp(t0 + delay_seconds)
}

blast_seconds_until <- function(timestamp_text, now = Sys.time()) {
  t1 <- blast_parse_timestamp(timestamp_text)
  if (is.na(t1)) return(NA_real_)
  as.numeric(difftime(t1, now, units = "secs"))
}

# RTOE floor: manual Check now (and auto) must not contact NCBI before this.
blast_manual_retrieve_allowed <- function(job_row, now = Sys.time()) {
  job_row <- job_row[1, , drop = FALSE]
  earliest <- ""
  if ("earliest_retrieve_at" %in% names(job_row)) {
    earliest <- as.character(job_row$earliest_retrieve_at[1])
  }
  if (!nzchar(earliest)) {
    # Legacy rows: derive from submitted_at + RTOE floor.
    submitted <- blast_parse_timestamp(job_row$submitted_at[1])
    if (is.na(submitted)) return(TRUE)
    floor_secs <- blast_rtoe_floor_seconds(job_row$rtoe[1])
    earliest_t <- submitted + floor_secs
    return(isTRUE(now >= earliest_t))
  }
  t0 <- blast_parse_timestamp(earliest)
  if (is.na(t0)) return(TRUE)
  isTRUE(now >= t0)
}

blast_auto_poll_due <- function(job_row, now = Sys.time()) {
  job_row <- job_row[1, , drop = FALSE]
  if (!isTRUE(as.logical(job_row$auto_poll_enabled[1]))) return(FALSE)
  status <- as.character(job_row$status[1])
  if (!status %in% c("SUBMITTED", "WAITING")) return(FALSE)
  if (!blast_manual_retrieve_allowed(job_row, now)) return(FALSE)
  next_at <- if ("next_poll_at" %in% names(job_row)) as.character(job_row$next_poll_at[1]) else ""
  if (!nzchar(next_at)) return(FALSE)
  t1 <- blast_parse_timestamp(next_at)
  if (is.na(t1)) return(FALSE)
  isTRUE(now >= t1)
}

# Pure state transition after one automatic poll result. Does not mutate callers.
blast_job_after_auto_poll <- function(job_row, poll_status, now = Sys.time()) {
  job_row <- job_row[1, , drop = FALSE]
  attempts <- suppressWarnings(as.integer(job_row$auto_poll_attempts[1]))
  if (!is.finite(attempts)) attempts <- 0L
  status <- toupper(as.character(poll_status)[1])

  if (identical(status, "TOO_SOON")) {
    return(list())
  }

  if (identical(status, "STALE")) {
    return(list(
      auto_poll_enabled = FALSE,
      next_poll_at = "",
      manual_retrieval_required = FALSE
    ))
  }

  if (status %in% c("READY", "FAILED", "UNKNOWN", "NO_HITS")) {
    return(list(
      auto_poll_enabled = FALSE,
      auto_poll_attempts = attempts + 1L,
      next_poll_at = "",
      manual_retrieval_required = FALSE
    ))
  }

  # WAITING or network/contact ERROR: count toward the 3-attempt budget.
  if (status %in% c("WAITING", "ERROR", "NETWORK")) {
    new_attempts <- attempts + 1L
    if (new_attempts >= 3L) {
      return(list(
        auto_poll_enabled = FALSE,
        auto_poll_attempts = new_attempts,
        next_poll_at = "",
        manual_retrieval_required = TRUE
      ))
    }
    delay <- if (identical(status, "WAITING")) {
      blast_next_wait_after_attempt(new_attempts)
    } else {
      60L
    }
    return(list(
      auto_poll_enabled = TRUE,
      auto_poll_attempts = new_attempts,
      next_poll_at = blast_schedule_next_poll_at(now, delay),
      manual_retrieval_required = FALSE
    ))
  }

  list()
}

blast_disable_auto_poll_fields <- function() {
  list(
    auto_poll_enabled = FALSE,
    next_poll_at = "",
    manual_retrieval_required = FALSE
  )
}

# UI copy helpers (no HTML). Elapsed/countdown formatting is left to the stage.
blast_status_banner_parts <- function(job_row, now = Sys.time(), format_elapsed = NULL) {
  job_row <- job_row[1, , drop = FALSE]
  if (is.null(format_elapsed)) {
    format_elapsed <- function(secs) {
      secs <- suppressWarnings(as.integer(round(secs)))
      if (!length(secs) || !is.finite(secs[1]) || secs[1] < 0) return("unknown")
      secs <- secs[1]
      if (secs < 60L) return(paste0(secs, "s"))
      mins <- secs %/% 60L
      rem <- secs %% 60L
      if (mins < 60L) return(sprintf("%dm %ds", mins, rem))
      hours <- mins %/% 60L
      mins <- mins %% 60L
      sprintf("%dh %dm", hours, mins)
    }
  }

  status <- as.character(job_row$status[1])
  rid <- as.character(job_row$rid[1])
  rtoe <- blast_parse_rtoe_seconds(job_row$rtoe[1])
  attempts <- suppressWarnings(as.integer(job_row$auto_poll_attempts[1]))
  if (!is.finite(attempts)) attempts <- 0L
  auto_on <- isTRUE(as.logical(job_row$auto_poll_enabled[1]))
  manual_req <- isTRUE(as.logical(job_row$manual_retrieval_required[1]))
  submitted <- as.character(job_row$submitted_at[1])
  checked <- as.character(job_row$last_checked_at[1])
  if (length(checked) != 1L || is.na(checked)) checked <- ""
  elapsed_secs <- as.numeric(difftime(now, blast_parse_timestamp(submitted), units = "secs"))

  lines <- c(
    paste0("Status: ", status),
    paste0("RID: ", rid),
    paste0("Elapsed: ", format_elapsed(elapsed_secs)),
    paste0("NCBI RTOE: ~", rtoe, "s")
  )
  if ("database" %in% names(job_row) && nzchar(as.character(job_row$database[1]))) {
    lines <- c(lines, paste0("DB: ", as.character(job_row$database[1])))
  }
  if ("hitlist_size" %in% names(job_row)) {
    lines <- c(lines, paste0("Hits requested: ", as.character(job_row$hitlist_size[1])))
  }
  if (nzchar(checked)) lines <- c(lines, paste0("Last NCBI check: ", checked))

  if (manual_req && status %in% c("SUBMITTED", "WAITING")) {
    lines <- c(
      lines,
      paste0("Automatic checks completed: ", attempts, " of 3"),
      "NCBI has not returned the BLAST result after 3 automatic checks.",
      "The request may still be running at NCBI.",
      "Use Check now to retrieve the result manually."
    )
    return(list(class = "status-warning", lines = lines, kind = "manual_required"))
  }

  if (auto_on && status %in% c("SUBMITTED", "WAITING")) {
    if (attempts < 1L) {
      lines <- c(lines, "Waiting for the first automatic NCBI check (after RTOE).")
    } else {
      lines <- c(
        lines,
        paste0("Automatic NCBI checks so far: ", attempts, " of 3"),
        "Waiting for the next automatic NCBI check."
      )
    }
    return(list(class = "status-warning", lines = lines, kind = "auto_active"))
  }

  if (identical(status, "READY")) {
    return(list(class = "status-ok", lines = c(lines, "Results available."), kind = "ready"))
  }
  if (grepl("FAILED|UNKNOWN", status)) {
    return(list(class = "status-warning", lines = c(lines, "NCBI reports failed or unknown job."), kind = "failed"))
  }
  if (identical(status, "STALE")) {
    return(list(class = "status-warning", lines = c(lines, "RID is stale for the active sequence revision."), kind = "stale"))
  }

  list(class = "status-warning", lines = lines, kind = "other")
}
