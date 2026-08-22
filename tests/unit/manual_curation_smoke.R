# Smoke tests for v2.11.0 manual curation, bulk-correction proposal and undo/redo.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "domain", "sanger", "sequence_tools.R"))

trace <- matrix(0, nrow = 70, ncol = 4)
peak_pos <- c(10L, 20L, 30L, 40L, 50L, 60L)
channel_map <- c(A=1L, C=2L, G=3L, T=4L)
trace[10,] <- c(100,10,5,5)
trace[20,] <- c(5,100,10,5)
# Position 3: called G=40, T=100, third channel=10. Alternative is centered exactly at the call.
trace[30,] <- c(5,10,40,100)
trace[40,] <- c(5,5,100,10)
trace[50,] <- c(10,5,5,100)
trace[60,] <- c(100,5,5,10)

result <- list(
  sample_id = "TEST",
  raw_seq = "ACGTTA",
  seq = "ACGTTA",
  trace = trace,
  peak_pos = peak_pos,
  channel_map = channel_map,
  metrics = data.frame(called_signal=c(100,100,40,100,100,100)),
  summary = data.frame(trim_start=1L, trim_end=6L, trimmed_length=6L, status="OK", stringsAsFactors=FALSE)
)
settings <- list(min_usable_len=1L)
result <- ensure_curation_state(result)
result <- curation_rebuild(result, settings)

proposals <- high_confidence_autocorrections(result, settings)
idx3 <- which(proposals$Position == 3L & proposals$Competing_channel == "T")
if (!length(idx3)) {
  all_flags <- ambiguous_peak_flags(result, scope = "trimmed", params = ambiguous_peak_params_from_settings(settings))
  print(all_flags)
  stop("Expected a high-confidence G->T proposal at raw base position 3, but it was not returned.")
}
stopifnot(isTRUE(proposals$Auto_correct_candidate[idx3[1]]))

# Auto-correction proposal thresholds are configurable through settings.
strict_settings <- settings
strict_settings$auto_correct_min_alt_to_called <- 3.00
strict_proposals <- high_confidence_autocorrections(result, strict_settings)
stopifnot(!any(strict_proposals$Position == 3L, na.rm = TRUE))
relaxed_settings <- settings
relaxed_settings$auto_correct_min_alt_to_called <- 1.20
relaxed_settings$auto_correct_min_alt_to_third <- 1.20
relaxed_settings$auto_correct_max_peak_offset <- 3L
relaxed_settings$auto_correct_min_relative_signal <- 0.20
relaxed_proposals <- high_confidence_autocorrections(result, relaxed_settings)
stopifnot(any(relaxed_proposals$Position == 3L, na.rm = TRUE))

snap <- curation_set_base_snapshot(result, 3L, "T")
row <- data.frame(Action="Base edit", Position=3L, Before="G", After="T", Method="Manual base correction", Evidence="test", Details="test", stringsAsFactors=FALSE)
edited <- curation_commit(result, snap, row, settings, "Test edit")
stopifnot(edited$seq == "ACTTTA")
stopifnot(nrow(edited$curation$base_edits) == 1L)
stopifnot(nrow(edited$curation$audit_log) == 1L)

undo <- curation_undo(edited, settings)
stopifnot(isTRUE(undo$changed))
stopifnot(undo$result$seq == "ACGTTA")
redo <- curation_redo(undo$result, settings)
stopifnot(isTRUE(redo$changed))
stopifnot(redo$result$seq == "ACTTTA")

trim_snap <- curation_trim_snapshot(redo$result, 5L, "right")
trim_row <- data.frame(Action="Manual trim right", Position=5L, Before="1-6", After="1-4", Method="Manual chromatogram curation", Evidence="test", Details="test", stringsAsFactors=FALSE)
trimmed <- curation_commit(redo$result, trim_snap, trim_row, settings, "Test trim")
stopifnot(trimmed$seq == "ACTT")
stopifnot(trimmed$summary$trim_end == 4L)

cat("v2.14.2 manual curation smoke tests passed.\n")
