# Stored AB1 re-trim regression. No Shiny session and no AB1 temp files.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "domain", "sanger", "core_sanger.R"))
source(file.path(app_dir, "R", "domain", "sanger", "sequence_tools.R"))
source(file.path(app_dir, "R", "domain", "sanger", "retrim_saved.R"))

make_saved_read <- function(id) {
  raw <- paste0(strrep("N", 8), strrep("A", 12))
  pos <- seq.int(10L, by = 10L, length.out = nchar(raw))
  trace <- matrix(0, nrow = max(pos) + 10L, ncol = 4)
  trace[pos[9:20], 1] <- 100
  trace[pos[9:20], 2] <- 5
  settings <- list(
    window = 10L, min_peak_ratio = 1.5,
    min_len_before_collapse = 15L, min_relative_signal = 0.2,
    bad_run_windows = 3L, absolute_max_base_index = 50L,
    min_usable_len = 1L, target = "TEF1", forward_primer = "F",
    reverse_primer = "R", sequencing_primer = "Forward"
  )
  result <- list(
    sample_id = id, raw_seq = raw, seq = raw,
    trace = trace, peak_pos = pos,
    channel_map = c(A = 1L, C = 2L, G = 3L, T = 4L),
    ab1_evidence = list(detail = data.frame(Basecaller_quality = rep(25, 20))),
    processing_settings = settings,
    read_assignment = list(Source_ID = id, Isolate = id, Assay_ID = "assay-tef1",
                           Locus = "TEF1", Direction = "Forward",
                           Final_Name = paste0(id, "_TEF1_F")),
    summary = data.frame(
      sample_id = id, target = "TEF1",
      raw_length = 20L, trimmed_length = 20L,
      trim_start = 1L, trim_end = 20L,
      collapse_index = NA_integer_, reason = "ABSOLUTE_CAP",
      mean_called_signal_trimmed = 100,
      median_peak_ratio_trimmed = 20,
      status = "OK", error_message = "", stringsAsFactors = FALSE
    )
  )
  ensure_curation_state(result)
}

one <- make_saved_read("FB001")
# Preserve prior audited manual base edits and reflect them in the new trim.
snap <- curation_set_base_snapshot(one, 12L, "G")
edit <- data.frame(
  Action = "Base edit", Position = 12L, Before = "A", After = "G",
  Method = "Manual", Evidence = "trace", Details = "approved",
  stringsAsFactors = FALSE
)
one <- curation_commit(one, snap, edit, one$processing_settings, "manual edit")
old_curated <- one$seq
stopifnot(nchar(old_curated) == 20L, substring(old_curated, 12L, 12L) == "G")
out <- pitax_retrim_saved_read(one)
stopifnot(out$curation$auto_trim_start == 7L)
stopifnot(out$curation$auto_trim_end == 20L)
stopifnot(out$summary$trim_start == 7L, out$summary$trim_end == 20L)
stopifnot(nchar(out$seq) == 14L)
stopifnot(substring(out$seq, 6L, 6L) == "G")
stopifnot(identical(out$read_assignment, one$read_assignment))
stopifnot(identical(out$processing_settings, one$processing_settings))
stopifnot(identical(out$raw_seq, one$raw_seq))
stopifnot(identical(out$ab1_evidence, one$ab1_evidence))
stopifnot(identical(out$trace, one$trace))
stopifnot(nrow(out$curation$base_edits) == 1L)
stopifnot(nrow(out$curation$audit_log) == 2L)
stopifnot(length(out$curation$undo_stack) == 0L)
stopifnot(length(out$retrim_history) == 1L)
stopifnot(out$curation$revision == one$curation$revision + 1L)
stopifnot(nchar(one$seq) == 20L) # original in-memory input not mutated

# Repeated re-trim does not accumulate revisions or destroy edits.
again <- pitax_retrim_saved_read(out)
stopifnot(identical(out$seq, again$seq))
stopifnot(identical(out$curation$revision, again$curation$revision))
stopifnot(length(again$retrim_history) == 1L)

# Whole project (50 reads) inherits all IDs and source read names without
# touching upload or assignment editor state.
ids <- sprintf("FB%03d", seq_len(50))
saved <- stats::setNames(lapply(ids, make_saved_read), ids)
batch <- pitax_retrim_saved_batch(saved)
stopifnot(nrow(batch$summary) == 50L)
stopifnot(length(batch$changed_sources) == 50L)
stopifnot(identical(names(batch$results), ids))
stopifnot(all(vapply(ids, function(id) identical(
  saved[[id]]$read_assignment, batch$results[[id]]$read_assignment
), logical(1))))
stopifnot(all(vapply(batch$results, function(x) x$summary$trim_start[1] == 7L, logical(1))))

# One missing chromatogram aborts the entire preview (no partial writes).
broken <- saved
broken[["FB025"]]$trace <- NULL
error <- tryCatch(pitax_retrim_saved_batch(broken), error = function(e) conditionMessage(e))
stopifnot(is.character(error), grepl("FB025", error, fixed = TRUE))
stopifnot(identical(saved[["FB001"]]$summary$trim_start, 1L))

# If the stored AB1 cannot be fully reconstructed, explain re-upload.
missing <- make_saved_read("legacy")
missing$processing_settings <- NULL
stopifnot(grepl("re-upload", pitax_saved_retrim_requirements(missing), fixed = TRUE))

cat("saved_project_retrim_smoke: OK\n")
