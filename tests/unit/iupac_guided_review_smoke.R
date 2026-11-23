# Explicit IUPAC ambiguity curation, not inferred from chromatogram dominance.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "domain", "sanger", "sequence_tools.R"))
source(file.path(root, "R", "domain", "consensus", "stage3_consensus.R"))

stopifnot(setequal(pitax_iupac_call_codes(),
                   c("A", "C", "G", "T", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")))
choices <- pitax_iupac_ambiguity_choices()
stopifnot(identical(unname(choices["M (A/C)"]), "M"))
stopifnot(identical(unname(choices["V (A/C/G)"]), "V"))
# Regression: R locales on Posit Connect may stringify Unicode punctuation
# as literal UTF-8 bytes, e.g. <e2><80><94>, in HTML select labels.
# Verify both dropdowns contain ASCII-only labels and unchanged DNA values.
all_options <- pitax_iupac_confirmation_choices()
stopifnot(identical(unname(all_options["G (guanine)"]), "G"))
stopifnot(identical(unname(all_options["N (unresolved A/C/G/T)"]), "N"))
stopifnot(setequal(unname(all_options), pitax_iupac_call_codes()))
stopifnot(all(vapply(names(choices), function(x) all(as.integer(charToRaw(x)) < 128L), logical(1))))
stopifnot(all(vapply(names(all_options), function(x) all(as.integer(charToRaw(x)) < 128L), logical(1))))
stopifnot(setequal(unname(choices), c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")))
stopifnot(identical(unname(iupac_for_pair("A", "C")), "M"))
stopifnot(identical(stage3_reverse_complement("M"), "K"))
stopifnot(identical(stage3_reverse_complement("V"), "B"))

pos <- seq.int(10L, by = 10L, length.out = 9)
trace <- matrix(0, nrow = 102L, ncol = 4)
trace[pos, 1L] <- 90
trace[pos, 2L] <- 75
trace[pos, 3L] <- 65
result <- list(
  raw_seq = "AACNNTGCA", seq = "AACNNTGCA",
  trace = trace, peak_pos = pos,
  channel_map = c(A = 1L, C = 2L, G = 3L, T = 4L),
  summary = data.frame(trim_start = 1L, trim_end = 9L,
                       trimmed_length = 9L, status = "OK")
)
result <- curation_rebuild(ensure_curation_state(result), list(min_usable_len = 1L))
original <- result$raw_seq

for (choice in c("M", "V", "R", "Y", "S", "W", "K", "B", "D", "H", "N")) {
  snap <- curation_set_base_snapshot(result, 4L, choice)
  stopifnot(is.list(snap))
  res <- curation_restore_snapshot(result, snap, list(min_usable_len = 1L))
  stopifnot(identical(substr(res$seq, 4L, 4L), choice))
  stopifnot(identical(res$raw_seq, original))
}
stopifnot(is.null(curation_set_base_snapshot(result, 4L, "Z")))
stopifnot(is.null(curation_set_base_snapshot(result, 4L, "AC")))
stopifnot(is.null(curation_set_base_snapshot(result, 4L, "")))

snap <- curation_set_base_snapshot(result, 4L, "M")
row <- data.frame(Action = "Base edit", Position = 4L, Before = "N", After = "M",
                  Method = "Manual IUPAC ambiguity", Evidence = "operator confirmed",
                  Details = "Confirmed by user", stringsAsFactors = FALSE)
result_m <- curation_commit(result, snap, row, list(min_usable_len = 1L))
stopifnot(substr(result_m$seq, 4L, 4L) == "M")
stopifnot(result_m$curation$revision == 1L)
stopifnot(result_m$curation$audit_log$Method[1] == "Manual IUPAC ambiguity")
result_v <- curation_commit(
  result_m, curation_set_base_snapshot(result_m, 4L, "V"),
  transform(row, Before = "M", After = "V"), list(min_usable_len = 1L)
)
stopifnot(substr(result_v$seq, 4L, 4L) == "V")
stopifnot(result_v$curation$revision == 2L)
stopifnot(identical(result_v$raw_seq, original))
stopifnot(identical(stage3_sanitize_dna(result_v$seq), result_v$seq))

saved <- tempfile(fileext = ".rds")
saveRDS(result_v, saved)
reloaded <- readRDS(saved)
unlink(saved)
stopifnot(identical(reloaded$seq, result_v$seq))
stopifnot(reloaded$curation$base_edits$Edited == "V")
undo <- curation_undo(result_v, list(min_usable_len = 1L))
stopifnot(isTRUE(undo$changed), substr(undo$result$seq, 4L, 4L) == "M")

# Local focus preview uses lightweight R graphics, respects curated IUPAC
# labels and is safe in headless R without a browser or Plotly widget.
plot_path <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_path, width = 6, height = 4)
plotted <- draw_chromatogram_focus_fast(result_v, 4L, window = 4L)
grDevices::dev.off()
stopifnot(isTRUE(plotted), file.exists(plot_path), file.info(plot_path)$size > 1000L)
unlink(plot_path)

# Repeated "mark reviewed" must NOT re-run expensive N peak metrics.
# A true base edit or trim change must invalidate the cached evidence.
peak_cache <- new.env(parent = emptyenv())
cached_1 <- pitax_cached_peak_flags(peak_cache, "FB001", result, scope = "raw")
stopifnot(peak_cache$compute_count == 1L)
result_reviewed <- result
result_reviewed$curation$reviewed_positions <- c(4L)
cached_2 <- pitax_cached_peak_flags(peak_cache, "FB001", result_reviewed, scope = "raw")
stopifnot(peak_cache$compute_count == 1L)
if (nrow(cached_2)) {
  stopifnot(any(cached_2$Position == 4L))
  stopifnot(all(cached_2$Review_status[cached_2$Position == 4L] == "Reviewed"))
}
changed_base <- result_reviewed
changed_base$curation$base_edits <- curation_set_base_snapshot(changed_base, 4L, "M")$base_edits
pitax_cached_peak_flags(peak_cache, "FB001", changed_base, scope = "raw")
stopifnot(peak_cache$compute_count == 2L)
changed_window <- changed_base
changed_window$summary$trim_start <- 2L
pitax_cached_peak_flags(peak_cache, "FB001", changed_window, scope = "raw")
stopifnot(peak_cache$compute_count == 3L)
pitax_cached_peak_flags(peak_cache, "DIFFERENT_READ", changed_window, scope = "raw")
stopifnot(peak_cache$compute_count == 4L)

# Production server wiring: the ambiguity dropdown must reach a dedicated
# explicit-operator event and its confirmation must accept full IUPAC codes.
server_code <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"),
                               warn = FALSE), collapse = "\n")
stopifnot(grepl('"manual_iupac_choice"', server_code, fixed = TRUE))
stopifnot(grepl('observeEvent(input$peak_review_apply_iupac', server_code, fixed = TRUE))
stopifnot(grepl('pitax_iupac_confirmation_choices()', server_code, fixed = TRUE))
stopifnot(grepl('plotOutput("peak_review_local_plot"', server_code, fixed = TRUE))
stopifnot(grepl('draw_chromatogram_focus_fast(', server_code, fixed = TRUE))
stopifnot(grepl('cancelOutput = TRUE', server_code, fixed = TRUE))

cat("iupac_guided_review_smoke: OK\n")
