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
stopifnot(identical(unname(choices["M — A/C"]), "M"))
stopifnot(identical(unname(choices["V — A/C/G"]), "V"))
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

cat("iupac_guided_review_smoke: OK\n")
