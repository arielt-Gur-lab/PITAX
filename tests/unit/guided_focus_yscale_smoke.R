# Display-only Focus / Full Y-scale for guided local chromatograms.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "domain", "sanger", "sequence_tools.R"))
values <- matrix(c(
  100, 80, 60, 40, 700, 520, 180, 90, 820, 610, 200, 110,
  2800, 140, 90, 70, 3000, 140, 90, 70, 2900, 140, 90, 70,
  760, 580, 210, 100, 690, 500, 170, 80, 120, 90, 50, 30
), ncol = 4L, byrow = TRUE)
focus <- pitax_chromatogram_focus_yscale(values, mode = "focus")
full <- pitax_chromatogram_focus_yscale(values, mode = "full")
stopifnot(isTRUE(abs(full$display_max - 3000) < 1e-9))
stopifnot(focus$display_max < full$display_max * 0.95)
stopifnot(focus$display_max >= 820, isTRUE(focus$clipped))
uniform <- matrix(rep(c(500, 480, 460, 440), each = 8L), ncol = 4L)
stopifnot(isFALSE(pitax_chromatogram_focus_yscale(uniform, "focus")$clipped))
n_bases <- 9L; peak_pos <- seq.int(10L, by = 10L, length.out = n_bases)
trace <- matrix(120, nrow = 110L, ncol = 4L)
for (i in seq_len(n_bases)) trace[peak_pos[i], ] <- c(800, 600, 200, 100)
trace[peak_pos[5L], 1L] <- 3200
result <- list(raw_seq="ACGTNACGT", seq="ACGTNACGT", trace=trace, peak_pos=peak_pos,
  channel_map=c(A=1L,C=2L,G=3L,T=4L),
  summary=data.frame(trim_start=1L,trim_end=n_bases,trimmed_length=n_bases,status="OK"))
result <- ensure_curation_state(result); before <- result$trace
pdf(tempfile(fileext=".pdf")); stopifnot(draw_chromatogram_focus_fast(result,5L,4L,"focus")); stopifnot(draw_chromatogram_focus_fast(result,5L,4L,"full")); dev.off()
stopifnot(identical(result$trace, before))
server <- paste(readLines(file.path(root,"R/server/stages/50_evidence_review.R"), warn=FALSE), collapse="\n")
stopifnot(grepl('"peak_review_yscale_focus"', server, fixed=TRUE))
stopifnot(grepl('"peak_review_yscale_full"', server, fixed=TRUE))
stopifnot(grepl("rv$peak_review_yscale", server, fixed=TRUE))
stopifnot(grepl('height = "340px"', server, fixed=TRUE))
state <- paste(readLines(file.path(root,"R/server/stages/00_state.R"), warn=FALSE), collapse="\n")
stopifnot(grepl('peak_review_yscale = "focus"', state, fixed=TRUE))
cat("guided_focus_yscale_smoke: OK\n")