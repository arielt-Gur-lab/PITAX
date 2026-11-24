# Display-only Focus / Full Y-scale for guided local chromatograms.
# Outlier spikes must not flatten competing peaks; stored traces stay intact.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "domain", "sanger", "sequence_tools.R"))

# Synthetic local window: typical competing peaks ~700–900, one spike at 3000.
values <- matrix(c(
  100, 80, 60, 40,
  700, 520, 180, 90,
  820, 610, 200, 110,
  3000, 140, 90, 70,
  760, 580, 210, 100,
  690, 500, 170, 80,
  120, 90, 50, 30
), ncol = 4L, byrow = TRUE)

focus <- pitax_chromatogram_focus_yscale(values, mode = "focus")
full <- pitax_chromatogram_focus_yscale(values, mode = "full")
stopifnot(identical(full$mode, "full"))
stopifnot(isTRUE(abs(full$full_max - 3000) < 1e-9))
stopifnot(isTRUE(abs(full$display_max - 3000) < 1e-9))
stopifnot(isFALSE(full$clipped))
stopifnot(!length(full$clip_rows))

stopifnot(identical(focus$mode, "focus"))
stopifnot(isTRUE(abs(focus$full_max - 3000) < 1e-9))
stopifnot(focus$display_max < 2000)
stopifnot(focus$display_max > 600)
stopifnot(isTRUE(focus$clipped))
stopifnot(4L %in% focus$clip_rows)
# Competing peaks around 700–820 must remain inside the Focus ceiling.
stopifnot(focus$display_max >= 820)

# Uniform window: Focus and Full should agree (no artificial clipping).
uniform <- matrix(rep(c(500, 480, 460, 440), each = 8L), ncol = 4L)
u_focus <- pitax_chromatogram_focus_yscale(uniform, mode = "focus")
u_full <- pitax_chromatogram_focus_yscale(uniform, mode = "full")
stopifnot(isFALSE(u_focus$clipped))
stopifnot(isTRUE(abs(u_focus$display_max - u_full$display_max) < 1e-9))

# Drawing must not mutate the stored chromatogram matrix.
n_bases <- 9L
peak_pos <- seq.int(10L, by = 10L, length.out = n_bases)
trace <- matrix(120, nrow = 110L, ncol = 4L)
for (i in seq_len(n_bases)) {
  p <- peak_pos[i]
  trace[p, ] <- c(800, 600, 200, 100)
}
trace[peak_pos[5L], 1L] <- 3200
seq_calls <- "ACGTNACGT"
result <- list(
  raw_seq = seq_calls, seq = seq_calls,
  trace = trace, peak_pos = peak_pos,
  channel_map = c(A = 1L, C = 2L, G = 3L, T = 4L),
  summary = data.frame(trim_start = 1L, trim_end = n_bases,
                       trimmed_length = n_bases, status = "OK")
)
result <- ensure_curation_state(result)
before <- result$trace
plot_path <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_path, width = 6, height = 3.4)
ok_focus <- draw_chromatogram_focus_fast(result, 5L, window = 4L, yscale = "focus")
ok_full <- draw_chromatogram_focus_fast(result, 5L, window = 4L, yscale = "full")
grDevices::dev.off()
stopifnot(isTRUE(ok_focus), isTRUE(ok_full))
stopifnot(identical(result$trace, before))
stopifnot(identical(result$raw_seq, seq_calls))
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 1000L)
unlink(plot_path)

# UI contract: toggle exists; plot height raised; legend/trim remain outside.
server <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"),
                          warn = FALSE), collapse = "\n")
stopifnot(grepl('"peak_review_yscale"', server, fixed = TRUE))
stopifnot(grepl('"Focus Scale" = "focus"', server, fixed = TRUE))
stopifnot(grepl('"Full Scale" = "full"', server, fixed = TRUE))
stopifnot(grepl('height = "260px"', server, fixed = TRUE))
stopifnot(grepl("yscale = yscale", server, fixed = TRUE))
stopifnot(grepl('class = "peak-review-dye-legend"', server, fixed = TRUE))
stopifnot(grepl('class = "peak-review-trim-inline"', server, fixed = TRUE))
css <- paste(readLines(file.path(root, "www", "pitax.css"), warn = FALSE), collapse = "\n")
stopifnot(grepl(".peak-review-yscale-toggle", css, fixed = TRUE))

cat("guided_focus_yscale_smoke: OK\n")
