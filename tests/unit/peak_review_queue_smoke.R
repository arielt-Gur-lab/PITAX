# Smoke tests for guided ambiguous-peak review queue helpers.

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

df <- data.frame(
  Position = c(10L, 25L, 40L, 55L),
  Review_status = c("Active", "Reviewed", "Active", "Active"),
  stringsAsFactors = FALSE
)

stopifnot(identical(pitax_next_active_peak_position(df, NA_integer_), 10L))
stopifnot(identical(pitax_next_active_peak_position(df, 10L), 40L))
stopifnot(identical(pitax_next_active_peak_position(df, 40L), 55L))
stopifnot(identical(pitax_next_active_peak_position(df, 55L), NA_integer_))
stopifnot(identical(pitax_next_active_peak_position(df, 10L, deferred = 40L), 55L))
stopifnot(identical(pitax_next_active_peak_position(df, NA_integer_, deferred = c(10L, 40L)), 55L))
stopifnot(identical(pitax_next_active_peak_position(df, NA_integer_, deferred = c(10L, 40L, 55L)), NA_integer_))
stopifnot(identical(pitax_next_active_peak_position(df[0, , drop = FALSE], 10L), NA_integer_))

cat("peak_review_queue_smoke: OK\n")
