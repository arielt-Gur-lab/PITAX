# BLAST batch status stays compact even when many RID jobs are WAITING.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "services", "blast_polling.R"))

message <- blast_batch_summary_text(
  "Batch retrieval", ready = 10, waiting = 8,
  too_soon = 0, no_hits = 0, stale = 0, failed = 0,
  elapsed = "1m 10s"
)
stopifnot(grepl("10 ready", message, fixed = TRUE))
stopifnot(grepl("8 still running at NCBI", message, fixed = TRUE))
stopifnot(grepl("1m 10s", message, fixed = TRUE))
stopifnot(!grepl("RID", message, fixed = TRUE))
stopifnot(!grepl("too soon", message, fixed = TRUE))
stopifnot(nchar(message) < 120L)
with_errors <- blast_batch_summary_text("Selected retrieval", ready = 2,
                                       waiting = 3, too_soon = 1,
                                       stale = 1, failed = 1)
stopifnot(grepl("1 stale", with_errors, fixed = TRUE))
stopifnot(grepl("1 failed", with_errors, fixed = TRUE))
stopifnot(grepl("1 too soon", with_errors, fixed = TRUE))
cat("blast_batch_summary_smoke: OK\n")
