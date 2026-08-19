# Smoke tests for v2.11.0 ambiguous-peak / channel-competition review.

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
app_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "sequence_tools.R"))

trace <- matrix(0, nrow = 60, ncol = 4)
peak_pos <- c(10L, 20L, 30L, 40L, 50L)
channel_map <- c(A=1L, C=2L, G=3L, T=4L)

# Base 1: clear A, should not flag.
trace[10,] <- c(100,20,10,5)
# Base 2: C/G nearly tied -> Strong competition.
trace[20,] <- c(5,100,90,5)
# Base 3: G/T moderately close -> Moderate competition.
trace[30,] <- c(5,5,100,65)
# Base 4: called T but A is stronger -> Called base not dominant.
trace[40,] <- c(100,5,5,80)
# Base 5: N call with two strong channels -> always reviewable.
trace[50,] <- c(5,90,100,5)

result <- list(
  raw_seq = "ACGTN",
  seq = "ACGTN",
  trace = trace,
  peak_pos = peak_pos,
  channel_map = channel_map,
  metrics = data.frame(called_signal=c(100,100,100,80,NA_real_)),
  summary = data.frame(trim_start=1L, trim_end=5L)
)

flags <- ambiguous_peak_flags(result, scope="trimmed")
stopifnot(identical(flags$Position, c(2L,3L,4L,5L)))
stopifnot(flags$Flag[flags$Position==2] == "Strong channel competition")
stopifnot(flags$Flag[flags$Position==3] == "Moderate channel competition")
stopifnot(flags$Flag[flags$Position==4] == "Called base not dominant")
stopifnot(flags$Flag[flags$Position==5] == "Ambiguous N base call")

cat("v2.11.0 QC peak-flag smoke tests passed.\n")
