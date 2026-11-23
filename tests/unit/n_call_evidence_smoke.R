# N-call evidence is manual-review assistance, never a second basecaller.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "domain", "sanger", "sequence_tools.R"))

trace <- matrix(0, nrow = 100, ncol = 4)
peak_pos <- c(10L, 25L, 40L, 55L, 70L)
# A well-separated C-shaped summit at the unresolved second base.
trace[18:32, 2] <- 120 * exp(-0.5 * ((18:32 - 25) / 2)^2)
trace[, 1] <- 4
trace[25, 4] <- 12
result <- list(
  raw_seq = "ANAAA", seq = "ANAAA",
  trace = trace, peak_pos = peak_pos,
  channel_map = c(A = 1L, C = 2L, G = 3L, T = 4L),
  summary = data.frame(trim_start = 1L, trim_end = 5L),
  metrics = data.frame(called_signal = rep(50, 5)),
  ab1_evidence = list(detail = data.frame(Basecaller_quality = c(33, 8, 32, 30, 31)))
)

e <- pitax_n_call_evidence(result, 2L)
stopifnot(identical(e$Top_channel, "C"))
stopifnot(identical(e$Basecaller_quality, 8))
stopifnot(identical(e$Peak_separation, "Distinct local summit"))
stopifnot(is.finite(e$Top_to_second_ratio) && e$Top_to_second_ratio > 5)
stopifnot(e$Peak_prominence_percent >= 20)
stopifnot(grepl("low basecaller Q", e$Review_hint, fixed = TRUE))

# N is still N: neither evidence generation nor peak flags modify the read.
flags <- ambiguous_peak_flags(result, scope = "raw")
n_flag <- flags[flags$Position == 2L, , drop = FALSE]
stopifnot(nrow(n_flag) == 1L, identical(n_flag$Top_channel, "C"))
stopifnot(identical(result$raw_seq, "ANAAA"))

# A C plateau over multiple basecall intervals cannot be interpreted as
# multiple resolved C bases simply because its dye dominates.
flat <- result
flat$trace <- matrix(0, nrow = 100, ncol = 4)
flat$trace[, 1] <- 5
flat$trace[17:33, 2] <- 120
plateau <- pitax_n_call_evidence(flat, 2L)
stopifnot(identical(plateau$Top_channel, "C"))
stopifnot(identical(plateau$Peak_separation, "Broad/overlapping or unseparated"))
stopifnot(grepl("base count is unresolved", plateau$Review_hint, fixed = TRUE))

# Graceful quality fallback for legacy projects; no Q invented.
legacy <- result
legacy$ab1_evidence <- NULL
missing_q <- pitax_n_call_evidence(legacy, 2L)
stopifnot(is.na(missing_q$Basecaller_quality))
stopifnot(identical(missing_q$Top_channel, "C"))

cat("n_call_evidence_smoke: OK\n")
