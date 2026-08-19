# PITAX v3.0.0-alpha.1.1 - Stage 1 AB1 evidence helper tests.
# These are synthetic tests; the Stage 1 manual checklist validates real AB1 files.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "ab1_evidence.R"))

# 1. Canonical A/C/G/T mapping is explicit and does not need permutation inference.
x <- matrix(0, nrow = 4, ncol = 4)
colnames(x) <- c("A", "C", "G", "T")
map <- pitax_canonical_channel_map(x)
stopifnot(identical(unname(as.integer(map[c("A","C","G","T")])), 1:4))

# 2. A called base selects its own column from peakPosMatrix.
seq_string <- "ACGTN"
pm <- matrix(
  c(
    10,11,12,13,
    20,21,22,23,
    30,31,32,33,
    40,41,42,43,
    50,51,52,53
  ),
  nrow = 5, byrow = TRUE
)
colnames(pm) <- c("A","C","G","T")
pos <- pitax_called_peak_positions(seq_string, pm)
stopifnot(identical(as.numeric(pos[1:4]), c(10,21,32,43)))
stopifnot(is.na(pos[5]))

# 3. Signal evidence is evaluated at the supplied called-base-specific position.
trace <- matrix(1, nrow = 60, ncol = 4)
colnames(trace) <- c("A","C","G","T")
trace[10,] <- c(100,10,8,5)
trace[21,] <- c(8,110,9,7)
trace[32,] <- c(6,8,120,9)
trace[43,] <- c(7,6,8,130)
met <- pitax_signal_metrics(seq_string, trace, pos, map)
stopifnot(all(met$called_is_max[1:4] %in% TRUE))
stopifnot(identical(met$best_channel[1:4], c("A","C","G","T")))
stopifnot(all(met$called_to_alt_ratio[1:4] >= 10))

# 4. The legacy first-column peak positions are demonstrably different for C/G/T.
legacy_pos <- as.numeric(pm[,1])
comparable <- is.finite(legacy_pos) & is.finite(pos)
stopifnot(sum(legacy_pos[comparable] != pos[comparable]) == 3L)

# 5. peakAmpMatrix is interpreted in the same documented A/C/G/T column model.
pam <- matrix(
  c(
    100,10,8,5,
    8,110,9,7,
    6,8,120,9,
    7,6,8,130,
    1,1,1,1
  ),
  nrow = 5, byrow = TRUE
)
colnames(pam) <- c("A","C","G","T")
amp <- pitax_peak_amplitude_metrics(seq_string, pam)
stopifnot(all(amp$peakamp_called_is_max[1:4] %in% TRUE))
stopifnot(identical(amp$peakamp_best_channel[1:4], c("A","C","G","T")))

# 6. Quality vectors are padded/truncated without shifting the base index.
stopifnot(identical(pitax_normalize_numeric(c(10,20), 4), c(10,20,NA_real_,NA_real_)))
stopifnot(identical(pitax_normalize_numeric(c(10,20,30), 2), c(10,20)))

cat("v3.0.0-alpha.1.1 AB1 evidence helper tests passed.\n")
