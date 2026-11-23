# Regression: unresolved AB1 base calls must count as bad positions when
# evaluating the minimum good-base fraction at an automatic trim start.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "domain", "sanger", "core_sanger.R"))

make_metrics <- function(ratios, signals = ifelse(is.na(ratios), NA_real_, 100)) {
  data.frame(peak_ratio = as.numeric(ratios), called_signal = as.numeric(signals))
}

# Previously this returned 1: mean(c(NA, ..., TRUE, TRUE), na.rm=TRUE)
# incorrectly treated only two usable calls in a 10-base window as 100% good.
n_prefix <- make_metrics(c(rep(NA_real_, 8), rep(2, 12)))
stopifnot(identical(find_start(n_prefix, window = 10L, min_peak_ratio = 1.5), 7L))

# Seven resolved calls and three N calls are 70%, NOT an 80% pass.
below_cutoff <- make_metrics(c(rep(NA_real_, 3), rep(2, 7)))
stopifnot(is.na(find_start(below_cutoff, window = 10L, min_peak_ratio = 1.5)))

# Exactly 80% good calls are acceptable, and a fully clean window starts at 1.
at_cutoff <- make_metrics(c(rep(NA_real_, 2), rep(2, 8)))
stopifnot(identical(find_start(at_cutoff, window = 10L, min_peak_ratio = 1.5), 1L))
all_good <- make_metrics(rep(2, 10))
stopifnot(identical(find_start(all_good, window = 10L, min_peak_ratio = 1.5), 1L))

# A window with only unresolved positions must still be rejected.
all_n <- make_metrics(rep(NA_real_, 10))
stopifnot(is.na(find_start(all_n, window = 10L, min_peak_ratio = 1.5)))

# A genuine signal without a detectable alternative may have an infinite
# peak ratio; this remains a valid quality pass when the called signal is >0.
unopposed <- make_metrics(rep(Inf, 10), rep(100, 10))
stopifnot(identical(find_start(unopposed, window = 10L, min_peak_ratio = 1.5), 1L))

cat("qc_n_trim_window_smoke: OK\n")
