# ============================================================
# Sanger Sequence Pipeline - Shiny App
#
# Pipeline:
#   1. Upload AB1 files
#   2. Define target/gene, primers and trimming parameters
#   3. Trim + QC
#   4. Optional sample renaming
#   5. Export FASTA / BLAST FASTA / CSV / ZIP
#
# Uses AB1:
#   primarySeq
#   traceMatrix
#   peakPosMatrix[,1]
#
# Does NOT use PHD.1
# ============================================================


# ============================================================
# Package installation
# ============================================================

required_cran <- c(
  "shiny",
  "DT",
  "zip",
  "readxl"
)

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("sangerseqR", quietly = TRUE)) {
  BiocManager::install("sangerseqR")
}


# ============================================================
# Libraries
# ============================================================

library(shiny)
library(DT)
library(sangerseqR)


# Allow relatively large AB1 batches
options(
  shiny.maxRequestSize = 500 * 1024^2
)


# ============================================================
# Helper functions
# ============================================================

wrap_sequence <- function(seq, width = 80) {
  
  if (is.null(seq) || nchar(seq) == 0) {
    return("")
  }
  
  starts <- seq(
    1,
    nchar(seq),
    by = width
  )
  
  paste(
    substring(
      seq,
      starts,
      pmin(
        starts + width - 1,
        nchar(seq)
      )
    ),
    collapse = "\n"
  )
}


rolling_median <- function(x, window) {
  
  n <- length(x)
  
  out <- rep(
    NA_real_,
    n
  )
  
  if (n < window) {
    return(out)
  }
  
  for (i in seq_len(n - window + 1)) {
    
    out[i] <- median(
      x[i:(i + window - 1)],
      na.rm = TRUE
    )
  }
  
  out
}


make_permutations <- function(x) {
  
  if (length(x) == 1) {
    return(
      matrix(
        x,
        nrow = 1
      )
    )
  }
  
  out <- NULL
  
  for (i in seq_along(x)) {
    
    rest <- x[-i]
    
    sub <- make_permutations(
      rest
    )
    
    out <- rbind(
      out,
      cbind(
        x[i],
        sub
      )
    )
  }
  
  out
}


# ============================================================
# Infer A/C/G/T channel order
# ============================================================

infer_channel_map <- function(
    seq_string,
    trace,
    peak_pos,
    max_index
) {
  
  bases <- strsplit(
    seq_string,
    ""
  )[[1]]
  
  bases4 <- c(
    "A",
    "C",
    "G",
    "T"
  )
  
  perms <- make_permutations(
    1:4
  )
  
  n_check <- min(
    length(bases),
    length(peak_pos),
    max_index
  )
  
  idx <- which(
    bases[1:n_check] %in% bases4 &
      !is.na(
        peak_pos[1:n_check]
      )
  )
  
  if (length(idx) < 50) {
    
    stop(
      "Too few usable bases to infer trace channel order."
    )
  }
  
  scores <- data.frame(
    
    perm_id =
      seq_len(
        nrow(perms)
      ),
    
    pct_called_is_max =
      NA_real_,
    
    median_called_signal =
      NA_real_
  )
  
  
  for (p in seq_len(nrow(perms))) {
    
    map <- setNames(
      perms[p, ],
      bases4
    )
    
    called_is_max <- rep(
      NA,
      length(idx)
    )
    
    called_signal <- rep(
      NA_real_,
      length(idx)
    )
    
    
    for (j in seq_along(idx)) {
      
      i <- idx[j]
      
      b <- bases[i]
      
      pos <- peak_pos[i]
      
      
      if (
        is.na(pos) ||
        pos < 1 ||
        pos > nrow(trace)
      ) {
        next
      }
      
      
      vals <- trace[
        pos,
      ]
      
      
      called <- vals[
        map[[b]]
      ]
      
      
      called_signal[j] <- called
      
      
      called_is_max[j] <-
        called ==
        max(
          vals,
          na.rm = TRUE
        )
    }
    
    
    scores$pct_called_is_max[p] <-
      mean(
        called_is_max,
        na.rm = TRUE
      )
    
    
    scores$median_called_signal[p] <-
      median(
        called_signal,
        na.rm = TRUE
      )
  }
  
  
  best <- scores[
    order(
      -scores$pct_called_is_max,
      -scores$median_called_signal
    ),
  ][1, ]
  
  
  best_map <- setNames(
    perms[
      best$perm_id,
    ],
    bases4
  )
  
  
  list(
    map = best_map,
    score = best
  )
}


# ============================================================
# Calculate AB1 quality metrics
# ============================================================

calculate_metrics <- function(
    seq_string,
    trace,
    peak_pos,
    channel_map
) {
  
  bases <- strsplit(
    seq_string,
    ""
  )[[1]]
  
  
  n <- min(
    length(bases),
    length(peak_pos)
  )
  
  
  bases <- bases[
    seq_len(n)
  ]
  
  
  peak_pos <- peak_pos[
    seq_len(n)
  ]
  
  
  called_signal <- rep(
    NA_real_,
    n
  )
  
  
  second_signal <- rep(
    NA_real_,
    n
  )
  
  
  peak_ratio <- rep(
    NA_real_,
    n
  )
  
  
  for (i in seq_len(n)) {
    
    b <- bases[i]
    
    pos <- peak_pos[i]
    
    
    if (
      !(b %in% names(channel_map))
    ) {
      next
    }
    
    
    if (
      is.na(pos) ||
      pos < 1 ||
      pos > nrow(trace)
    ) {
      next
    }
    
    
    vals <- trace[
      pos,
    ]
    
    
    called_col <-
      channel_map[[b]]
    
    
    called_signal[i] <-
      vals[
        called_col
      ]
    
    
    ordered_vals <- sort(
      vals,
      decreasing = TRUE,
      na.last = NA
    )
    
    
    if (length(ordered_vals) >= 2) {
      
      second_signal[i] <-
        ordered_vals[2]
      
    } else {
      
      second_signal[i] <-
        NA_real_
    }
    
    
    if (
      !is.na(second_signal[i]) &&
      second_signal[i] > 0
    ) {
      
      peak_ratio[i] <-
        called_signal[i] /
        second_signal[i]
      
    } else if (
      !is.na(called_signal[i]) &&
      called_signal[i] > 0
    ) {
      
      peak_ratio[i] <- Inf
    }
  }
  
  
  data.frame(
    
    index =
      seq_len(n),
    
    base =
      bases,
    
    peak_pos =
      peak_pos,
    
    called_signal =
      called_signal,
    
    second_signal =
      second_signal,
    
    peak_ratio =
      peak_ratio,
    
    stringsAsFactors =
      FALSE
  )
}


# ============================================================
# Find good sequence start
# ============================================================

find_start <- function(
    metrics,
    window,
    min_peak_ratio
) {
  
  good <-
    metrics$peak_ratio >= min_peak_ratio &
    metrics$called_signal > 0
  
  
  n <- length(good)
  
  
  if (n < window) {
    return(
      NA_integer_
    )
  }
  
  
  for (
    i in seq_len(
      n - window + 1
    )
  ) {
    
    idx <-
      i:(i + window - 1)
    
    
    good_fraction <- mean(
      good[idx],
      na.rm = TRUE
    )
    
    
    if (
      !is.na(good_fraction) &&
      good_fraction >= 0.80
    ) {
      
      return(i)
    }
  }
  
  
  NA_integer_
}


# ============================================================
# Find signal collapse
# ============================================================

find_collapse <- function(
    metrics,
    start_index,
    window,
    min_len_before_collapse,
    min_relative_signal,
    min_peak_ratio,
    bad_run_windows,
    absolute_max_base_index
) {
  
  if (is.na(start_index)) {
    
    return(
      NA_integer_
    )
  }
  
  
  n <- nrow(metrics)
  
  
  search_from <-
    start_index +
    min_len_before_collapse -
    1
  
  
  search_to <- min(
    n - window + 1,
    absolute_max_base_index
  )
  
  
  if (
    search_from >= search_to
  ) {
    
    return(
      NA_integer_
    )
  }
  
  
  baseline_end <- min(
    start_index + 250,
    n,
    absolute_max_base_index
  )
  
  
  baseline_signal <- median(
    metrics$called_signal[
      start_index:baseline_end
    ],
    na.rm = TRUE
  )
  
  
  if (
    is.na(baseline_signal) ||
    baseline_signal <= 0
  ) {
    
    return(
      NA_integer_
    )
  }
  
  
  roll_signal <- rolling_median(
    metrics$called_signal,
    window
  )
  
  
  roll_ratio <- rolling_median(
    metrics$peak_ratio,
    window
  )
  
  
  relative_signal <-
    roll_signal /
    baseline_signal
  
  
  bad <-
    (
      !is.na(relative_signal) &
        relative_signal <
        min_relative_signal
    ) |
    (
      !is.na(roll_ratio) &
        roll_ratio <
        min_peak_ratio
    )
  
  
  last_start <-
    search_to -
    bad_run_windows +
    1
  
  
  if (
    last_start < search_from
  ) {
    
    return(
      NA_integer_
    )
  }
  
  
  for (
    i in seq(
      search_from,
      last_start
    )
  ) {
    
    idx <-
      i:(i + bad_run_windows - 1)
    
    
    this_bad <- bad[idx]
    
    
    if (
      length(this_bad) ==
      bad_run_windows &&
      all(
        this_bad,
        na.rm = FALSE
      )
    ) {
      
      return(i)
    }
  }
  
  
  NA_integer_
}


# ============================================================
# Trim one AB1
# ============================================================

trim_one_ab1 <- function(
    file,
    sample_id,
    settings
) {
  
  s <- readsangerseq(
    file
  )
  
  
  seq_string <- toupper(
    as.character(
      s@primarySeq
    )
  )
  
  
  seq_string <- gsub(
    "[^ACGTN]",
    "N",
    seq_string
  )
  
  
  trace <- as.matrix(
    s@traceMatrix
  )
  
  
  peak_pos <- as.integer(
    s@peakPosMatrix[
      ,
      1
    ]
  )
  
  
  if (
    length(peak_pos) !=
    nchar(seq_string)
  ) {
    
    stop(
      sample_id,
      ": primarySeq length and peakPosMatrix rows do not match."
    )
  }
  
  
  inferred <- infer_channel_map(
    
    seq_string =
      seq_string,
    
    trace =
      trace,
    
    peak_pos =
      peak_pos,
    
    max_index =
      settings$absolute_max_base_index
  )
  
  
  metrics <- calculate_metrics(
    
    seq_string =
      seq_string,
    
    trace =
      trace,
    
    peak_pos =
      peak_pos,
    
    channel_map =
      inferred$map
  )
  
  
  start_index <- find_start(
    
    metrics =
      metrics,
    
    window =
      settings$window,
    
    min_peak_ratio =
      settings$min_peak_ratio
  )
  
  
  collapse_index <- find_collapse(
    
    metrics =
      metrics,
    
    start_index =
      start_index,
    
    window =
      settings$window,
    
    min_len_before_collapse =
      settings$min_len_before_collapse,
    
    min_relative_signal =
      settings$min_relative_signal,
    
    min_peak_ratio =
      settings$min_peak_ratio,
    
    bad_run_windows =
      settings$bad_run_windows,
    
    absolute_max_base_index =
      settings$absolute_max_base_index
  )
  
  
  raw_length <-
    nchar(
      seq_string
    )
  
  
  if (
    is.na(start_index)
  ) {
    
    end_index <-
      NA_integer_
    
    trimmed_seq <-
      ""
    
    reason <-
      "NO_GOOD_START"
    
  } else {
    
    end_by_collapse <-
      if (
        !is.na(collapse_index)
      ) {
        
        collapse_index - 1
        
      } else {
        
        raw_length
      }
    
    
    end_by_cap <- min(
      raw_length,
      settings$absolute_max_base_index
    )
    
    
    end_index <- min(
      end_by_collapse,
      end_by_cap
    )
    
    
    reason <-
      if (
        !is.na(collapse_index) &&
        end_by_collapse <= end_by_cap
      ) {
        
        "AB1_SIGNAL_COLLAPSE"
        
      } else {
        
        "ABSOLUTE_CAP"
      }
    
    
    if (
      end_index >= start_index
    ) {
      
      trimmed_seq <- substr(
        seq_string,
        start_index,
        end_index
      )
      
    } else {
      
      trimmed_seq <- ""
    }
  }
  
  
  trimmed_length <-
    nchar(
      trimmed_seq
    )
  
  
  trimmed_metrics <-
    if (
      !is.na(start_index) &&
      !is.na(end_index) &&
      end_index >= start_index
    ) {
      
      metrics[
        start_index:end_index,
      ]
      
    } else {
      
      metrics[
        0,
      ]
    }
  
  
  status <- "OK"
  
  
  if (
    trimmed_length == 0
  ) {
    
    status <-
      "FAILED_TRIMMING"
  }
  
  
  if (
    trimmed_length > 0 &&
    trimmed_length <
    settings$min_usable_len
  ) {
    
    status <-
      "SHORT_AFTER_TRIMMING"
  }
  
  
  mean_signal <-
    if (
      trimmed_length > 0
    ) {
      
      value <- mean(
        trimmed_metrics$called_signal,
        na.rm = TRUE
      )
      
      if (
        is.nan(value)
      ) {
        NA_real_
      } else {
        round(
          value,
          2
        )
      }
      
    } else {
      
      NA_real_
    }
  
  
  median_ratio <-
    if (
      trimmed_length > 0
    ) {
      
      finite_values <-
        trimmed_metrics$peak_ratio[
          is.finite(
            trimmed_metrics$peak_ratio
          )
        ]
      
      
      if (
        length(finite_values) > 0
      ) {
        
        round(
          median(
            finite_values,
            na.rm = TRUE
          ),
          2
        )
        
      } else {
        
        NA_real_
      }
      
    } else {
      
      NA_real_
    }
  
  
  summary <- data.frame(
    
    sample_id =
      sample_id,
    
    target =
      settings$target,
    
    forward_primer =
      settings$forward_primer,
    
    reverse_primer =
      settings$reverse_primer,
    
    raw_length =
      raw_length,
    
    trimmed_length =
      trimmed_length,
    
    trim_start =
      start_index,
    
    trim_end =
      end_index,
    
    collapse_index =
      collapse_index,
    
    reason =
      reason,
    
    expected_amplicon_len =
      settings$expected_amplicon_len,
    
    absolute_max_base_index =
      settings$absolute_max_base_index,
    
    channel_A =
      inferred$map[["A"]],
    
    channel_C =
      inferred$map[["C"]],
    
    channel_G =
      inferred$map[["G"]],
    
    channel_T =
      inferred$map[["T"]],
    
    channel_map_pct_called_is_max =
      round(
        inferred$score$pct_called_is_max,
        4
      ),
    
    mean_called_signal_trimmed =
      mean_signal,
    
    median_peak_ratio_trimmed =
      median_ratio,
    
    status =
      status,
    
    error_message =
      "",
    
    stringsAsFactors =
      FALSE
  )
  
  
  list(
    
    sample_id =
      sample_id,
    
    raw_seq =
      seq_string,
    
    seq =
      trimmed_seq,
    
    metrics =
      metrics,
    
    summary =
      summary
  )
}


# ============================================================
# Failure summary
# ============================================================

make_failure_summary <- function(
    sample_id,
    settings,
    error_message
) {
  
  data.frame(
    
    sample_id =
      sample_id,
    
    target =
      settings$target,
    
    forward_primer =
      settings$forward_primer,
    
    reverse_primer =
      settings$reverse_primer,
    
    raw_length =
      NA_integer_,
    
    trimmed_length =
      0,
    
    trim_start =
      NA_integer_,
    
    trim_end =
      NA_integer_,
    
    collapse_index =
      NA_integer_,
    
    reason =
      "READ_OR_PROCESSING_ERROR",
    
    expected_amplicon_len =
      settings$expected_amplicon_len,
    
    absolute_max_base_index =
      settings$absolute_max_base_index,
    
    channel_A =
      NA_integer_,
    
    channel_C =
      NA_integer_,
    
    channel_G =
      NA_integer_,
    
    channel_T =
      NA_integer_,
    
    channel_map_pct_called_is_max =
      NA_real_,
    
    mean_called_signal_trimmed =
      NA_real_,
    
    median_peak_ratio_trimmed =
      NA_real_,
    
    status =
      "ERROR",
    
    error_message =
      error_message,
    
    stringsAsFactors =
      FALSE
  )
}


# ============================================================
# FASTA helpers
# ============================================================

clean_fasta_name <- function(x) {
  
  x <- trimws(x)
  
  x <- gsub(
    "\\s+",
    "_",
    x
  )
  
  x <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    x
  )
  
  x
}


make_fasta <- function(
    records,
    include_metadata = FALSE,
    summary_df = NULL
) {
  
  out <- character()
  
  
  if (
    length(records) == 0
  ) {
    
    return("")
  }
  
  
  for (
    original_name in names(records)
  ) {
    
    record <- records[[original_name]]
    
    
    if (
      is.null(record) ||
      nchar(record$seq) == 0
    ) {
      
      next
    }
    
    
    final_name <-
      clean_fasta_name(
        record$final_name
      )
    
    
    header <- paste0(
      ">",
      final_name
    )
    
    
    if (
      include_metadata &&
      !is.null(summary_df)
    ) {
      
      sm <- summary_df[
        summary_df$sample_id ==
          original_name,
        ,
        drop = FALSE
      ]
      
      
      if (
        nrow(sm) == 1
      ) {
        
        header <- paste0(
          
          header,
          
          " target=",
          gsub(
            "\\s+",
            "_",
            sm$target
          ),
          
          " length=",
          sm$trimmed_length,
          
          " trim=",
          sm$trim_start,
          "-",
          sm$trim_end,
          
          " reason=",
          sm$reason
        )
      }
    }
    
    
    out <- c(
      out,
      header,
      wrap_sequence(
        record$seq,
        80
      )
    )
  }
  
  
  paste(
    out,
    collapse = "\n"
  )
}


# ============================================================
# Plot helper
# ============================================================

draw_result_plot <- function(
    result,
    settings
) {
  
  metrics <-
    result$metrics
  
  
  sm <-
    result$summary
  
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  
  on.exit(
    par(old_par)
  )
  
  
  par(
    mfrow = c(2, 1),
    mar = c(4, 5, 3, 2)
  )
  
  
  plot(
    
    metrics$index,
    
    metrics$called_signal,
    
    type = "l",
    
    xlab =
      "Base-call index",
    
    ylab =
      "Called-base signal",
    
    main =
      paste0(
        result$sample_id,
        " - AB1 called-base signal"
      )
  )
  
  
  if (
    !is.na(sm$trim_start)
  ) {
    
    abline(
      v = sm$trim_start,
      lty = 2
    )
  }
  
  
  if (
    !is.na(sm$trim_end)
  ) {
    
    abline(
      v = sm$trim_end,
      lty = 2
    )
  }
  
  
  abline(
    
    v =
      settings$absolute_max_base_index,
    
    lty =
      3
  )
  
  
  if (
    !is.na(sm$collapse_index)
  ) {
    
    abline(
      v = sm$collapse_index,
      lty = 4
    )
  }
  
  
  finite_ratio <-
    metrics$peak_ratio[
      is.finite(
        metrics$peak_ratio
      )
    ]
  
  
  ymax <-
    if (
      length(finite_ratio) > 0
    ) {
      
      max_value <- max(
        finite_ratio,
        na.rm = TRUE
      )
      
      max(
        5,
        min(
          10,
          max_value
        )
      )
      
    } else {
      
      5
    }
  
  
  plot(
    
    metrics$index,
    
    metrics$peak_ratio,
    
    type = "l",
    
    xlab =
      "Base-call index",
    
    ylab =
      "Called peak / second peak",
    
    main =
      paste0(
        result$sample_id,
        " - AB1 peak ratio"
      ),
    
    ylim =
      c(
        0,
        ymax
      )
  )
  
  
  abline(
    
    h =
      settings$min_peak_ratio,
    
    lty =
      2
  )
  
  
  if (
    !is.na(sm$trim_start)
  ) {
    
    abline(
      v = sm$trim_start,
      lty = 2
    )
  }
  
  
  if (
    !is.na(sm$trim_end)
  ) {
    
    abline(
      v = sm$trim_end,
      lty = 2
    )
  }
  
  
  abline(
    
    v =
      settings$absolute_max_base_index,
    
    lty =
      3
  )
  
  
  if (
    !is.na(sm$collapse_index)
  ) {
    
    abline(
      v = sm$collapse_index,
      lty = 4
    )
  }
}


# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  
  tags$head(
    
    tags$style(
      HTML("

        body {
          background-color: #f5f7fa;
        }

        .pipeline-container {
          max-width: 1400px;
          margin: auto;
        }

        .app-header {
          background: white;
          padding: 22px 28px;
          margin-top: 20px;
          margin-bottom: 18px;
          border-radius: 12px;
          border: 1px solid #e5e7eb;
        }

        .app-header h2 {
          margin-top: 0;
          margin-bottom: 5px;
        }

        .app-header p {
          margin: 0;
          color: #6b7280;
        }

        .panel-box {
          background: white;
          padding: 22px;
          border-radius: 12px;
          border: 1px solid #e5e7eb;
          margin-bottom: 18px;
        }

        .section-title {
          font-size: 20px;
          font-weight: 600;
          margin-bottom: 16px;
        }

        .summary-card {
          background: #f8fafc;
          border: 1px solid #e5e7eb;
          border-radius: 10px;
          padding: 15px;
          margin-bottom: 10px;
        }

        .summary-number {
          font-size: 28px;
          font-weight: 700;
        }

        .summary-label {
          color: #6b7280;
        }

        .btn {
          border-radius: 7px;
        }

        .btn-primary {
          font-weight: 600;
        }

        .download-group .btn {
          margin-right: 8px;
          margin-bottom: 8px;
        }

        .settings-note {
          color: #6b7280;
          font-size: 13px;
        }

        .status-ok {
          color: #15803d;
          font-weight: 600;
        }

        .status-warning {
          color: #b45309;
          font-weight: 600;
        }

        .status-error {
          color: #b91c1c;
          font-weight: 600;
        }
        
        .table-legend {
  margin-top: 12px;
  padding: 10px 14px;
  background: #f8fafc;
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  font-size: 12px;
  color: #4b5563;
}

.table-legend ul {
  margin-top: 7px;
  margin-bottom: 0;
  padding-left: 20px;
}

.table-legend li {
  margin-bottom: 3px;
}

table.dataTable td input {
  color: #111827 !important;
  background-color: #ffffff !important;
}

table.dataTable td input:focus {
  color: #111827 !important;
  background-color: #ffffff !important;
}

      ")
    )
  ),
  
  
  div(
    
    class =
      "pipeline-container",
    
    
    div(
      
      class =
        "app-header",
      
      h2(
        "Sanger Sequence Pipeline"
      ),
      
      p(
        "AB1 trimming, quality control, renaming and BLAST FASTA export"
      )
    ),
    
    
    tabsetPanel(
      
      id =
        "pipeline_step",
      
      type =
        "tabs",
      
      
      # ======================================================
      # STEP 1 - UPLOAD
      # ======================================================
      
      tabPanel(
        
        title =
          "1 · Upload",
        
        value =
          "upload",
        
        
        div(
          
          class =
            "panel-box",
          
          div(
            class =
              "section-title",
            
            "Upload raw Sanger sequences"
          ),
          
          
          fileInput(
            
            inputId =
              "ab1_files",
            
            label =
              "Select AB1 files",
            
            multiple =
              TRUE,
            
            accept =
              c(
                ".ab1",
                ".AB1"
              ),
            
            buttonLabel =
              "Select files",
            
            placeholder =
              "No AB1 files selected"
          ),
          
          
          p(
            class =
              "settings-note",
            
            "You can select or drag multiple .ab1 chromatogram files at once."
          ),
          
          
          br(),
          
          
          DTOutput(
            "uploaded_files_table"
          ),
          
          
          br(),
          
          
          actionButton(
            
            inputId =
              "to_settings",
            
            label =
              "Continue to Settings →",
            
            class =
              "btn-primary"
          )
        )
      ),
      
      
      # ======================================================
      # STEP 2 - SETTINGS
      # ======================================================
      
      tabPanel(
        
        title =
          "2 · Settings",
        
        value =
          "settings",
        
        
        fluidRow(
          
          column(
            
            width =
              6,
            
            
            div(
              
              class =
                "panel-box",
              
              div(
                class =
                  "section-title",
                
                "Assay information"
              ),
              
              
              selectInput(
                
                inputId =
                  "target",
                
                label =
                  "Target / Gene",
                
                choices =
                  c(
                    "ITS",
                    "LSU",
                    "TEF1 / EF1-alpha",
                    "RPB2",
                    "Beta-tubulin",
                    "CYP51",
                    "SDHB",
                    "IGS",
                    "Other"
                  ),
                
                selected =
                  "ITS"
              ),
              
              
              textInput(
                
                inputId =
                  "forward_primer",
                
                label =
                  "Forward primer",
                
                value =
                  ""
              ),
              
              
              textInput(
                
                inputId =
                  "reverse_primer",
                
                label =
                  "Reverse primer",
                
                value =
                  ""
              ),
              
              
              numericInput(
                
                inputId =
                  "expected_amplicon_len",
                
                label =
                  "Expected amplicon length (bp)",
                
                value =
                  650,
                
                min =
                  1,
                
                step =
                  1
              ),
              
              
              numericInput(
                
                inputId =
                  "absolute_max_base_index",
                
                label =
                  "Maximum sequence position (bp)",
                
                value =
                  680,
                
                min =
                  50,
                
                step =
                  1
              )
            )
          ),
          
          
          column(
            
            width =
              6,
            
            
            div(
              
              class =
                "panel-box",
              
              div(
                class =
                  "section-title",
                
                "Advanced trimming settings"
              ),
              
              
              numericInput(
                
                inputId =
                  "window",
                
                label =
                  "Window size",
                
                value =
                  25,
                
                min =
                  5,
                
                step =
                  1
              ),
              
              
              numericInput(
                
                inputId =
                  "min_peak_ratio",
                
                label =
                  "Minimum peak ratio",
                
                value =
                  3,
                
                min =
                  0,
                
                step =
                  0.1
              ),
              
              
              numericInput(
                
                inputId =
                  "min_relative_signal",
                
                label =
                  "Minimum relative signal",
                
                value =
                  0.20,
                
                min =
                  0,
                
                max =
                  1,
                
                step =
                  0.01
              ),
              
              
              numericInput(
                
                inputId =
                  "min_len_before_collapse",
                
                label =
                  "Minimum length before collapse search",
                
                value =
                  350,
                
                min =
                  1,
                
                step =
                  1
              ),
              
              
              numericInput(
                
                inputId =
                  "bad_run_windows",
                
                label =
                  "Consecutive bad windows",
                
                value =
                  12,
                
                min =
                  1,
                
                step =
                  1
              ),
              
              
              numericInput(
                
                inputId =
                  "min_usable_len",
                
                label =
                  "Minimum usable trimmed length",
                
                value =
                  400,
                
                min =
                  1,
                
                step =
                  1
              )
            )
          )
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          actionButton(
            
            inputId =
              "back_upload",
            
            label =
              "← Back"
          ),
          
          actionButton(
            
            inputId =
              "run_trimming",
            
            label =
              "Run trimming",
            
            class =
              "btn-primary"
          )
        )
      ),
      
      
      # ======================================================
      # STEP 3 - TRIM + QC
      # ======================================================
      
      tabPanel(
        
        title =
          "3 · Trim & QC",
        
        value =
          "qc",
        
        
        uiOutput(
          "qc_summary_cards"
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          div(
            class =
              "section-title",
            
            "Trimming results"
          ),
          
          DTOutput(
            "summary_table"
          ),
          
          br(),
          
          div(
            class = "table-legend",
            
            tags$strong("Table legend"),
            
            tags$ul(
              tags$li(
                tags$strong("Sample"),
                " – sample / AB1 file name."
              ),
              tags$li(
                tags$strong("Target"),
                " – amplified target or gene selected for this run."
              ),
              tags$li(
                tags$strong("Raw length"),
                " – number of bases in the original AB1 base-called sequence before trimming."
              ),
              tags$li(
                tags$strong("Trimmed length"),
                " – number of bases retained after trimming."
              ),
              tags$li(
                tags$strong("Start"),
                " – first base retained in the trimmed sequence."
              ),
              tags$li(
                tags$strong("End"),
                " – last base retained in the trimmed sequence."
              ),
              tags$li(
                tags$strong("Collapse"),
                " – first position where sustained signal / quality collapse was detected."
              ),
              tags$li(
                tags$strong("Reason"),
                " – criterion that determined the sequence end, such as AB1_SIGNAL_COLLAPSE or ABSOLUTE_CAP."
              ),
              tags$li(
                tags$strong("Median peak ratio"),
                " – median ratio between the called-base peak and the second-highest peak in the retained region."
              ),
              tags$li(
                tags$strong("Status"),
                " – final trimming status: OK, SHORT_AFTER_TRIMMING, FAILED_TRIMMING or ERROR."
              )
            )
          )
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          div(
            class =
              "section-title",
            
            "Sequence quality plot"
          ),
          
          
          selectInput(
            
            inputId =
              "plot_sample",
            
            label =
              "Sample",
            
            choices =
              NULL
          ),
          
          
          plotOutput(
            
            outputId =
              "qc_plot",
            
            height =
              "650px"
          )
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          actionButton(
            
            inputId =
              "back_settings",
            
            label =
              "← Back to Settings"
          ),
          
          actionButton(
            
            inputId =
              "to_rename",
            
            label =
              "Continue to Rename →",
            
            class =
              "btn-primary"
          )
        )
      ),
      
      
      # ======================================================
      # STEP 4 - RENAME
      # ======================================================
      
      tabPanel(
        
        title =
          "4 · Rename",
        
        value =
          "rename",
        
        
        div(
          
          class =
            "panel-box",
          
          div(
            class =
              "section-title",
            
            "Optional sequence renaming"
          ),
          
          
          p(
            "Edit the New name column manually, or upload a rename key."
          ),
          
          fluidRow(
            
            column(
              width = 6,
              
              fileInput(
                inputId = "rename_key_file",
                label = "Upload rename key",
                multiple = FALSE,
                accept = c(".xlsx", ".csv"),
                buttonLabel = "Select file",
                placeholder = "No rename key selected"
              )
            ),
            
            column(
              width = 6,
              
              br(),
              
              actionButton(
                inputId = "apply_rename_key",
                label = "Apply rename key",
                class = "btn-primary"
              )
            )
          ),
          
          p(
            class = "settings-note",
            "Required columns: old_id and new_id. Example: 265DMAA000 → FB26."
          ),
          
          uiOutput(
            "rename_key_status"
          ),
          
          br(),
          
          DTOutput(
            "rename_table"
          ),
          
          
          br(),
          
          
          uiOutput(
            "rename_validation"
          )
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          actionButton(
            
            inputId =
              "back_qc",
            
            label =
              "← Back to QC"
          ),
          
          actionButton(
            
            inputId =
              "to_export",
            
            label =
              "Continue to Export →",
            
            class =
              "btn-primary"
          )
        )
      ),
      
      
      # ======================================================
      # STEP 5 - EXPORT
      # ======================================================
      
      tabPanel(
        
        title =
          "5 · Export",
        
        value =
          "export",
        
        
        div(
          
          class =
            "panel-box",
          
          div(
            class =
              "section-title",
            
            "Run summary"
          ),
          
          uiOutput(
            "export_summary"
          )
        ),
        
        
        div(
          
          class =
            "panel-box download-group",
          
          div(
            class =
              "section-title",
            
            "Download results"
          ),
          
          
          uiOutput(
            "download_buttons"
          )
        ),
        
        
        div(
          
          class =
            "panel-box",
          
          actionButton(
            
            inputId =
              "back_rename",
            
            label =
              "← Back to Rename"
          ),
          
          actionButton(
            
            inputId =
              "reset_pipeline",
            
            label =
              "Start new run"
          )
        )
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================

server <- function(
    input,
    output,
    session
) {
  
  
  # ==========================================================
  # Reactive storage
  # ==========================================================
  
  rv <- reactiveValues(
    
    results =
      list(),
    
    summary =
      NULL,
    
    rename =
      NULL,
    
    settings =
      NULL
  )
  
  
  # ==========================================================
  # Upload table
  # ==========================================================
  
  output$uploaded_files_table <- renderDT({
    
    req(
      input$ab1_files
    )
    
    
    df <- data.frame(
      
      File =
        input$ab1_files$name,
      
      Size_KB =
        round(
          input$ab1_files$size /
            1024,
          1
        ),
      
      stringsAsFactors =
        FALSE
    )
    
    
    datatable(
      
      df,
      
      rownames =
        FALSE,
      
      options =
        list(
          pageLength = 15,
          dom = "tip"
        )
    )
  })
  
  
  # ==========================================================
  # Navigation
  # ==========================================================
  
  observeEvent(
    input$to_settings,
    {
      
      if (
        is.null(input$ab1_files) ||
        nrow(input$ab1_files) == 0
      ) {
        
        showNotification(
          "Please upload at least one AB1 file.",
          type = "error"
        )
        
        return()
      }
      
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "settings"
      )
    }
  )
  
  
  observeEvent(
    input$back_upload,
    {
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "upload"
      )
    }
  )
  
  
  observeEvent(
    input$back_settings,
    {
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "settings"
      )
    }
  )
  
  
  observeEvent(
    input$to_rename,
    {
      
      req(
        rv$summary
      )
      
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "rename"
      )
    }
  )
  
  
  observeEvent(
    input$back_qc,
    {
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "qc"
      )
    }
  )
  
  
  observeEvent(
    input$back_rename,
    {
      
      updateTabsetPanel(
        session,
        "pipeline_step",
        selected = "rename"
      )
    }
  )
  
  
  # ==========================================================
  # Run trimming
  # ==========================================================
  
  observeEvent(
    input$run_trimming,
    {
      
      req(
        input$ab1_files
      )
      
      
      settings <- list(
        
        target =
          input$target,
        
        forward_primer =
          input$forward_primer,
        
        reverse_primer =
          input$reverse_primer,
        
        expected_amplicon_len =
          as.integer(
            input$expected_amplicon_len
          ),
        
        absolute_max_base_index =
          as.integer(
            input$absolute_max_base_index
          ),
        
        window =
          as.integer(
            input$window
          ),
        
        min_peak_ratio =
          as.numeric(
            input$min_peak_ratio
          ),
        
        min_relative_signal =
          as.numeric(
            input$min_relative_signal
          ),
        
        min_len_before_collapse =
          as.integer(
            input$min_len_before_collapse
          ),
        
        bad_run_windows =
          as.integer(
            input$bad_run_windows
          ),
        
        min_usable_len =
          as.integer(
            input$min_usable_len
          )
      )
      
      
      rv$settings <-
        settings
      
      
      all_results <- list()
      
      summaries <-
        list()
      
      
      files <-
        input$ab1_files
      
      
      withProgress(
        
        message =
          "Processing AB1 files",
        
        value =
          0,
        
        {
          
          
          for (
            i in seq_len(
              nrow(files)
            )
          ) {
            
            
            sample_id <- sub(
              
              "\\.ab1$",
              
              "",
              
              files$name[i],
              
              ignore.case =
                TRUE
            )
            
            
            incProgress(
              
              1 / nrow(files),
              
              detail =
                paste(
                  "Processing",
                  files$name[i]
                )
            )
            
            
            result <- tryCatch(
              
              {
                
                trim_one_ab1(
                  
                  file =
                    files$datapath[i],
                  
                  sample_id =
                    sample_id,
                  
                  settings =
                    settings
                )
              },
              
              error = function(e) {
                structure(
                  
                  list(
                    error =
                      conditionMessage(e)
                  ),
                  
                  class =
                    "ab1_error"
                )
              }
            )
            
            
            if (
              inherits(
                result,
                "ab1_error"
              )
            ) {
              
              summaries[[sample_id]] <- make_failure_summary(
                
                sample_id =
                  sample_id,
                
                settings =
                  settings,
                
                error_message =
                  result$error
              )
              
            } else {
              
              all_results[[sample_id]] <- result
              
              
              summaries[[sample_id]] <- result$summary
            }
          }
        }
      )
      
      
      rv$results <-
        all_results
      
      
      rv$summary <- do.call(
        rbind,
        summaries
      )
      
      
      rownames(
        rv$summary
      ) <- NULL
      
      
      rv$rename <- data.frame(
        
        Original_name =
          rv$summary$sample_id,
        
        New_name =
          rv$summary$sample_id,
        
        stringsAsFactors =
          FALSE
      )
      
      
      successful_samples <-
        names(
          rv$results
        )
      
      
      updateSelectInput(
        
        session,
        
        "plot_sample",
        
        choices =
          successful_samples,
        
        selected =
          if (
            length(successful_samples) > 0
          ) {
            successful_samples[1]
          } else {
            character(0)
          }
      )
      
      
      updateTabsetPanel(
        
        session,
        
        "pipeline_step",
        
        selected =
          "qc"
      )
      
      
      showNotification(
        
        paste(
          "Processing complete:",
          nrow(rv$summary),
          "samples"
        ),
        
        type =
          "message"
      )
    }
  )
  
  
  # ==========================================================
  # QC summary cards
  # ==========================================================
  
  output$qc_summary_cards <- renderUI({
    
    req(
      rv$summary
    )
    
    
    total <-
      nrow(
        rv$summary
      )
    
    
    ok <- sum(
      rv$summary$status == "OK",
      na.rm = TRUE
    )
    
    
    warning <- sum(
      rv$summary$status ==
        "SHORT_AFTER_TRIMMING",
      na.rm = TRUE
    )
    
    
    failed <- sum(
      rv$summary$status %in%
        c(
          "FAILED_TRIMMING",
          "ERROR"
        ),
      na.rm = TRUE
    )
    
    
    fluidRow(
      
      column(
        
        3,
        
        div(
          
          class =
            "summary-card",
          
          div(
            class =
              "summary-number",
            total
          ),
          
          div(
            class =
              "summary-label",
            "Total samples"
          )
        )
      ),
      
      
      column(
        
        3,
        
        div(
          
          class =
            "summary-card",
          
          div(
            class =
              "summary-number",
            ok
          ),
          
          div(
            class =
              "summary-label",
            "OK"
          )
        )
      ),
      
      
      column(
        
        3,
        
        div(
          
          class =
            "summary-card",
          
          div(
            class =
              "summary-number",
            warning
          ),
          
          div(
            class =
              "summary-label",
            "Warnings"
          )
        )
      ),
      
      
      column(
        
        3,
        
        div(
          
          class =
            "summary-card",
          
          div(
            class =
              "summary-number",
            failed
          ),
          
          div(
            class =
              "summary-label",
            "Failed"
          )
        )
      )
    )
  })
  
  
  # ==========================================================
  # QC summary table
  # ==========================================================
  
  output$summary_table <- renderDT({
    
    req(
      rv$summary
    )
    
    
    display_df <- rv$summary[
      ,
      c(
        "sample_id",
        "target",
        "raw_length",
        "trimmed_length",
        "trim_start",
        "trim_end",
        "collapse_index",
        "reason",
        "median_peak_ratio_trimmed",
        "status"
      )
    ]
    
    
    names(display_df) <- c(
      
      "Sample",
      
      "Target",
      
      "Raw length",
      
      "Trimmed length",
      
      "Start",
      
      "End",
      
      "Collapse",
      
      "Reason",
      
      "Median peak ratio",
      
      "Status"
    )
    
    
    datatable(
      
      display_df,
      
      rownames =
        FALSE,
      
      filter =
        "top",
      
      options =
        list(
          pageLength = 15,
          scrollX = TRUE
        )
    )
  })
  
  
  # ==========================================================
  # QC plot
  # ==========================================================
  
  output$qc_plot <- renderPlot({
    
    req(
      input$plot_sample
    )
    
    
    req(
      rv$results[[input$plot_sample]]
    )
    
    
    draw_result_plot(
      
      result =
        rv$results[[input$plot_sample]],
      
      settings =
        rv$settings
    )
  })
  
  # ==========================================================
  # Rename key import
  # ==========================================================
  
  rename_key_data <- reactive({
    
    req(input$rename_key_file)
    
    file <- input$rename_key_file
    ext <- tolower(tools::file_ext(file$name))
    
    if (ext == "xlsx") {
      
      key <- readxl::read_excel(file$datapath)
      
    } else if (ext == "csv") {
      
      key <- read.csv(
        file$datapath,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    } else {
      
      stop("Rename key must be an XLSX or CSV file.")
    }
    
    key <- as.data.frame(
      key,
      stringsAsFactors = FALSE
    )
    
    required_columns <- c(
      "old_id",
      "new_id"
    )
    
    if (!all(required_columns %in% names(key))) {
      
      stop(
        "Rename key must contain columns named old_id and new_id."
      )
    }
    
    key$old_id <- trimws(
      as.character(key$old_id)
    )
    
    key$new_id <- trimws(
      as.character(key$new_id)
    )
    
    key <- key[
      key$old_id != "" &
        key$new_id != "",
      ,
      drop = FALSE
    ]
    
    key
  })
  
  
  observeEvent(
    input$apply_rename_key,
    {
      
      req(
        rv$rename
      )
      
      key <- tryCatch(
        
        rename_key_data(),
        
        error = function(e) {
          
          showNotification(
            conditionMessage(e),
            type = "error"
          )
          
          return(NULL)
        }
      )
      
      if (is.null(key)) {
        return()
      }
      
      matched_count <- 0
      
      for (i in seq_len(nrow(rv$rename))) {
        
        original <- rv$rename$Original_name[i]
        
        # First try exact match
        idx <- which(
          key$old_id == original
        )
        
        # If exact match fails, also allow filenames containing
        # extra text after the sample ID, e.g. 265DMAA000_customer
        if (length(idx) == 0) {
          
          idx <- which(
            startsWith(
              original,
              key$old_id
            )
          )
        }
        
        if (length(idx) == 1) {
          
          rv$rename$New_name[i] <-
            key$new_id[idx]
          
          matched_count <-
            matched_count + 1
        }
      }
      
      showNotification(
        paste(
          matched_count,
          "of",
          nrow(rv$rename),
          "sample names matched."
        ),
        type = "message"
      )
    }
  )
  
  # ==========================================================
  # Rename table
  # ==========================================================
  
  output$rename_table <- renderDT({
    
    req(
      rv$rename
    )
    
    
    datatable(
      
      rv$rename,
      
      rownames =
        FALSE,
      
      editable =
        list(
          
          target =
            "cell",
          
          disable =
            list(
              columns =
                c(0)
            )
        ),
      
      options =
        list(
          pageLength = 20,
          dom = "tip"
        )
    )
  })
  
  
  observeEvent(
    input$rename_table_cell_edit,
    {
      
      info <-
        input$rename_table_cell_edit
      
      
      row <-
        info$row
      
      
      col <-
        info$col
      
      
      value <-
        as.character(
          info$value
        )
      
      
      # Only New_name should be editable
      if (
        col == 1
      ) {
        
        rv$rename[
          row,
          "New_name"
        ] <- value
      }
    }
  )
  
  
  # ==========================================================
  # Validate renamed sequence names
  # ==========================================================
  
  rename_error <- reactive({
    
    req(
      rv$rename
    )
    
    
    new_names <-
      trimws(
        rv$rename$New_name
      )
    
    
    if (
      any(
        new_names == ""
      )
    ) {
      
      return(
        "One or more sequence names are empty."
      )
    }
    
    
    cleaned_names <-
      clean_fasta_name(
        new_names
      )
    
    
    if (
      anyDuplicated(
        cleaned_names
      ) > 0
    ) {
      
      duplicated_values <- unique(
        
        cleaned_names[
          duplicated(
            cleaned_names
          ) |
            duplicated(
              cleaned_names,
              fromLast = TRUE
            )
        ]
      )
      
      
      return(
        
        paste0(
          "Duplicate sequence name(s): ",
          paste(
            duplicated_values,
            collapse = ", "
          )
        )
      )
    }
    
    
    NULL
  })
  
  
  output$rename_validation <- renderUI({
    
    req(
      rv$rename
    )
    
    
    error <-
      rename_error()
    
    
    if (
      is.null(error)
    ) {
      
      div(
        
        class =
          "status-ok",
        
        "✓ Sequence names are valid."
      )
      
    } else {
      
      div(
        
        class =
          "status-error",
        
        paste0(
          "⚠ ",
          error
        )
      )
    }
  })
  
  
  # ==========================================================
  # Continue to export
  # ==========================================================
  
  observeEvent(
    input$to_export,
    {
      
      error <-
        rename_error()
      
      
      if (
        !is.null(error)
      ) {
        
        showNotification(
          
          error,
          
          type =
            "error"
        )
        
        return()
      }
      
      
      updateTabsetPanel(
        
        session,
        
        "pipeline_step",
        
        selected =
          "export"
      )
    }
  )
  
  
  # ==========================================================
  # Build export records
  # ==========================================================
  
  export_records <- reactive({
    
    req(
      rv$results,
      rv$rename
    )
    
    
    records <-
      list()
    
    
    for (
      original_name in names(
        rv$results
      )
    ) {
      
      
      result <-
        rv$results[[original_name]]
      
      
      rename_row <- rv$rename[
        rv$rename$Original_name ==
          original_name,
        ,
        drop = FALSE
      ]
      
      
      if (
        nrow(rename_row) == 1
      ) {
        
        final_name <-
          rename_row$New_name
        
      } else {
        
        final_name <-
          original_name
      }
      
      
      result$final_name <-
        final_name
      
      
      records[[original_name]] <- result
    }
    
    
    records
  })
  
  
  # ==========================================================
  # Export summary
  # ==========================================================
  
  export_summary_df <- reactive({
    
    req(
      rv$summary,
      rv$rename
    )
    
    
    df <-
      rv$summary
    
    
    match_index <- match(
      
      df$sample_id,
      
      rv$rename$Original_name
    )
    
    
    df$final_name <-
      rv$rename$New_name[
        match_index
      ]
    
    
    df
  })
  
  
  output$export_summary <- renderUI({
    
    req(
      rv$summary,
      rv$settings
    )
    
    
    valid_sequences <- sum(
      rv$summary$trimmed_length > 0,
      na.rm = TRUE
    )
    
    
    tagList(
      
      p(
        strong(
          "Target: "
        ),
        rv$settings$target
      ),
      
      
      p(
        strong(
          "Forward primer: "
        ),
        ifelse(
          rv$settings$forward_primer == "",
          "Not specified",
          rv$settings$forward_primer
        )
      ),
      
      
      p(
        strong(
          "Reverse primer: "
        ),
        ifelse(
          rv$settings$reverse_primer == "",
          "Not specified",
          rv$settings$reverse_primer
        )
      ),
      
      
      p(
        strong(
          "Samples processed: "
        ),
        nrow(rv$summary)
      ),
      
      
      p(
        strong(
          "Sequences available for export: "
        ),
        valid_sequences
      ),
      
      
      p(
        strong(
          "Expected amplicon: "
        ),
        paste0(
          rv$settings$expected_amplicon_len,
          " bp"
        )
      )
    )
  })
  
  
  # ==========================================================
  # Download buttons
  # ==========================================================
  
  output$download_buttons <- renderUI({
    
    req(
      rv$summary
    )
    
    
    error <-
      rename_error()
    
    
    if (
      !is.null(error)
    ) {
      
      return(
        
        div(
          
          class =
            "status-error",
          
          paste(
            "Downloads disabled:",
            error
          )
        )
      )
    }
    
    
    tagList(
      
      downloadButton(
        "download_blast_fasta",
        "BLAST FASTA"
      ),
      
      downloadButton(
        "download_full_fasta",
        "FASTA + metadata"
      ),
      
      downloadButton(
        "download_summary_csv",
        "Summary CSV"
      ),
      
      downloadButton(
        "download_all_zip",
        "Download all as ZIP"
      )
    )
  })
  
  
  # ==========================================================
  # BLAST FASTA
  # ==========================================================
  
  output$download_blast_fasta <- downloadHandler(
    
    filename = function() {
      
      target_name <- clean_fasta_name(
        rv$settings$target
      )
      
      
      paste0(
        target_name,
        "_BLAST.fasta"
      )
    },
    
    
    content = function(file) {
      
      fasta_text <- make_fasta(
        
        records =
          export_records(),
        
        include_metadata =
          FALSE
      )
      
      
      writeLines(
        
        fasta_text,
        
        con =
          file,
        
        useBytes =
          TRUE
      )
    }
  )
  
  
  # ==========================================================
  # Full FASTA
  # ==========================================================
  
  output$download_full_fasta <- downloadHandler(
    
    filename = function() {
      
      target_name <- clean_fasta_name(
        rv$settings$target
      )
      
      
      paste0(
        target_name,
        "_trimmed_sequences.fasta"
      )
    },
    
    
    content = function(file) {
      
      fasta_text <- make_fasta(
        
        records =
          export_records(),
        
        include_metadata =
          TRUE,
        
        summary_df =
          export_summary_df()
      )
      
      
      writeLines(
        
        fasta_text,
        
        con =
          file,
        
        useBytes =
          TRUE
      )
    }
  )
  
  
  # ==========================================================
  # Summary CSV
  # ==========================================================
  
  output$download_summary_csv <- downloadHandler(
    
    filename = function() {
      
      target_name <- clean_fasta_name(
        rv$settings$target
      )
      
      
      paste0(
        target_name,
        "_trim_summary.csv"
      )
    },
    
    
    content = function(file) {
      
      write.csv(
        
        export_summary_df(),
        
        file,
        
        row.names =
          FALSE,
        
        fileEncoding =
          "UTF-8"
      )
    }
  )
  
  
  # ==========================================================
  # ZIP with everything
  # ==========================================================
  
  output$download_all_zip <- downloadHandler(
    
    filename = function() {
      
      target_name <- clean_fasta_name(
        rv$settings$target
      )
      
      
      paste0(
        target_name,
        "_Sanger_pipeline_results.zip"
      )
    },
    
    
    content = function(file) {
      
      temp_dir <- tempfile(
        "sanger_pipeline_"
      )
      
      
      dir.create(
        temp_dir,
        recursive = TRUE
      )
      
      
      on.exit(
        
        unlink(
          temp_dir,
          recursive = TRUE,
          force = TRUE
        ),
        
        add = TRUE
      )
      
      
      target_name <- clean_fasta_name(
        rv$settings$target
      )
      
      
      # ------------------------------------------
      # FASTA for BLAST
      # ------------------------------------------
      
      blast_path <- file.path(
        
        temp_dir,
        
        paste0(
          target_name,
          "_BLAST.fasta"
        )
      )
      
      
      writeLines(
        
        make_fasta(
          
          records =
            export_records(),
          
          include_metadata =
            FALSE
        ),
        
        blast_path
      )
      
      
      # ------------------------------------------
      # FASTA with metadata
      # ------------------------------------------
      
      full_fasta_path <- file.path(
        
        temp_dir,
        
        paste0(
          target_name,
          "_trimmed_sequences.fasta"
        )
      )
      
      
      writeLines(
        
        make_fasta(
          
          records =
            export_records(),
          
          include_metadata =
            TRUE,
          
          summary_df =
            export_summary_df()
        ),
        
        full_fasta_path
      )
      
      
      # ------------------------------------------
      # CSV
      # ------------------------------------------
      
      csv_path <- file.path(
        
        temp_dir,
        
        paste0(
          target_name,
          "_trim_summary.csv"
        )
      )
      
      
      write.csv(
        
        export_summary_df(),
        
        csv_path,
        
        row.names =
          FALSE,
        
        fileEncoding =
          "UTF-8"
      )
      
      
      # ------------------------------------------
      # Plots
      # ------------------------------------------
      
      plot_dir <- file.path(
        
        temp_dir,
        
        "QC_plots"
      )
      
      
      dir.create(
        plot_dir
      )
      
      
      for (
        sample_name in names(
          rv$results
        )
      ) {
        
        result <-
          rv$results[[sample_name]]
        
        
        png_path <- file.path(
          
          plot_dir,
          
          paste0(
            clean_fasta_name(
              sample_name
            ),
            "_QC.png"
          )
        )
        
        
        png(
          
          filename =
            png_path,
          
          width =
            1500,
          
          height =
            900,
          
          res =
            140
        )
        
        
        draw_result_plot(
          
          result =
            result,
          
          settings =
            rv$settings
        )
        
        
        dev.off()
      }
      
      
      # ------------------------------------------
      # Run settings
      # ------------------------------------------
      
      settings_path <- file.path(
        
        temp_dir,
        
        "run_settings.txt"
      )
      
      
      settings_lines <- c(
        
        paste0(
          "Target: ",
          rv$settings$target
        ),
        
        paste0(
          "Forward primer: ",
          rv$settings$forward_primer
        ),
        
        paste0(
          "Reverse primer: ",
          rv$settings$reverse_primer
        ),
        
        paste0(
          "Expected amplicon length: ",
          rv$settings$expected_amplicon_len
        ),
        
        paste0(
          "Absolute max base index: ",
          rv$settings$absolute_max_base_index
        ),
        
        paste0(
          "Window: ",
          rv$settings$window
        ),
        
        paste0(
          "Minimum peak ratio: ",
          rv$settings$min_peak_ratio
        ),
        
        paste0(
          "Minimum relative signal: ",
          rv$settings$min_relative_signal
        ),
        
        paste0(
          "Minimum length before collapse: ",
          rv$settings$min_len_before_collapse
        ),
        
        paste0(
          "Bad run windows: ",
          rv$settings$bad_run_windows
        ),
        
        paste0(
          "Minimum usable length: ",
          rv$settings$min_usable_len
        )
      )
      
      
      writeLines(
        
        settings_lines,
        
        settings_path
      )
      
      
      # ------------------------------------------
      # ZIP
      # ------------------------------------------
      
      files_to_zip <- list.files(
        
        temp_dir,
        
        recursive =
          TRUE,
        
        full.names =
          TRUE
      )
      
      
      zip::zipr(
        
        zipfile =
          file,
        
        files =
          files_to_zip,
        
        root =
          temp_dir
      )
    }
  )
  
  
  # ==========================================================
  # Reset
  # ==========================================================
  
  observeEvent(
    input$reset_pipeline,
    {
      
      rv$results <-
        list()
      
      rv$summary <-
        NULL
      
      rv$rename <-
        NULL
      
      rv$settings <-
        NULL
      
      
      updateTabsetPanel(
        
        session,
        
        "pipeline_step",
        
        selected =
          "upload"
      )
      
      
      showNotification(
        
        "Ready for a new run.",
        
        type =
          "message"
      )
    }
  )
}


# ============================================================
# Launch application
# ============================================================

shinyApp(
  ui = ui,
  server = server
)