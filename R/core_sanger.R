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
      all(this_bad %in% TRUE)
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

    trace =
      trace,

    peak_pos =
      peak_pos,

    channel_map =
      inferred$map,
    
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
