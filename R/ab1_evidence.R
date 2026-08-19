# ============================================================
# PITAX 3.0 Stage 1 - AB1 evidence audit helpers
# ============================================================
#
# Stage 1 is deliberately observational. These helpers capture additional
# evidence that is already present in the AB1 / sangerseqR object so we can
# validate the assumptions used by the v2 trimming/QC engine on real lab
# chromatograms before changing any sequence decision rule.
#
# Nothing in this file changes trimming boundaries, curated bases, BLAST input
# or taxonomic interpretation in 3.0.0-alpha.1.1.
# ============================================================

pitax_safe_slot <- function(object, slot_name) {
  if (is.null(object)) return(NULL)
  slots <- tryCatch(methods::slotNames(object), error = function(e) character())
  if (!(slot_name %in% slots)) return(NULL)
  tryCatch(methods::slot(object, slot_name), error = function(e) NULL)
}

pitax_normalize_numeric <- function(x, n) {
  n <- max(0L, as.integer(n))
  if (n == 0L) return(numeric())
  if (is.null(x) || !length(x)) return(rep(NA_real_, n))
  out <- suppressWarnings(as.numeric(x))
  if (length(out) >= n) return(out[seq_len(n)])
  c(out, rep(NA_real_, n - length(out)))
}

pitax_abif_quality <- function(abif_object, n_bases) {
  # read.abif() returns an `abif` S4 object with a `data` list. Keep a
  # secondary compatibility path for wrappers that retain the ABIF object in
  # an `abifRawData` slot (for example some higher-level Sanger classes).
  dat <- pitax_safe_slot(abif_object, "data")
  if (is.null(dat)) {
    nested_abif <- pitax_safe_slot(abif_object, "abifRawData")
    dat <- pitax_safe_slot(nested_abif, "data")
  }
  if (is.null(dat)) {
    return(list(tag = "Unavailable", values = rep(NA_real_, n_bases)))
  }

  tags <- c("PCON.2", "PCON.1", "PCON")
  for (tag in tags) {
    value <- tryCatch(dat[[tag]], error = function(e) NULL)
    if (!is.null(value) && length(value)) {
      return(list(tag = tag, values = pitax_normalize_numeric(value, n_bases)))
    }
  }

  list(tag = "Unavailable", values = rep(NA_real_, n_bases))
}

pitax_canonical_channel_map <- function(x) {
  bases <- c("A", "C", "G", "T")
  if (is.null(x) || length(dim(x)) != 2L || ncol(x) < 4L) {
    return(setNames(rep(NA_integer_, 4L), bases))
  }

  cn <- colnames(x)
  if (!is.null(cn)) {
    cn_clean <- toupper(gsub("[^ACGT]", "", trimws(as.character(cn))))
    if (all(bases %in% cn_clean)) {
      return(setNames(match(bases, cn_clean), bases))
    }
  }

  # sangerseqR documents the four signal/peak columns as A, C, G, T.
  setNames(1:4, bases)
}

pitax_map_text <- function(channel_map) {
  bases <- c("A", "C", "G", "T")
  vals <- suppressWarnings(as.integer(channel_map[bases]))
  paste(paste0(bases, "=", ifelse(is.na(vals), "NA", vals)), collapse = ", ")
}

pitax_called_peak_positions <- function(seq_string, peak_pos_matrix) {
  bases <- strsplit(toupper(as.character(seq_string)), "", fixed = TRUE)[[1]]
  n <- length(bases)
  out <- rep(NA_real_, n)
  if (is.null(peak_pos_matrix) || length(dim(peak_pos_matrix)) != 2L || !nrow(peak_pos_matrix)) return(out)

  pm <- as.matrix(peak_pos_matrix)
  cmap <- pitax_canonical_channel_map(pm)
  rows <- seq_len(min(n, nrow(pm)))
  cols <- suppressWarnings(as.integer(cmap[bases[rows]]))
  valid <- !is.na(cols) & cols >= 1L & cols <= ncol(pm)
  if (any(valid)) {
    rr <- rows[valid]
    cc <- cols[valid]
    out[rr] <- suppressWarnings(as.numeric(pm[cbind(rr, cc)]))
  }
  out
}

pitax_signal_metrics <- function(seq_string, trace, peak_positions, channel_map = NULL) {
  bases <- strsplit(toupper(as.character(seq_string)), "", fixed = TRUE)[[1]]
  n <- min(length(bases), length(peak_positions))
  bases <- bases[seq_len(n)]
  pos <- suppressWarnings(as.integer(peak_positions[seq_len(n)]))

  out <- data.frame(
    index = seq_len(n),
    base = bases,
    peak_pos = pos,
    called_signal = rep(NA_real_, n),
    best_alt_signal = rep(NA_real_, n),
    called_to_alt_ratio = rep(NA_real_, n),
    best_channel = rep("", n),
    called_is_max = rep(NA, n),
    stringsAsFactors = FALSE
  )

  if (is.null(trace) || length(dim(trace)) != 2L || nrow(trace) < 1L || ncol(trace) < 4L) return(out)
  trace <- as.matrix(trace)
  if (is.null(channel_map)) channel_map <- pitax_canonical_channel_map(trace)
  canonical <- c("A", "C", "G", "T")
  canonical_cols <- suppressWarnings(as.integer(channel_map[canonical]))
  usable_channels <- which(!is.na(canonical_cols) & canonical_cols >= 1L & canonical_cols <= ncol(trace))
  if (!length(usable_channels)) return(out)

  for (i in seq_len(n)) {
    b <- bases[i]
    p <- pos[i]
    if (!(b %in% names(channel_map))) next
    called_col <- suppressWarnings(as.integer(channel_map[[b]]))
    if (is.na(called_col) || is.na(p) || p < 1L || p > nrow(trace) || called_col < 1L || called_col > ncol(trace)) next

    vals <- suppressWarnings(as.numeric(trace[p, canonical_cols[usable_channels], drop = TRUE]))
    names(vals) <- canonical[usable_channels]
    if (!length(vals) || all(!is.finite(vals))) next

    called <- suppressWarnings(as.numeric(trace[p, called_col]))
    alt <- vals[names(vals) != b]
    best_alt <- if (length(alt) && any(is.finite(alt))) max(alt, na.rm = TRUE) else NA_real_
    best_name <- if (any(is.finite(vals))) names(vals)[which.max(replace(vals, !is.finite(vals), -Inf))] else ""

    out$called_signal[i] <- called
    out$best_alt_signal[i] <- best_alt
    if (is.finite(called) && is.finite(best_alt)) {
      out$called_to_alt_ratio[i] <- if (best_alt > 0) called / best_alt else if (called > 0) Inf else NA_real_
      out$called_is_max[i] <- called >= max(vals, na.rm = TRUE)
    }
    out$best_channel[i] <- best_name
  }

  out
}

pitax_peak_amplitude_metrics <- function(seq_string, peak_amp_matrix) {
  bases <- strsplit(toupper(as.character(seq_string)), "", fixed = TRUE)[[1]]
  n <- length(bases)
  out <- data.frame(
    peakamp_called_signal = rep(NA_real_, n),
    peakamp_best_alt_signal = rep(NA_real_, n),
    peakamp_called_to_alt_ratio = rep(NA_real_, n),
    peakamp_best_channel = rep("", n),
    peakamp_called_is_max = rep(NA, n),
    stringsAsFactors = FALSE
  )
  if (is.null(peak_amp_matrix) || length(dim(peak_amp_matrix)) != 2L || !nrow(peak_amp_matrix) || ncol(peak_amp_matrix) < 4L) return(out)

  pam <- as.matrix(peak_amp_matrix)
  cmap <- pitax_canonical_channel_map(pam)
  canonical <- c("A", "C", "G", "T")
  canonical_cols <- suppressWarnings(as.integer(cmap[canonical]))
  n_use <- min(n, nrow(pam))

  for (i in seq_len(n_use)) {
    b <- bases[i]
    if (!(b %in% names(cmap))) next
    called_col <- suppressWarnings(as.integer(cmap[[b]]))
    if (is.na(called_col) || called_col < 1L || called_col > ncol(pam)) next
    vals <- suppressWarnings(as.numeric(pam[i, canonical_cols, drop = TRUE]))
    names(vals) <- canonical
    if (!length(vals) || all(!is.finite(vals))) next
    called <- suppressWarnings(as.numeric(pam[i, called_col]))
    alt <- vals[names(vals) != b]
    best_alt <- if (length(alt) && any(is.finite(alt))) max(alt, na.rm = TRUE) else NA_real_

    out$peakamp_called_signal[i] <- called
    out$peakamp_best_alt_signal[i] <- best_alt
    if (is.finite(called) && is.finite(best_alt)) {
      out$peakamp_called_to_alt_ratio[i] <- if (best_alt > 0) called / best_alt else if (called > 0) Inf else NA_real_
      out$peakamp_called_is_max[i] <- called >= max(vals, na.rm = TRUE)
    }
    out$peakamp_best_channel[i] <- if (any(is.finite(vals))) names(vals)[which.max(replace(vals, !is.finite(vals), -Inf))] else ""
  }

  out
}

pitax_pct_true <- function(x) {
  usable <- !is.na(x)
  if (!any(usable)) return(NA_real_)
  round(100 * mean(x[usable] %in% TRUE), 2)
}

pitax_median_finite <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  round(stats::median(x), digits)
}

build_ab1_evidence <- function(sanger_object, seq_string, trace, legacy_peak_pos, legacy_channel_map, abif_object = NULL) {
  n <- nchar(seq_string)
  peak_pos_matrix <- pitax_safe_slot(sanger_object, "peakPosMatrix")
  peak_amp_matrix <- pitax_safe_slot(sanger_object, "peakAmpMatrix")
  quality_source <- if (!is.null(abif_object)) abif_object else sanger_object
  quality <- pitax_abif_quality(quality_source, n)

  documented_peak_pos <- pitax_called_peak_positions(seq_string, peak_pos_matrix)
  documented_channel_map <- pitax_canonical_channel_map(trace)

  legacy_metrics <- pitax_signal_metrics(seq_string, trace, legacy_peak_pos, legacy_channel_map)
  documented_metrics <- pitax_signal_metrics(seq_string, trace, documented_peak_pos, documented_channel_map)
  amp_metrics <- pitax_peak_amplitude_metrics(seq_string, peak_amp_matrix)

  n_rows <- min(n, nrow(legacy_metrics), nrow(documented_metrics), nrow(amp_metrics))
  if (n_rows < 1L) {
    detail <- data.frame()
  } else {
    detail <- data.frame(
      Position = seq_len(n_rows),
      Base = legacy_metrics$base[seq_len(n_rows)],
      Basecaller_quality = quality$values[seq_len(n_rows)],
      Legacy_peak_pos = suppressWarnings(as.numeric(legacy_peak_pos[seq_len(n_rows)])),
      Called_base_peak_pos = documented_peak_pos[seq_len(n_rows)],
      Peak_pos_delta = documented_peak_pos[seq_len(n_rows)] - suppressWarnings(as.numeric(legacy_peak_pos[seq_len(n_rows)])),
      Legacy_called_signal = legacy_metrics$called_signal[seq_len(n_rows)],
      Legacy_called_is_max = legacy_metrics$called_is_max[seq_len(n_rows)],
      Called_base_signal = documented_metrics$called_signal[seq_len(n_rows)],
      Called_base_best_alt_signal = documented_metrics$best_alt_signal[seq_len(n_rows)],
      Called_base_to_alt_ratio = documented_metrics$called_to_alt_ratio[seq_len(n_rows)],
      Called_base_best_channel = documented_metrics$best_channel[seq_len(n_rows)],
      Called_base_is_max = documented_metrics$called_is_max[seq_len(n_rows)],
      PeakAmp_called_signal = amp_metrics$peakamp_called_signal[seq_len(n_rows)],
      PeakAmp_best_alt_signal = amp_metrics$peakamp_best_alt_signal[seq_len(n_rows)],
      PeakAmp_called_to_alt_ratio = amp_metrics$peakamp_called_to_alt_ratio[seq_len(n_rows)],
      PeakAmp_best_channel = amp_metrics$peakamp_best_channel[seq_len(n_rows)],
      PeakAmp_called_is_max = amp_metrics$peakamp_called_is_max[seq_len(n_rows)],
      stringsAsFactors = FALSE
    )
  }

  comparable <- if (nrow(detail)) is.finite(detail$Legacy_peak_pos) & is.finite(detail$Called_base_peak_pos) else logical()
  changed <- if (any(comparable)) detail$Legacy_peak_pos[comparable] != detail$Called_base_peak_pos[comparable] else logical()
  delta <- if (any(comparable)) abs(detail$Peak_pos_delta[comparable]) else numeric()
  quality_available <- if (nrow(detail)) sum(is.finite(detail$Basecaller_quality)) else 0L
  call_pos_available <- if (nrow(detail)) sum(is.finite(detail$Called_base_peak_pos)) else 0L

  trace_cols <- colnames(trace)
  if (is.null(trace_cols) || !length(trace_cols)) trace_cols <- paste0("column", seq_len(ncol(trace)))
  peak_pos_cols <- if (!is.null(peak_pos_matrix)) colnames(peak_pos_matrix) else NULL
  if (is.null(peak_pos_cols) && !is.null(peak_pos_matrix)) peak_pos_cols <- paste0("column", seq_len(ncol(as.matrix(peak_pos_matrix))))
  peak_amp_cols <- if (!is.null(peak_amp_matrix)) colnames(peak_amp_matrix) else NULL
  if (is.null(peak_amp_cols) && !is.null(peak_amp_matrix)) peak_amp_cols <- paste0("column", seq_len(ncol(as.matrix(peak_amp_matrix))))

  summary <- data.frame(
    quality_tag = quality$tag,
    quality_scores_available = quality_available,
    quality_coverage_percent = if (nrow(detail)) round(100 * quality_available / nrow(detail), 2) else NA_real_,
    median_quality_raw = if (nrow(detail)) pitax_median_finite(detail$Basecaller_quality) else NA_real_,
    trace_columns = paste(trace_cols, collapse = " | "),
    peak_pos_columns = if (length(peak_pos_cols)) paste(peak_pos_cols, collapse = " | ") else "Unavailable",
    peak_amp_columns = if (length(peak_amp_cols)) paste(peak_amp_cols, collapse = " | ") else "Unavailable",
    legacy_inferred_channel_map = pitax_map_text(legacy_channel_map),
    documented_channel_map = pitax_map_text(documented_channel_map),
    called_base_peak_positions_available = call_pos_available,
    called_base_peak_position_coverage_percent = if (nrow(detail)) round(100 * call_pos_available / nrow(detail), 2) else NA_real_,
    positions_compared = sum(comparable),
    positions_with_different_peak_position = if (any(comparable)) sum(changed, na.rm = TRUE) else 0L,
    different_peak_position_percent = if (any(comparable)) round(100 * mean(changed, na.rm = TRUE), 2) else NA_real_,
    median_absolute_peak_position_delta = if (length(delta)) pitax_median_finite(delta) else NA_real_,
    legacy_called_is_max_percent = if (nrow(detail)) pitax_pct_true(detail$Legacy_called_is_max) else NA_real_,
    called_base_position_called_is_max_percent = if (nrow(detail)) pitax_pct_true(detail$Called_base_is_max) else NA_real_,
    peak_amp_called_is_max_percent = if (nrow(detail)) pitax_pct_true(detail$PeakAmp_called_is_max) else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    schema = "ab1-evidence-audit-v1",
    mode = "observational-only",
    quality_tag = quality$tag,
    quality = quality$values,
    peak_pos_matrix = peak_pos_matrix,
    peak_amp_matrix = peak_amp_matrix,
    documented_peak_pos = documented_peak_pos,
    documented_channel_map = documented_channel_map,
    detail = detail,
    summary = summary
  )
}

ab1_evidence_result_summary <- function(result) {
  base <- data.frame(
    Sample = if (!is.null(result$sample_id)) as.character(result$sample_id) else "",
    Evidence = "Unavailable",
    Quality_tag = "Unavailable",
    Quality_coverage_percent = NA_real_,
    Median_quality_auto_trim = NA_real_,
    Legacy_map = "",
    Documented_map = "",
    Peak_position_difference_percent = NA_real_,
    Median_peak_position_delta = NA_real_,
    Legacy_called_is_max_percent = NA_real_,
    Called_base_called_is_max_percent = NA_real_,
    PeakAmp_called_is_max_percent = NA_real_,
    stringsAsFactors = FALSE
  )

  ev <- result$ab1_evidence
  if (is.null(ev) || !is.list(ev)) return(base)
  if (!is.null(ev$error) && nzchar(as.character(ev$error)[1])) {
    base$Evidence <- paste0("Audit error: ", as.character(ev$error)[1])
    return(base)
  }
  if (!is.data.frame(ev$summary) || !nrow(ev$summary)) return(base)
  sm <- ev$summary[1, , drop = FALSE]
  base$Evidence <- "Captured"
  base$Quality_tag <- as.character(sm$quality_tag)
  base$Quality_coverage_percent <- suppressWarnings(as.numeric(sm$quality_coverage_percent))
  base$Legacy_map <- as.character(sm$legacy_inferred_channel_map)
  base$Documented_map <- as.character(sm$documented_channel_map)
  base$Peak_position_difference_percent <- suppressWarnings(as.numeric(sm$different_peak_position_percent))
  base$Median_peak_position_delta <- suppressWarnings(as.numeric(sm$median_absolute_peak_position_delta))
  base$Legacy_called_is_max_percent <- suppressWarnings(as.numeric(sm$legacy_called_is_max_percent))
  base$Called_base_called_is_max_percent <- suppressWarnings(as.numeric(sm$called_base_position_called_is_max_percent))
  base$PeakAmp_called_is_max_percent <- suppressWarnings(as.numeric(sm$peak_amp_called_is_max_percent))

  detail <- ev$detail
  if (is.data.frame(detail) && nrow(detail) && "Basecaller_quality" %in% names(detail)) {
    st <- NA_integer_
    en <- NA_integer_
    if (!is.null(result$curation) && is.list(result$curation)) {
      if (!is.null(result$curation$auto_trim_start)) st <- suppressWarnings(as.integer(result$curation$auto_trim_start[1]))
      if (!is.null(result$curation$auto_trim_end)) en <- suppressWarnings(as.integer(result$curation$auto_trim_end[1]))
    }
    if (!is.finite(st) || !is.finite(en)) {
      st <- suppressWarnings(as.integer(result$summary$trim_start[1]))
      en <- suppressWarnings(as.integer(result$summary$trim_end[1]))
    }
    if (is.finite(st) && is.finite(en) && en >= st) {
      idx <- seq.int(max(1L, st), min(nrow(detail), en))
      if (length(idx)) base$Median_quality_auto_trim <- pitax_median_finite(detail$Basecaller_quality[idx])
    }
  }

  base
}

ab1_evidence_run_summary <- function(results) {
  if (is.null(results) || !length(results)) return(data.frame())
  rows <- lapply(results, ab1_evidence_result_summary)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
