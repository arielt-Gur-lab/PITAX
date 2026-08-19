# ============================================================
# PITAX 3.0 Stage 1 - AB1 evidence audit helpers
# v3.0.0-alpha.3
# ============================================================
#
# Stage 1 remains observational.  alpha.2 corrects an assumption made in
# alpha.1: for a freshly read ABIF object, sangerseqR's peakPosMatrix[,1]
# contains the primary ABIF base-call positions (PLOC.2 + 1).  The other
# columns are not an A/C/G/T peak-position matrix until makeBaseCalls() is
# explicitly run.  PITAX does not run makeBaseCalls() in its established v2
# processing path.
#
# Therefore this audit treats:
#   * ABIF PLOC.2 (+1, matching sangerseqR's 1-based conversion) as the raw
#     primary base-call position source;
#   * ABIF PCON.2/PCON.1 as optional basecaller quality evidence;
#   * sangerseqR traceMatrix as canonical A,C,G,T signal columns;
#   * the existing PITAX peakPosMatrix[,1] and inferred channel map as the
#     legacy path to be validated, not changed.
#
# alpha.3 adds a same-length, PCON-only window comparison.  It is diagnostic:
# it never changes trimming boundaries, curated bases, BLAST input or taxonomic
# interpretation.
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

pitax_abif_data <- function(abif_object) {
  dat <- pitax_safe_slot(abif_object, "data")
  if (!is.null(dat)) return(dat)
  nested_abif <- pitax_safe_slot(abif_object, "abifRawData")
  pitax_safe_slot(nested_abif, "data")
}

pitax_abif_quality <- function(abif_object, n_bases) {
  dat <- pitax_abif_data(abif_object)
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

pitax_abif_primary_positions <- function(abif_object, n_bases, fallback = NULL) {
  dat <- pitax_abif_data(abif_object)
  raw <- if (!is.null(dat)) tryCatch(dat[["PLOC.2"]], error = function(e) NULL) else NULL

  if (!is.null(raw) && length(raw)) {
    # sangerseqR converts ABIF PLOC.2 to 1-based R coordinates with +1.
    return(list(
      tag = "ABIF PLOC.2 + 1",
      values = pitax_normalize_numeric(raw, n_bases) + 1
    ))
  }

  if (!is.null(fallback) && length(fallback)) {
    return(list(
      tag = "sangerseqR peakPosMatrix[,1] fallback",
      values = pitax_normalize_numeric(fallback, n_bases)
    ))
  }

  list(tag = "Unavailable", values = rep(NA_real_, n_bases))
}

pitax_abif_dye_order <- function(abif_object) {
  dat <- pitax_abif_data(abif_object)
  if (is.null(dat)) return("Unavailable")
  x <- tryCatch(dat[["FWO.1"]], error = function(e) NULL)
  if (is.null(x) || !length(x)) return("Unavailable")
  paste(as.character(x), collapse = "")
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

  # sangerseqR constructs traceMatrix in A,C,G,T order for ABIF input.
  setNames(1:4, bases)
}

pitax_map_text <- function(channel_map) {
  bases <- c("A", "C", "G", "T")
  vals <- suppressWarnings(as.integer(channel_map[bases]))
  paste(paste0(bases, "=", ifelse(is.na(vals), "NA", vals)), collapse = ", ")
}

pitax_maps_identical <- function(a, b) {
  bases <- c("A", "C", "G", "T")
  aa <- suppressWarnings(as.integer(a[bases]))
  bb <- suppressWarnings(as.integer(b[bases]))
  length(aa) == 4L && length(bb) == 4L && all(!is.na(aa)) && all(!is.na(bb)) && identical(aa, bb)
}

pitax_signal_metrics <- function(seq_string, trace, peak_positions, channel_map = NULL) {
  bases <- strsplit(toupper(as.character(seq_string)), "", fixed = TRUE)[[1]]
  n <- min(length(bases), length(peak_positions))
  if (n < 1L) return(data.frame())
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
  valid_cols <- !is.na(canonical_cols) & canonical_cols >= 1L & canonical_cols <= ncol(trace)
  if (!all(valid_cols)) return(out)

  for (i in seq_len(n)) {
    b <- bases[i]
    p <- pos[i]
    if (!(b %in% canonical) || is.na(p) || p < 1L || p > nrow(trace)) next

    vals <- suppressWarnings(as.numeric(trace[p, canonical_cols, drop = TRUE]))
    names(vals) <- canonical
    if (!length(vals) || all(!is.finite(vals))) next

    called <- vals[[b]]
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

pitax_pct_true <- function(x) {
  usable <- !is.na(x)
  if (!any(usable)) return(NA_real_)
  round(100 * mean(x[usable] %in% TRUE), 2)
}

pitax_pct_threshold <- function(x, threshold) {
  x <- suppressWarnings(as.numeric(x))
  usable <- is.finite(x)
  if (!any(usable)) return(NA_real_)
  round(100 * mean(x[usable] >= threshold), 2)
}

pitax_median_finite <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  round(stats::median(x), digits)
}

pitax_mean_finite <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  round(mean(x), digits)
}

pitax_auto_trim_bounds <- function(result) {
  st <- NA_integer_
  en <- NA_integer_
  if (!is.null(result$curation) && is.list(result$curation)) {
    if (!is.null(result$curation$auto_trim_start)) st <- suppressWarnings(as.integer(result$curation$auto_trim_start[1]))
    if (!is.null(result$curation$auto_trim_end)) en <- suppressWarnings(as.integer(result$curation$auto_trim_end[1]))
  }
  if ((!is.finite(st) || !is.finite(en)) && is.data.frame(result$summary) && nrow(result$summary)) {
    st <- suppressWarnings(as.integer(result$summary$trim_start[1]))
    en <- suppressWarnings(as.integer(result$summary$trim_end[1]))
  }
  len <- if (is.finite(st) && is.finite(en) && en >= st) en - st + 1L else NA_integer_
  list(start = st, end = en, length = len)
}

pitax_window_quality_metrics <- function(quality, start, end) {
  quality <- suppressWarnings(as.numeric(quality))
  start <- suppressWarnings(as.integer(start[1]))
  end <- suppressWarnings(as.integer(end[1]))
  if (!length(quality) || !is.finite(start) || !is.finite(end) || end < start || start < 1L || end > length(quality)) {
    return(data.frame(
      Coverage_percent = NA_real_, Median_quality = NA_real_, Mean_quality = NA_real_,
      Q20_percent = NA_real_, Q30_percent = NA_real_, stringsAsFactors = FALSE
    ))
  }
  x <- quality[seq.int(start, end)]
  usable <- is.finite(x)
  data.frame(
    Coverage_percent = round(100 * mean(usable), 2),
    Median_quality = pitax_median_finite(x),
    Mean_quality = pitax_mean_finite(x),
    Q20_percent = pitax_pct_threshold(x, 20),
    Q30_percent = pitax_pct_threshold(x, 30),
    stringsAsFactors = FALSE
  )
}

# Compare every contiguous PCON window with the same length as the established
# automatic trim.  Candidates require at least 90% quality coverage. Ranking is
# deterministic and deliberately simple: Q20, Q30, median, mean, coverage, then
# the smallest shift from the established start. Missing scores count against a
# candidate during ranking. The result is audit evidence only.
pitax_quality_window_proposal <- function(quality, auto_start, auto_end, min_coverage = 0.90) {
  quality <- suppressWarnings(as.numeric(quality))
  auto_start <- suppressWarnings(as.integer(auto_start[1]))
  auto_end <- suppressWarnings(as.integer(auto_end[1]))
  n <- length(quality)
  win_len <- if (is.finite(auto_start) && is.finite(auto_end) && auto_end >= auto_start) auto_end - auto_start + 1L else NA_integer_
  unavailable <- list(
    available = FALSE, method = "PCON same-length window; observational only",
    start = NA_integer_, end = NA_integer_, length = win_len,
    coverage_percent = NA_real_, median_quality = NA_real_, mean_quality = NA_real_,
    q20_percent = NA_real_, q30_percent = NA_real_, start_shift = NA_integer_
  )
  if (!is.finite(win_len) || win_len < 1L || n < win_len || !any(is.finite(quality))) return(unavailable)

  starts <- seq_len(n - win_len + 1L)
  candidates <- lapply(starts, function(st) {
    en <- st + win_len - 1L
    x <- quality[seq.int(st, en)]
    usable <- is.finite(x)
    coverage <- mean(usable)
    finite_x <- x[usable]
    data.frame(
      start = st,
      end = en,
      coverage = coverage,
      q20_score = sum(finite_x >= 20) / win_len,
      q30_score = sum(finite_x >= 30) / win_len,
      median_quality = if (length(finite_x)) stats::median(finite_x) else -Inf,
      mean_quality = if (length(finite_x)) mean(finite_x) else -Inf,
      shift = abs(st - auto_start),
      stringsAsFactors = FALSE
    )
  })
  candidates <- do.call(rbind, candidates)
  candidates <- candidates[candidates$coverage >= min_coverage, , drop = FALSE]
  if (!nrow(candidates)) return(unavailable)

  ord <- order(
    -candidates$q20_score, -candidates$q30_score,
    -candidates$median_quality, -candidates$mean_quality,
    -candidates$coverage, candidates$shift, candidates$start
  )
  best <- candidates[ord[1], , drop = FALSE]
  metrics <- pitax_window_quality_metrics(quality, best$start, best$end)
  list(
    available = TRUE, method = "PCON same-length window; observational only",
    start = as.integer(best$start), end = as.integer(best$end), length = as.integer(win_len),
    coverage_percent = metrics$Coverage_percent[1],
    median_quality = metrics$Median_quality[1], mean_quality = metrics$Mean_quality[1],
    q20_percent = metrics$Q20_percent[1], q30_percent = metrics$Q30_percent[1],
    start_shift = as.integer(best$start - auto_start)
  )
}

pitax_result_quality_window <- function(result) {
  bounds <- pitax_auto_trim_bounds(result)
  detail <- if (!is.null(result$ab1_evidence)) result$ab1_evidence$detail else NULL
  quality <- if (is.data.frame(detail) && "Basecaller_quality" %in% names(detail)) detail$Basecaller_quality else numeric()
  pitax_quality_window_proposal(quality, bounds$start, bounds$end)
}

pitax_add_trim_membership <- function(detail, result) {
  if (!is.data.frame(detail) || !nrow(detail)) return(detail)
  bounds <- pitax_auto_trim_bounds(result)
  proposal <- pitax_result_quality_window(result)
  pos <- if ("Position" %in% names(detail)) suppressWarnings(as.integer(detail$Position)) else seq_len(nrow(detail))
  detail$In_auto_trim <- if (is.finite(bounds$start) && is.finite(bounds$end)) pos >= bounds$start & pos <= bounds$end else NA
  detail$In_quality_proposed_window <- if (isTRUE(proposal$available)) pos >= proposal$start & pos <= proposal$end else NA
  detail
}

pitax_trim_window_comparison <- function(result) {
  bounds <- pitax_auto_trim_bounds(result)
  detail <- if (!is.null(result$ab1_evidence)) result$ab1_evidence$detail else NULL
  quality <- if (is.data.frame(detail) && "Basecaller_quality" %in% names(detail)) detail$Basecaller_quality else numeric()
  legacy <- pitax_window_quality_metrics(quality, bounds$start, bounds$end)
  proposal <- pitax_quality_window_proposal(quality, bounds$start, bounds$end)
  data.frame(
    Window = c("Legacy v2 auto trim", "PCON-only proposal"),
    Status = c("Active output", if (isTRUE(proposal$available)) "Observational only" else "Unavailable"),
    Start = c(bounds$start, proposal$start),
    End = c(bounds$end, proposal$end),
    Length = c(bounds$length, proposal$length),
    Quality_coverage_percent = c(legacy$Coverage_percent[1], proposal$coverage_percent),
    Median_quality = c(legacy$Median_quality[1], proposal$median_quality),
    Q20_percent = c(legacy$Q20_percent[1], proposal$q20_percent),
    Q30_percent = c(legacy$Q30_percent[1], proposal$q30_percent),
    Start_shift = c(0L, proposal$start_shift),
    stringsAsFactors = FALSE
  )
}

pitax_result_sample_id <- function(result, fallback = "") {
  x <- if (!is.null(result$sample_id)) as.character(result$sample_id)[1] else ""
  if (is.na(x) || !nzchar(x)) x <- as.character(fallback)[1]
  if (is.na(x)) "" else x
}

pitax_assert_export_identity <- function(result, selected_key) {
  selected_key <- as.character(selected_key)[1]
  actual <- pitax_result_sample_id(result, selected_key)
  if (!nzchar(selected_key) || !identical(actual, selected_key)) {
    stop(
      "Selected sample/export mismatch. Download blocked: selected key '",
      selected_key, "', result sample '", actual, "'."
    )
  }
  actual
}

pitax_evidence_detail_export <- function(result, selected_key) {
  sid <- pitax_assert_export_identity(result, selected_key)
  ev <- result$ab1_evidence
  if (is.null(ev) || !is.data.frame(ev$detail) || !nrow(ev$detail)) {
    return(data.frame(Sample_ID = sid, Message = "No Stage 1 evidence available.", stringsAsFactors = FALSE))
  }
  detail <- pitax_add_trim_membership(ev$detail, result)
  preferred <- c("Position", "Base", "In_auto_trim", "In_quality_proposed_window")
  detail <- detail[, c(intersect(preferred, names(detail)), setdiff(names(detail), preferred)), drop = FALSE]
  data.frame(Sample_ID = rep(sid, nrow(detail)), detail, check.names = FALSE, stringsAsFactors = FALSE)
}

pitax_result_key_for_sample <- function(results, sample_id) {
  if (is.null(results) || !length(results)) return("")
  sample_id <- as.character(sample_id)[1]
  keys <- names(results)
  if (sample_id %in% keys) return(sample_id)
  ids <- vapply(seq_along(results), function(i) pitax_result_sample_id(results[[i]], keys[i]), character(1))
  idx <- which(ids == sample_id)
  if (!length(idx)) "" else keys[idx[1]]
}

build_ab1_evidence <- function(sanger_object, seq_string, trace, legacy_peak_pos, legacy_channel_map, abif_object = NULL) {
  n <- nchar(seq_string)
  quality_source <- if (!is.null(abif_object)) abif_object else sanger_object
  quality <- pitax_abif_quality(quality_source, n)
  primary_pos <- pitax_abif_primary_positions(abif_object, n, fallback = legacy_peak_pos)
  canonical_channel_map <- pitax_canonical_channel_map(trace)

  legacy_metrics <- pitax_signal_metrics(seq_string, trace, legacy_peak_pos, legacy_channel_map)
  canonical_metrics <- pitax_signal_metrics(seq_string, trace, primary_pos$values, canonical_channel_map)

  n_rows <- min(n, nrow(legacy_metrics), nrow(canonical_metrics))
  if (n_rows < 1L) {
    detail <- data.frame()
  } else {
    legacy_pos <- suppressWarnings(as.numeric(legacy_peak_pos[seq_len(n_rows)]))
    raw_pos <- suppressWarnings(as.numeric(primary_pos$values[seq_len(n_rows)]))
    detail <- data.frame(
      Position = seq_len(n_rows),
      Base = legacy_metrics$base[seq_len(n_rows)],
      Basecaller_quality = quality$values[seq_len(n_rows)],
      Legacy_primary_peak_pos = legacy_pos,
      Raw_ABIF_primary_peak_pos = raw_pos,
      Primary_peak_pos_delta = raw_pos - legacy_pos,
      Legacy_called_signal = legacy_metrics$called_signal[seq_len(n_rows)],
      Legacy_called_is_max = legacy_metrics$called_is_max[seq_len(n_rows)],
      Canonical_called_signal = canonical_metrics$called_signal[seq_len(n_rows)],
      Canonical_best_alt_signal = canonical_metrics$best_alt_signal[seq_len(n_rows)],
      Canonical_called_to_alt_ratio = canonical_metrics$called_to_alt_ratio[seq_len(n_rows)],
      Canonical_best_channel = canonical_metrics$best_channel[seq_len(n_rows)],
      Canonical_called_is_max = canonical_metrics$called_is_max[seq_len(n_rows)],
      stringsAsFactors = FALSE
    )
  }

  comparable <- if (nrow(detail)) is.finite(detail$Legacy_primary_peak_pos) & is.finite(detail$Raw_ABIF_primary_peak_pos) else logical()
  changed <- if (any(comparable)) detail$Legacy_primary_peak_pos[comparable] != detail$Raw_ABIF_primary_peak_pos[comparable] else logical()
  delta <- if (any(comparable)) abs(detail$Primary_peak_pos_delta[comparable]) else numeric()
  quality_available <- if (nrow(detail)) sum(is.finite(detail$Basecaller_quality)) else 0L
  primary_available <- if (nrow(detail)) sum(is.finite(detail$Raw_ABIF_primary_peak_pos)) else 0L
  maps_match <- pitax_maps_identical(legacy_channel_map, canonical_channel_map)

  trace_cols <- colnames(trace)
  if (is.null(trace_cols) || !length(trace_cols)) trace_cols <- paste0("column", seq_len(ncol(trace)))

  summary <- data.frame(
    quality_tag = quality$tag,
    quality_scores_available = quality_available,
    quality_coverage_percent = if (nrow(detail)) round(100 * quality_available / nrow(detail), 2) else NA_real_,
    median_quality_raw = if (nrow(detail)) pitax_median_finite(detail$Basecaller_quality) else NA_real_,
    primary_position_source = primary_pos$tag,
    primary_positions_available = primary_available,
    primary_position_coverage_percent = if (nrow(detail)) round(100 * primary_available / nrow(detail), 2) else NA_real_,
    positions_compared = sum(comparable),
    positions_with_different_primary_position = if (any(comparable)) sum(changed, na.rm = TRUE) else 0L,
    different_primary_position_percent = if (any(comparable)) round(100 * mean(changed, na.rm = TRUE), 2) else NA_real_,
    median_absolute_primary_position_delta = if (length(delta)) pitax_median_finite(delta) else NA_real_,
    legacy_inferred_channel_map = pitax_map_text(legacy_channel_map),
    canonical_channel_map = pitax_map_text(canonical_channel_map),
    channel_maps_match = maps_match,
    abif_dye_order = pitax_abif_dye_order(abif_object),
    trace_columns = paste(trace_cols, collapse = " | "),
    legacy_called_is_max_percent = if (nrow(detail)) pitax_pct_true(detail$Legacy_called_is_max) else NA_real_,
    canonical_called_is_max_percent = if (nrow(detail)) pitax_pct_true(detail$Canonical_called_is_max) else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    schema = "ab1-evidence-audit-v2",
    mode = "observational-only",
    quality_tag = quality$tag,
    quality = quality$values,
    primary_position_source = primary_pos$tag,
    raw_primary_positions = primary_pos$values,
    canonical_channel_map = canonical_channel_map,
    detail = detail,
    summary = summary
  )
}

ab1_evidence_result_summary <- function(result) {
  base <- data.frame(
    Sample = pitax_result_sample_id(result),
    Evidence = "Unavailable",
    Quality_tag = "Unavailable",
    Quality_coverage_percent = NA_real_,
    Auto_trim_start = NA_integer_,
    Auto_trim_end = NA_integer_,
    Auto_trim_length = NA_integer_,
    Median_quality_auto_trim = NA_real_,
    Q20_auto_trim_percent = NA_real_,
    Q30_auto_trim_percent = NA_real_,
    Quality_window_start = NA_integer_,
    Quality_window_end = NA_integer_,
    Quality_window_length = NA_integer_,
    Quality_window_start_shift = NA_integer_,
    Median_quality_quality_window = NA_real_,
    Q20_quality_window_percent = NA_real_,
    Q30_quality_window_percent = NA_real_,
    Primary_position_source = "",
    Primary_position_coverage_percent = NA_real_,
    Primary_position_difference_percent = NA_real_,
    Median_primary_position_delta = NA_real_,
    Legacy_map = "",
    Canonical_map = "",
    Channel_maps_match = NA,
    Legacy_called_is_max_percent = NA_real_,
    Canonical_called_is_max_percent = NA_real_,
    Canonical_called_is_max_auto_trim_percent = NA_real_,
    Median_called_to_alt_ratio_auto_trim = NA_real_,
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
  base$Primary_position_source <- as.character(sm$primary_position_source)
  base$Primary_position_coverage_percent <- suppressWarnings(as.numeric(sm$primary_position_coverage_percent))
  base$Primary_position_difference_percent <- suppressWarnings(as.numeric(sm$different_primary_position_percent))
  base$Median_primary_position_delta <- suppressWarnings(as.numeric(sm$median_absolute_primary_position_delta))
  base$Legacy_map <- as.character(sm$legacy_inferred_channel_map)
  base$Canonical_map <- as.character(sm$canonical_channel_map)
  base$Channel_maps_match <- as.logical(sm$channel_maps_match)
  base$Legacy_called_is_max_percent <- suppressWarnings(as.numeric(sm$legacy_called_is_max_percent))
  base$Canonical_called_is_max_percent <- suppressWarnings(as.numeric(sm$canonical_called_is_max_percent))

  bounds <- pitax_auto_trim_bounds(result)
  proposal <- pitax_result_quality_window(result)
  base$Auto_trim_start <- bounds$start
  base$Auto_trim_end <- bounds$end
  base$Auto_trim_length <- bounds$length
  base$Quality_window_start <- proposal$start
  base$Quality_window_end <- proposal$end
  base$Quality_window_length <- proposal$length
  base$Quality_window_start_shift <- proposal$start_shift
  base$Median_quality_quality_window <- proposal$median_quality
  base$Q20_quality_window_percent <- proposal$q20_percent
  base$Q30_quality_window_percent <- proposal$q30_percent

  detail <- ev$detail
  if (is.data.frame(detail) && nrow(detail)) {
    st <- bounds$start
    en <- bounds$end
    if (is.finite(st) && is.finite(en) && en >= st) {
      idx <- seq.int(max(1L, st), min(nrow(detail), en))
      if (length(idx)) {
        if ("Basecaller_quality" %in% names(detail)) {
          q <- detail$Basecaller_quality[idx]
          base$Median_quality_auto_trim <- pitax_median_finite(q)
          base$Q20_auto_trim_percent <- pitax_pct_threshold(q, 20)
          base$Q30_auto_trim_percent <- pitax_pct_threshold(q, 30)
        }
        if ("Canonical_called_is_max" %in% names(detail)) {
          base$Canonical_called_is_max_auto_trim_percent <- pitax_pct_true(detail$Canonical_called_is_max[idx])
        }
        if ("Canonical_called_to_alt_ratio" %in% names(detail)) {
          base$Median_called_to_alt_ratio_auto_trim <- pitax_median_finite(detail$Canonical_called_to_alt_ratio[idx], digits = 3L)
        }
      }
    }
  }

  base
}

ab1_evidence_run_summary <- function(results) {
  if (is.null(results) || !length(results)) return(data.frame())
  keys <- names(results)
  rows <- lapply(seq_along(results), function(i) {
    x <- ab1_evidence_result_summary(results[[i]])
    if (!nzchar(x$Sample[1])) x$Sample[1] <- keys[i]
    x
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
