# ============================================================
# Sequence / primer / chromatogram tools
# ============================================================

sanitize_dna <- function(x) {
  x <- toupper(gsub("\\s+", "", ifelse(is.null(x), "", x)))
  gsub("[^ACGTN]", "", x)
}

reverse_complement <- function(x) {
  x <- sanitize_dna(x)
  if (!nzchar(x)) return("")
  chars <- strsplit(x, "")[[1]]
  comp <- chartr("ACGTN", "TGCAN", chars)
  paste(rev(comp), collapse = "")
}

sequence_gc <- function(x) {
  x <- sanitize_dna(x)
  if (!nzchar(x)) return(NA_real_)
  b <- strsplit(x, "")[[1]]
  valid <- b %in% c("A", "C", "G", "T")
  if (!any(valid)) return(NA_real_)
  round(100 * mean(b[valid] %in% c("G", "C")), 2)
}

sequence_n_count <- function(x) {
  x <- sanitize_dna(x)
  if (!nzchar(x)) return(0L)
  sum(strsplit(x, "")[[1]] == "N")
}

empty_primer_match <- function(primer = "") {
  data.frame(
    orientation = NA_character_, start = NA_integer_, end = NA_integer_,
    identity = NA_real_, matches = NA_integer_, primer_length = nchar(primer),
    query = primer, target = NA_character_, search_from = NA_integer_,
    search_to = NA_integer_, stringsAsFactors = FALSE
  )
}

# Search a defined region for one defined sequence orientation.
# This avoids the old behaviour of choosing a high-scoring site anywhere in the read.
regional_primer_match <- function(sequence, query, search_from = 1L, search_to = NULL) {
  sequence <- sanitize_dna(sequence)
  query <- sanitize_dna(query)
  n <- nchar(sequence)
  L <- nchar(query)

  if (!nzchar(query) || !n || L > n) return(empty_primer_match(query))

  if (is.null(search_to) || is.na(search_to)) search_to <- n
  search_from <- max(1L, as.integer(search_from))
  search_to <- min(n, as.integer(search_to))
  last_start <- min(search_to - L + 1L, n - L + 1L)
  if (last_start < search_from) return(empty_primer_match(query))

  qchars <- strsplit(query, "")[[1]]
  best <- NULL
  for (i in seq.int(search_from, last_start)) {
    target <- substr(sequence, i, i + L - 1L)
    tchars <- strsplit(target, "")[[1]]
    comparable <- qchars != "N" & tchars != "N"
    denom <- sum(comparable)
    matches <- if (denom > 0) sum(qchars[comparable] == tchars[comparable]) else 0L
    identity <- if (denom > 0) 100 * matches / denom else 0
    row <- data.frame(
      orientation = "expected", start = i, end = i + L - 1L,
      identity = identity, matches = matches, primer_length = L,
      query = query, target = target, search_from = search_from,
      search_to = search_to, stringsAsFactors = FALSE
    )
    if (is.null(best) || row$identity > best$identity) best <- row
  }
  best$identity <- round(best$identity, 1)
  best
}

# Biological primer model for a single Sanger read.
# Forward read:
#   - Forward sequencing primer is the start anchor. It may be absent from the called bases
#     because extension begins immediately after the primer.
#   - The opposite Reverse primer site, if reached, is expected as reverse-complement(R).
# Reverse read:
#   - Reverse sequencing primer is the start anchor.
#   - The opposite Forward primer site, if reached, is expected as reverse-complement(F).
# Unknown:
#   - both biological models are evaluated and the model with stronger positional evidence wins.
primer_model_matches <- function(result, settings) {
  if (!isTRUE(settings$enable_primer_mapping)) return(data.frame())
  seq <- sanitize_dna(result$raw_seq)
  n <- nchar(seq)
  f <- sanitize_dna(settings$forward_primer_seq)
  r <- sanitize_dna(settings$reverse_primer_seq)
  direction <- ifelse(is.null(settings$sequencing_primer), "Forward", settings$sequencing_primer)
  expected_len <- suppressWarnings(as.integer(settings$expected_amplicon_len))
  if (is.na(expected_len) || expected_len < 1) expected_len <- n

  start_window <- min(n, max(80L, min(150L, round(expected_len * 0.20))))
  opposite_half_window <- max(100L, min(220L, round(expected_len * 0.25)))

  build_model <- function(dir) {
    rows <- list()
    if (dir == "Forward") {
      if (nzchar(f)) {
        m <- regional_primer_match(seq, f, 1L, start_window)
        rows[["Forward"]] <- cbind(
          data.frame(Primer="Forward", Role="Sequencing primer / read-start anchor",
                     Expected_sequence="as entered", stringsAsFactors=FALSE), m
        )
      }
      if (nzchar(r)) {
        q <- reverse_complement(r)
        center <- expected_len - nchar(r) + 1L
        lo <- max(start_window + 1L, center - opposite_half_window)
        hi <- min(n, center + opposite_half_window)
        # If the expected amplicon end is beyond the read, search the accessible downstream region.
        if (lo > hi) {
          lo <- max(start_window + 1L, floor(n * 0.35))
          hi <- n
        }
        m <- regional_primer_match(seq, q, lo, hi)
        rows[["Reverse"]] <- cbind(
          data.frame(Primer="Reverse", Role="Opposite primer site",
                     Expected_sequence="reverse-complement", stringsAsFactors=FALSE), m
        )
      }
    } else if (dir == "Reverse") {
      if (nzchar(r)) {
        m <- regional_primer_match(seq, r, 1L, start_window)
        rows[["Reverse"]] <- cbind(
          data.frame(Primer="Reverse", Role="Sequencing primer / read-start anchor",
                     Expected_sequence="as entered", stringsAsFactors=FALSE), m
        )
      }
      if (nzchar(f)) {
        q <- reverse_complement(f)
        center <- expected_len - nchar(f) + 1L
        lo <- max(start_window + 1L, center - opposite_half_window)
        hi <- min(n, center + opposite_half_window)
        if (lo > hi) {
          lo <- max(start_window + 1L, floor(n * 0.35))
          hi <- n
        }
        m <- regional_primer_match(seq, q, lo, hi)
        rows[["Forward"]] <- cbind(
          data.frame(Primer="Forward", Role="Opposite primer site",
                     Expected_sequence="reverse-complement", stringsAsFactors=FALSE), m
        )
      }
    }
    if (!length(rows)) return(data.frame())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }

  if (direction %in% c("Forward", "Reverse")) {
    out <- build_model(direction)
    attr(out, "model_direction") <- direction
    return(out)
  }

  fm <- build_model("Forward")
  rm <- build_model("Reverse")
  score_model <- function(x) {
    if (!nrow(x)) return(-Inf)
    ids <- x$identity[is.finite(x$identity)]
    if (!length(ids)) return(-Inf)
    # Prioritise the start-anchor evidence, then opposite-primer evidence.
    roles <- x$Role
    start_id <- ids[1]
    start_rows <- which(grepl("read-start", roles))
    if (length(start_rows) && is.finite(x$identity[start_rows[1]])) start_id <- x$identity[start_rows[1]]
    opp_rows <- which(roles == "Opposite primer site")
    opp_id <- if (length(opp_rows) && is.finite(x$identity[opp_rows[1]])) x$identity[opp_rows[1]] else 0
    start_id + 0.7 * opp_id
  }
  if (score_model(fm) >= score_model(rm)) {
    attr(fm, "model_direction") <- "Forward (inferred)"
    fm
  } else {
    attr(rm, "model_direction") <- "Reverse (inferred)"
    rm
  }
}

primer_match_table <- function(result, settings) {
  pm <- primer_model_matches(result, settings)
  if (!nrow(pm)) {
    return(data.frame(
      Primer=character(), Name=character(), Role=character(), Expected=character(),
      Start=integer(), End=integer(), Identity_percent=numeric(), Support=character(),
      In_trimmed=character(), Search_region=character(), stringsAsFactors=FALSE
    ))
  }

  name_for <- function(type) {
    val <- if (type == "Forward") settings$forward_primer else settings$reverse_primer
    if (is.null(val) || !nzchar(val)) type else val
  }
  sm <- result$summary
  out <- data.frame(
    Primer = pm$Primer,
    Name = vapply(pm$Primer, name_for, character(1)),
    Role = pm$Role,
    Expected = pm$Expected_sequence,
    Start = pm$start,
    End = pm$end,
    Identity_percent = pm$identity,
    Support = ifelse(is.na(pm$identity), "Not detected",
                     ifelse(pm$identity >= 95, "Strong",
                            ifelse(pm$identity >= 85, "Moderate", "Weak"))),
    In_trimmed = ifelse(is.na(pm$start), "No",
                        ifelse(!is.na(sm$trim_start) & !is.na(sm$trim_end) &
                                 pm$start >= sm$trim_start & pm$end <= sm$trim_end, "Yes", "No")),
    Search_region = ifelse(is.na(pm$search_from), NA_character_, paste0(pm$search_from, "-", pm$search_to)),
    stringsAsFactors = FALSE
  )
  attr(out, "model_direction") <- attr(pm, "model_direction")
  out
}

primer_alignment_text <- function(result, settings) {
  pm_raw <- primer_model_matches(result, settings)
  if (!nrow(pm_raw)) return("No primer sequences supplied.")
  model <- attr(pm_raw, "model_direction")
  pieces <- c(paste0("Sequencing model: ", model), "")
  for (i in seq_len(nrow(pm_raw))) {
    m <- pm_raw[i, , drop=FALSE]
    if (is.na(m$start)) {
      pieces <- c(pieces, paste0(m$Primer, ": no site evaluated/detected."), "")
      next
    }
    qchars <- strsplit(m$query, "")[[1]]
    tchars <- strsplit(m$target, "")[[1]]
    marks <- ifelse(qchars == tchars & qchars != "N" & tchars != "N", "|", " ")
    pieces <- c(
      pieces,
      paste0(m$Primer, " primer — ", m$Role),
      paste0("Search region: bases ", m$search_from, "-", m$search_to,
             " | best site: ", m$start, "-", m$end,
             " | ", m$identity, "% identity | expected ", m$Expected_sequence),
      paste0("Primer  ", m$query),
      paste0("        ", paste(marks, collapse="")),
      paste0("Read    ", m$target),
      ""
    )
  }
  paste(pieces, collapse="\n")
}

# Convert trace-sample coordinates to a continuous base-position coordinate.
# Peak positions map exactly to integer base positions; points between peaks are interpolated.
trace_to_base_x <- function(peak_pos, trace_idx) {
  good <- which(is.finite(peak_pos) & peak_pos >= 1)
  if (length(good) < 2) return(rep(NA_real_, length(trace_idx)))
  approx(x=peak_pos[good], y=good, xout=trace_idx, rule=2, ties="ordered")$y
}

# Main chromatogram: X axis is BASE POSITION, not raw trace sample index.
draw_chromatogram <- function(result, settings, start_base=1L, visible_bases=50L) {
  trace <- result$trace
  peak_pos <- result$peak_pos
  cmap <- result$channel_map
  n <- length(peak_pos)
  visible_bases <- max(10L, min(as.integer(visible_bases), n))
  from_base <- max(1L, min(as.integer(start_base), max(1L, n-visible_bases+1L)))
  to_base <- min(n, from_base + visible_bases - 1L)

  x_from <- max(1L, peak_pos[from_base])
  x_to <- min(nrow(trace), peak_pos[to_base])
  idx <- seq.int(x_from, x_to)
  x_base <- trace_to_base_x(peak_pos, idx)
  vals <- trace[idx, , drop=FALSE]
  ymax <- max(vals, na.rm=TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  plot(x_base, trace[idx, cmap[["A"]]], type="l", col="green3", lwd=1,
       xlab="Base position", ylab="Signal", xlim=c(from_base, to_base), ylim=c(0, ymax*1.20),
       xaxt="n", main=paste0(result$sample_id, " — bases ", from_base, "-", to_base))
  lines(x_base, trace[idx, cmap[["C"]]], col="blue3", lwd=1)
  lines(x_base, trace[idx, cmap[["G"]]], col="black", lwd=1)
  lines(x_base, trace[idx, cmap[["T"]]], col="red3", lwd=1)

  # X-axis resolution adapts to zoom.
  axis_step <- if (visible_bases <= 30) 1L else if (visible_bases <= 60) 5L else if (visible_bases <= 120) 10L else 20L
  at <- seq(ceiling(from_base/axis_step)*axis_step, to_base, by=axis_step)
  if (!length(at)) at <- c(from_base, to_base)
  axis(1, at=at, labels=at)
  legend("topright", legend=c("A","C","G","T"), col=c("green3","blue3","black","red3"), lty=1, horiz=TRUE, bty="n")

  # Nucleotide labels: every base when zoomed in, progressively sparser when zoomed out.
  bases <- strsplit(result$raw_seq, "")[[1]]
  label_step <- if (visible_bases <= 50) 1L else if (visible_bases <= 100) 2L else if (visible_bases <= 160) 5L else 10L
  label_idx <- seq(from_base, to_base, by=label_step)
  text(label_idx, rep(ymax*1.045, length(label_idx)), labels=bases[label_idx], cex=if (visible_bases<=50) 0.78 else 0.62)

  sm <- result$summary
  if (!is.na(sm$trim_start) && sm$trim_start >= from_base && sm$trim_start <= to_base) abline(v=sm$trim_start, lty=2, lwd=2)
  if (!is.na(sm$trim_end) && sm$trim_end >= from_base && sm$trim_end <= to_base) abline(v=sm$trim_end, lty=2, lwd=2)

  if (isTRUE(settings$enable_primer_mapping)) {
    pm <- primer_match_table(result, settings)
    if (nrow(pm)) {
      for (i in seq_len(nrow(pm))) {
        s <- pm$Start[i]; e <- pm$End[i]
        if (!is.na(s) && e >= from_base && s <= to_base) {
          ss <- max(s, from_base); ee <- min(e, to_base)
          col <- if (pm$Primer[i] == "Forward") "darkorange" else "purple"
          rect(ss, ymax*1.09, ee, ymax*1.17, border=col, lwd=2)
          text(mean(c(ss,ee)), ymax*1.13, labels=paste0(pm$Primer[i], " primer"), cex=0.68, col=col)
        }
      }
    }
  }
}

# Overview connects read length, trimming and biologically constrained primer matches.
draw_amplicon_overview <- function(result, settings) {
  n <- nchar(result$raw_seq)
  plot(c(1,n), c(0,1), type="n", xlab="Base position", ylab="", yaxt="n", ylim=c(0,1),
       main=paste0(result$sample_id, " — read / trim overview"))
  segments(1,0.78,n,0.78,lwd=5,col="grey70")
  text(1,0.88,"RAW READ",adj=0,cex=0.8)
  sm <- result$summary
  if (!is.na(sm$trim_start) && !is.na(sm$trim_end)) {
    segments(sm$trim_start,0.50,sm$trim_end,0.50,lwd=8,col="steelblue")
    text(sm$trim_start,0.61,"TRIMMED",adj=0,cex=0.8,col="steelblue4")
  }

  if (isTRUE(settings$enable_primer_mapping)) {
    pm <- primer_match_table(result, settings)
    if (nrow(pm)) {
      for (i in seq_len(nrow(pm))) {
        if (is.na(pm$Start[i])) next
        y <- if (pm$Primer[i]=="Forward") 0.28 else 0.15
        col <- if (pm$Primer[i]=="Forward") "darkorange" else "purple"
        segments(pm$Start[i], y, pm$End[i], y, lwd=8, col=col)
        text(pm$Start[i], y+0.065, paste0(pm$Primer[i], " (", pm$Identity_percent[i], "%)"), adj=0, cex=0.70, col=col)
      }
    }
  }

}

make_sequence_preview <- function(result) {
  result <- ensure_curation_state(result)
  data.frame(
    Metric=c("Raw length","Curated length","GC % (curated)","N count (curated)","Current trim start","Current trim end","Automatic trim","Manual base edits","Curation revision"),
    Value=c(nchar(result$raw_seq),nchar(result$seq),sequence_gc(result$seq),sequence_n_count(result$seq),
            result$summary$trim_start,result$summary$trim_end,
            paste0(result$curation$auto_trim_start, "-", result$curation$auto_trim_end),
            nrow(result$curation$base_edits), result$curation$revision),
    stringsAsFactors=FALSE
  )
}


# ============================================================
# Ambiguous-peak / channel-competition review
# ============================================================

ambiguous_peak_defaults <- function() {
  list(
    strong_ratio = 1.25,
    moderate_ratio = 1.75,
    min_relative_signal = 0.20,
    max_trace_radius = 5L,
    auto_min_alt_to_called = 1.80,
    auto_min_alt_to_third = 2.00,
    auto_max_peak_offset = 2L,
    auto_min_relative_signal = 0.50
  )
}

ambiguous_peak_params_from_settings <- function(settings = NULL) {
  p <- ambiguous_peak_defaults()
  if (!is.null(settings) && is.list(settings)) {
    if (!is.null(settings$ambiguous_peak_strong_ratio) && is.finite(as.numeric(settings$ambiguous_peak_strong_ratio))) p$strong_ratio <- as.numeric(settings$ambiguous_peak_strong_ratio)
    if (!is.null(settings$ambiguous_peak_moderate_ratio) && is.finite(as.numeric(settings$ambiguous_peak_moderate_ratio))) p$moderate_ratio <- as.numeric(settings$ambiguous_peak_moderate_ratio)
    if (!is.null(settings$ambiguous_peak_min_relative_signal) && is.finite(as.numeric(settings$ambiguous_peak_min_relative_signal))) p$min_relative_signal <- as.numeric(settings$ambiguous_peak_min_relative_signal)
    if (!is.null(settings$auto_correct_min_alt_to_called) && is.finite(as.numeric(settings$auto_correct_min_alt_to_called))) p$auto_min_alt_to_called <- as.numeric(settings$auto_correct_min_alt_to_called)
    if (!is.null(settings$auto_correct_min_alt_to_third) && is.finite(as.numeric(settings$auto_correct_min_alt_to_third))) p$auto_min_alt_to_third <- as.numeric(settings$auto_correct_min_alt_to_third)
    if (!is.null(settings$auto_correct_max_peak_offset) && is.finite(as.numeric(settings$auto_correct_max_peak_offset))) p$auto_max_peak_offset <- as.integer(settings$auto_correct_max_peak_offset)
    if (!is.null(settings$auto_correct_min_relative_signal) && is.finite(as.numeric(settings$auto_correct_min_relative_signal))) p$auto_min_relative_signal <- as.numeric(settings$auto_correct_min_relative_signal)
  }
  p
}

# ============================================================
# Manual sequence curation state
# ============================================================

empty_curation_log <- function() {
  data.frame(
    Timestamp = character(), Transaction_ID = character(), Revision = integer(),
    Action = character(), Position = integer(), Before = character(), After = character(),
    Method = character(), Evidence = character(), Details = character(),
    stringsAsFactors = FALSE
  )
}

empty_base_edits <- function() {
  data.frame(Position = integer(), Original = character(), Edited = character(), stringsAsFactors = FALSE)
}

ensure_curation_state <- function(result) {
  if (is.null(result) || !is.list(result)) return(result)
  sm <- result$summary
  auto_start <- suppressWarnings(as.integer(if (!is.null(sm$trim_start)) sm$trim_start[1] else NA_integer_))
  auto_end <- suppressWarnings(as.integer(if (!is.null(sm$trim_end)) sm$trim_end[1] else NA_integer_))
  if (is.null(result$curation) || !is.list(result$curation)) {
    result$curation <- list(
      auto_trim_start = auto_start,
      auto_trim_end = auto_end,
      auto_seq = if (!is.null(result$seq)) as.character(result$seq) else "",
      auto_status = if (!is.null(sm$status)) as.character(sm$status[1]) else "OK",
      trim_start = auto_start,
      trim_end = auto_end,
      base_edits = empty_base_edits(),
      reviewed_positions = integer(),
      revision = 0L,
      audit_log = empty_curation_log(),
      undo_stack = list(),
      redo_stack = list()
    )
  } else {
    if (is.null(result$curation$auto_trim_start)) result$curation$auto_trim_start <- auto_start
    if (is.null(result$curation$auto_trim_end)) result$curation$auto_trim_end <- auto_end
    if (is.null(result$curation$auto_seq)) result$curation$auto_seq <- if (!is.null(result$seq)) as.character(result$seq) else ""
    if (is.null(result$curation$auto_status)) result$curation$auto_status <- if (!is.null(sm$status)) as.character(sm$status[1]) else "OK"
    if (is.null(result$curation$trim_start)) result$curation$trim_start <- auto_start
    if (is.null(result$curation$trim_end)) result$curation$trim_end <- auto_end
    if (!is.data.frame(result$curation$base_edits)) result$curation$base_edits <- empty_base_edits()
    if (is.null(result$curation$reviewed_positions)) result$curation$reviewed_positions <- integer()
    if (is.null(result$curation$revision)) result$curation$revision <- 0L
    if (!is.data.frame(result$curation$audit_log)) result$curation$audit_log <- empty_curation_log()
    if (!is.list(result$curation$undo_stack)) result$curation$undo_stack <- list()
    if (!is.list(result$curation$redo_stack)) result$curation$redo_stack <- list()
  }
  result
}

curation_snapshot <- function(result) {
  result <- ensure_curation_state(result)
  list(
    trim_start = result$curation$trim_start,
    trim_end = result$curation$trim_end,
    base_edits = result$curation$base_edits,
    reviewed_positions = as.integer(result$curation$reviewed_positions)
  )
}

curated_raw_calls <- function(result) {
  result <- ensure_curation_state(result)
  chars <- strsplit(sanitize_dna(result$raw_seq), "")[[1]]
  edits <- result$curation$base_edits
  if (is.data.frame(edits) && nrow(edits)) {
    for (i in seq_len(nrow(edits))) {
      pos <- suppressWarnings(as.integer(edits$Position[i]))
      if (is.finite(pos) && pos >= 1 && pos <= length(chars)) chars[pos] <- toupper(as.character(edits$Edited[i]))
    }
  }
  chars
}

curation_rebuild <- function(result, settings = NULL) {
  result <- ensure_curation_state(result)
  chars <- curated_raw_calls(result)
  n <- length(chars)
  start <- suppressWarnings(as.integer(result$curation$trim_start))
  end <- suppressWarnings(as.integer(result$curation$trim_end))
  if (!is.finite(start) || !is.finite(end) || start < 1 || end < start || start > n) {
    seq_out <- ""
    current_len <- 0L
  } else {
    end <- min(end, n)
    result$curation$trim_end <- end
    seq_out <- paste0(chars[start:end], collapse = "")
    current_len <- nchar(seq_out)
  }
  result$seq <- seq_out
  if (!is.null(result$summary) && nrow(result$summary)) {
    if (!"auto_trim_start" %in% names(result$summary)) result$summary$auto_trim_start <- result$curation$auto_trim_start
    if (!"auto_trim_end" %in% names(result$summary)) result$summary$auto_trim_end <- result$curation$auto_trim_end
    if (!"auto_trimmed_length" %in% names(result$summary)) result$summary$auto_trimmed_length <- nchar(result$curation$auto_seq)
    result$summary$trim_start <- if (current_len) start else NA_integer_
    result$summary$trim_end <- if (current_len) end else NA_integer_
    result$summary$trimmed_length <- current_len
    result$summary$manual_curation <- if (result$curation$revision > 0L) "Yes" else "No"
    result$summary$curation_revision <- as.integer(result$curation$revision)
    result$summary$manual_base_edits <- if (is.data.frame(result$curation$base_edits)) nrow(result$curation$base_edits) else 0L
    if (!is.null(settings) && !is.null(settings$min_usable_len) && current_len > 0 && current_len < as.integer(settings$min_usable_len)) {
      result$summary$status <- "SHORT_AFTER_TRIMMING"
    } else if (current_len == 0L) {
      result$summary$status <- "FAILED_TRIMMING"
    } else {
      result$summary$status <- if (identical(result$curation$auto_status, "OK")) "OK" else result$curation$auto_status
    }
  }
  result
}

curation_restore_snapshot <- function(result, snapshot, settings = NULL) {
  result <- ensure_curation_state(result)
  result$curation$trim_start <- snapshot$trim_start
  result$curation$trim_end <- snapshot$trim_end
  result$curation$base_edits <- if (is.data.frame(snapshot$base_edits)) snapshot$base_edits else empty_base_edits()
  result$curation$reviewed_positions <- if (!is.null(snapshot$reviewed_positions)) as.integer(snapshot$reviewed_positions) else integer()
  curation_rebuild(result, settings)
}

curation_new_transaction_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "-", sprintf("%06d", sample.int(999999, 1)))
}

curation_commit <- function(result, new_snapshot, action_rows, settings = NULL, label = "Manual curation") {
  result <- ensure_curation_state(result)
  before <- curation_snapshot(result)
  txid <- curation_new_transaction_id()
  old_revision <- as.integer(result$curation$revision)
  result <- curation_restore_snapshot(result, new_snapshot, settings)
  result$curation$revision <- old_revision + 1L
  result <- curation_rebuild(result, settings)
  after <- curation_snapshot(result)

  if (is.null(action_rows) || !is.data.frame(action_rows) || !nrow(action_rows)) {
    action_rows <- data.frame(Action = label, Position = NA_integer_, Before = "", After = "", Method = "Manual", Evidence = "", Details = "", stringsAsFactors = FALSE)
  }
  required <- c("Action","Position","Before","After","Method","Evidence","Details")
  for (nm in required) if (!nm %in% names(action_rows)) action_rows[[nm]] <- if (nm == "Position") NA_integer_ else ""
  action_rows <- action_rows[, required, drop = FALSE]
  action_rows$Timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  action_rows$Transaction_ID <- txid
  action_rows$Revision <- as.integer(result$curation$revision)
  action_rows <- action_rows[, c("Timestamp","Transaction_ID","Revision", required), drop = FALSE]
  result$curation$audit_log <- rbind(result$curation$audit_log, action_rows)

  result$curation$undo_stack[[length(result$curation$undo_stack) + 1L]] <- list(
    id = txid, label = label, before = before, after = after
  )
  result$curation$redo_stack <- list()
  result
}

curation_undo <- function(result, settings = NULL) {
  result <- ensure_curation_state(result)
  if (!length(result$curation$undo_stack)) return(list(result = result, changed = FALSE, label = ""))
  tx <- result$curation$undo_stack[[length(result$curation$undo_stack)]]
  result$curation$undo_stack <- head(result$curation$undo_stack, -1)
  result <- curation_restore_snapshot(result, tx$before, settings)
  result$curation$revision <- as.integer(result$curation$revision) + 1L
  result <- curation_rebuild(result, settings)
  result$curation$redo_stack[[length(result$curation$redo_stack) + 1L]] <- tx
  row <- data.frame(
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Transaction_ID = tx$id,
    Revision = as.integer(result$curation$revision), Action = "UNDO", Position = NA_integer_,
    Before = "", After = "", Method = "Undo", Evidence = "", Details = tx$label,
    stringsAsFactors = FALSE
  )
  result$curation$audit_log <- rbind(result$curation$audit_log, row)
  list(result = result, changed = TRUE, label = tx$label)
}

curation_redo <- function(result, settings = NULL) {
  result <- ensure_curation_state(result)
  if (!length(result$curation$redo_stack)) return(list(result = result, changed = FALSE, label = ""))
  tx <- result$curation$redo_stack[[length(result$curation$redo_stack)]]
  result$curation$redo_stack <- head(result$curation$redo_stack, -1)
  result <- curation_restore_snapshot(result, tx$after, settings)
  result$curation$revision <- as.integer(result$curation$revision) + 1L
  result <- curation_rebuild(result, settings)
  result$curation$undo_stack[[length(result$curation$undo_stack) + 1L]] <- tx
  row <- data.frame(
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Transaction_ID = tx$id,
    Revision = as.integer(result$curation$revision), Action = "REDO", Position = NA_integer_,
    Before = "", After = "", Method = "Redo", Evidence = "", Details = tx$label,
    stringsAsFactors = FALSE
  )
  result$curation$audit_log <- rbind(result$curation$audit_log, row)
  list(result = result, changed = TRUE, label = tx$label)
}

curation_set_base_snapshot <- function(result, position, new_base) {
  result <- ensure_curation_state(result)
  snap <- curation_snapshot(result)
  pos <- as.integer(position)
  new_base <- toupper(as.character(new_base))
  raw_chars <- strsplit(sanitize_dna(result$raw_seq), "")[[1]]
  if (!is.finite(pos) || pos < 1 || pos > length(raw_chars) || !new_base %in% c("A","C","G","T","N","R","Y","S","W","K","M")) return(NULL)
  original <- raw_chars[pos]
  edits <- snap$base_edits
  edits <- edits[edits$Position != pos, , drop = FALSE]
  if (!identical(new_base, original)) {
    edits <- rbind(edits, data.frame(Position = pos, Original = original, Edited = new_base, stringsAsFactors = FALSE))
    edits <- edits[order(edits$Position), , drop = FALSE]
  }
  snap$base_edits <- edits
  snap$reviewed_positions <- sort(unique(c(snap$reviewed_positions, pos)))
  snap
}

curation_trim_snapshot <- function(result, position, side = c("left", "right")) {
  side <- match.arg(side)
  result <- ensure_curation_state(result)
  snap <- curation_snapshot(result)
  pos <- as.integer(position)
  start <- as.integer(snap$trim_start); end <- as.integer(snap$trim_end)
  if (!is.finite(pos) || !is.finite(start) || !is.finite(end) || pos < start || pos > end) return(NULL)
  if (side == "left") snap$trim_start <- pos + 1L else snap$trim_end <- pos - 1L
  if (!is.finite(snap$trim_start) || !is.finite(snap$trim_end) || snap$trim_end < snap$trim_start) return(NULL)
  snap$reviewed_positions <- sort(unique(c(snap$reviewed_positions, pos)))
  snap
}

curation_review_snapshot <- function(result, position) {
  result <- ensure_curation_state(result)
  snap <- curation_snapshot(result)
  snap$reviewed_positions <- sort(unique(c(snap$reviewed_positions, as.integer(position))))
  snap
}

curation_reset_snapshot <- function(result) {
  result <- ensure_curation_state(result)
  list(
    trim_start = result$curation$auto_trim_start,
    trim_end = result$curation$auto_trim_end,
    base_edits = empty_base_edits(),
    reviewed_positions = integer()
  )
}

iupac_for_pair <- function(a, b) {
  key <- paste(sort(unique(c(toupper(a), toupper(b)))), collapse = "")
  map <- c(AG="R", CT="Y", CG="S", AT="W", GT="K", AC="M")
  if (key %in% names(map)) unname(map[[key]]) else NA_character_
}

collect_curation_log <- function(results) {
  empty <- data.frame(
    Sample = character(), Timestamp = character(), Transaction_ID = character(), Revision = integer(),
    Action = character(), Position = integer(), Before = character(), After = character(), Method = character(),
    Evidence = character(), Details = character(), stringsAsFactors = FALSE
  )
  if (is.null(results) || !length(results)) return(empty)
  parts <- lapply(names(results), function(nm) {
    r <- ensure_curation_state(results[[nm]])
    lg <- r$curation$audit_log
    if (!is.data.frame(lg) || !nrow(lg)) return(NULL)
    data.frame(Sample = nm, lg, stringsAsFactors = FALSE)
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(empty)
  out <- do.call(rbind, parts); rownames(out) <- NULL; out
}

# Detect individual base positions where another dye channel competes strongly
# with the current curated base call. This is a REVIEW aid, not a biological
# diagnosis. Low-signal noise is suppressed using the retained-region median.
ambiguous_peak_flags <- function(result, scope = c("trimmed", "raw"), params = ambiguous_peak_defaults()) {
  scope <- match.arg(scope)
  result <- ensure_curation_state(result)

  trace <- result$trace
  peak_pos <- as.integer(result$peak_pos)
  cmap <- result$channel_map
  bases <- curated_raw_calls(result)
  n <- min(length(bases), length(peak_pos))

  empty <- data.frame(
    Position = integer(), Call = character(), Competing_channel = character(),
    Called_signal = numeric(), Competitor_signal = numeric(), Peak_ratio = numeric(),
    Competitor_percent = numeric(), Competitor_peak_offset = numeric(),
    Alternative_to_called_ratio = numeric(), Alternative_to_third_ratio = numeric(),
    Severity = character(), Flag = character(), Auto_correct_candidate = logical(),
    stringsAsFactors = FALSE
  )

  if (is.null(trace) || !nrow(trace) || n < 1 || is.null(cmap)) return(empty)

  sm <- result$summary
  trim_start <- suppressWarnings(as.integer(if (!is.null(sm$trim_start)) sm$trim_start[1] else NA_integer_))
  trim_end <- suppressWarnings(as.integer(if (!is.null(sm$trim_end)) sm$trim_end[1] else NA_integer_))

  if (scope == "trimmed" && is.finite(trim_start) && is.finite(trim_end) && trim_end >= trim_start) {
    positions <- seq.int(max(1L, trim_start), min(n, trim_end))
  } else {
    positions <- seq_len(n)
  }

  baseline_idx <- if (is.finite(trim_start) && is.finite(trim_end) && trim_end >= trim_start) {
    seq.int(max(1L, trim_start), min(n, trim_end))
  } else seq_len(n)
  metric_signal <- if (!is.null(result$metrics) && !is.null(result$metrics$called_signal)) result$metrics$called_signal else rep(NA_real_, n)
  baseline <- suppressWarnings(stats::median(metric_signal[baseline_idx], na.rm = TRUE))
  if (!is.finite(baseline) || baseline <= 0) baseline <- suppressWarnings(stats::median(metric_signal, na.rm = TRUE))
  if (!is.finite(baseline) || baseline <= 0) baseline <- 1
  min_signal <- baseline * as.numeric(params$min_relative_signal)
  auto_min_signal <- baseline * as.numeric(params$auto_min_relative_signal)

  channel_names <- c("A", "C", "G", "T")
  channel_cols <- unname(cmap[channel_names])
  valid_map <- is.finite(channel_cols) & channel_cols >= 1 & channel_cols <= ncol(trace)
  if (!all(valid_map)) return(empty)

  rows <- list()
  for (i in positions) {
    pos <- peak_pos[i]
    if (!is.finite(pos) || pos < 1 || pos > nrow(trace)) next
    gaps <- numeric()
    if (i > 1 && is.finite(peak_pos[i-1]) && peak_pos[i] > peak_pos[i-1]) gaps <- c(gaps, peak_pos[i] - peak_pos[i-1])
    if (i < n && is.finite(peak_pos[i+1]) && peak_pos[i+1] > peak_pos[i]) gaps <- c(gaps, peak_pos[i+1] - peak_pos[i])
    radius <- if (length(gaps)) max(1L, floor(min(gaps) / 3)) else 2L
    radius <- min(as.integer(params$max_trace_radius), radius)
    idx <- seq.int(max(1L, pos - radius), min(nrow(trace), pos + radius))

    local <- numeric(length(channel_names)); local_idx <- integer(length(channel_names))
    names(local) <- names(local_idx) <- channel_names
    for (j in seq_along(channel_names)) {
      v <- trace[idx, channel_cols[j]]
      if (!any(is.finite(v))) { local[j] <- NA_real_; local_idx[j] <- NA_integer_; next }
      k <- which.max(replace(v, !is.finite(v), -Inf))[1]
      local[j] <- v[k]
      local_idx[j] <- idx[k]
    }
    if (!any(is.finite(local))) next

    call <- bases[i]
    strongest_channel <- names(which.max(replace(local, !is.finite(local), -Inf)))[1]
    strongest_signal <- local[[strongest_channel]]
    is_n <- identical(call, "N")
    if (!is_n && (!is.finite(strongest_signal) || strongest_signal < min_signal)) next

    competitor_peak_offset <- NA_real_
    alt_to_called <- NA_real_
    alt_to_third <- NA_real_
    auto_candidate <- FALSE

    if (call %in% channel_names) {
      called_signal <- local[[call]]
      competitors <- local[setdiff(channel_names, call)]
      competitor_channel <- names(which.max(replace(competitors, !is.finite(competitors), -Inf)))[1]
      competitor_signal <- competitors[[competitor_channel]]
      if (!is.finite(called_signal)) called_signal <- 0
      if (!is.finite(competitor_signal)) competitor_signal <- 0
      ratio <- if (competitor_signal > 0) called_signal / competitor_signal else Inf
      competitor_pct <- if (called_signal > 0) 100 * competitor_signal / called_signal else Inf
      competitor_peak_offset <- if (is.finite(local_idx[[competitor_channel]])) local_idx[[competitor_channel]] - pos else NA_real_
      alt_to_called <- if (called_signal > 0) competitor_signal / called_signal else Inf
      third_channels <- setdiff(channel_names, c(call, competitor_channel))
      third_signal <- suppressWarnings(max(local[third_channels], na.rm = TRUE))
      if (!is.finite(third_signal)) third_signal <- 0
      alt_to_third <- if (third_signal > 0) competitor_signal / third_signal else Inf

      if (!identical(strongest_channel, call)) {
        severity <- "Strong"; flag <- "Called base not dominant"
      } else if (is.finite(ratio) && ratio <= as.numeric(params$strong_ratio)) {
        severity <- "Strong"; flag <- "Strong channel competition"
      } else if (is.finite(ratio) && ratio <= as.numeric(params$moderate_ratio)) {
        severity <- "Moderate"; flag <- "Moderate channel competition"
      } else next

      auto_candidate <- identical(flag, "Called base not dominant") &&
        identical(competitor_channel, strongest_channel) &&
        is.finite(alt_to_called) && alt_to_called >= as.numeric(params$auto_min_alt_to_called) &&
        is.finite(alt_to_third) && alt_to_third >= as.numeric(params$auto_min_alt_to_third) &&
        is.finite(competitor_peak_offset) && abs(competitor_peak_offset) <= as.integer(params$auto_max_peak_offset) &&
        is.finite(competitor_signal) && competitor_signal >= auto_min_signal
    } else {
      ord <- order(local, decreasing = TRUE, na.last = NA)
      if (!length(ord)) next
      first <- ord[1]; second <- if (length(ord) >= 2) ord[2] else ord[1]
      called_signal <- local[first]; competitor_signal <- local[second]
      competitor_channel <- paste0(names(local)[first], "/", names(local)[second])
      ratio <- if (is.finite(competitor_signal) && competitor_signal > 0) called_signal / competitor_signal else Inf
      competitor_pct <- if (is.finite(called_signal) && called_signal > 0) 100 * competitor_signal / called_signal else NA_real_
      severity <- if (is.finite(strongest_signal) && strongest_signal >= min_signal) "Strong" else "Review"
      flag <- if (identical(call, "N")) "Ambiguous N base call" else "IUPAC ambiguity call"
    }

    rows[[length(rows) + 1L]] <- data.frame(
      Position = i, Call = call, Competing_channel = competitor_channel,
      Called_signal = round(called_signal, 1), Competitor_signal = round(competitor_signal, 1),
      Peak_ratio = if (is.finite(ratio)) round(ratio, 2) else Inf,
      Competitor_percent = if (is.finite(competitor_pct)) round(competitor_pct, 1) else NA_real_,
      Competitor_peak_offset = if (is.finite(competitor_peak_offset)) round(competitor_peak_offset, 1) else NA_real_,
      Alternative_to_called_ratio = if (is.finite(alt_to_called)) round(alt_to_called, 2) else alt_to_called,
      Alternative_to_third_ratio = if (is.finite(alt_to_third)) round(alt_to_third, 2) else alt_to_third,
      Severity = severity, Flag = flag, Auto_correct_candidate = isTRUE(auto_candidate),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out[order(out$Position), , drop = FALSE]
}

high_confidence_autocorrections <- function(result, settings = NULL) {
  flags <- ambiguous_peak_flags(result, scope = "trimmed", params = ambiguous_peak_params_from_settings(settings))
  if (!nrow(flags)) return(flags[0, , drop = FALSE])
  flags[which(flags$Auto_correct_candidate %in% TRUE), , drop = FALSE]
}

collect_ambiguous_peak_flags <- function(results, scope = "trimmed", settings = NULL) {
  if (is.null(results) || !length(results)) return(data.frame())
  parts <- lapply(names(results), function(sample_name) {
    result <- results[[sample_name]]
    if (is.null(result)) return(NULL)
    result <- ensure_curation_state(result)
    df <- ambiguous_peak_flags(result, scope = scope, params = ambiguous_peak_params_from_settings(settings))
    if (!nrow(df)) return(NULL)
    df$Review_status <- ifelse(df$Position %in% as.integer(result$curation$reviewed_positions), "Reviewed", "Active")
    data.frame(Sample = sample_name, df, stringsAsFactors = FALSE)
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    return(data.frame(
      Sample=character(), Position=integer(), Call=character(), Competing_channel=character(),
      Called_signal=numeric(), Competitor_signal=numeric(), Peak_ratio=numeric(),
      Competitor_percent=numeric(), Competitor_peak_offset=numeric(),
      Alternative_to_called_ratio=numeric(), Alternative_to_third_ratio=numeric(),
      Severity=character(), Flag=character(), Auto_correct_candidate=logical(), Review_status=character(),
      stringsAsFactors=FALSE
    ))
  }
  out <- do.call(rbind, parts); rownames(out) <- NULL; out
}


# ============================================================
# Interactive chromatogram (client-side Plotly)
# ============================================================

# Build the complete chromatogram once and let Plotly perform zoom/pan/range-slider
# interactions in the browser. This avoids a Shiny server round-trip for every
# horizontal movement and keeps the X axis in called-base coordinates.
make_chromatogram_plotly <- function(result, settings, flags = NULL, show_flags = TRUE) {
  trace <- result$trace
  peak_pos <- as.integer(result$peak_pos)
  cmap <- result$channel_map
  n <- length(peak_pos)

  if (is.null(trace) || !nrow(trace) || n < 2) {
    return(plotly::plot_ly() |> plotly::layout(title = "No chromatogram data available"))
  }

  valid_peaks <- which(is.finite(peak_pos) & peak_pos >= 1 & peak_pos <= nrow(trace))
  if (length(valid_peaks) < 2) {
    return(plotly::plot_ly() |> plotly::layout(title = "Insufficient peak-position data"))
  }

  trace_idx <- seq.int(min(peak_pos[valid_peaks]), max(peak_pos[valid_peaks]))
  x_base <- trace_to_base_x(peak_pos, trace_idx)

  channel_values <- lapply(c("A", "C", "G", "T"), function(base) {
    col <- cmap[[base]]
    if (is.null(col) || is.na(col)) rep(NA_real_, length(trace_idx)) else trace[trace_idx, col]
  })
  names(channel_values) <- c("A", "C", "G", "T")

  all_vals <- unlist(channel_values, use.names = FALSE)
  all_vals <- all_vals[is.finite(all_vals)]
  ymax <- if (length(all_vals)) as.numeric(stats::quantile(all_vals, 0.995, na.rm = TRUE)) else 1
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(all_vals, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  # Use conventional chromatogram colors. The lines are WebGL traces, so zooming
  # and panning occur in the browser without repeatedly rendering PNG images.
  p <- plotly::plot_ly()
  colors <- c(A = "#2ca02c", C = "#1f77b4", G = "#111111", T = "#d62728")
  for (base in c("A", "C", "G", "T")) {
    p <- plotly::add_trace(
      p,
      x = x_base,
      y = channel_values[[base]],
      type = "scattergl",
      mode = "lines",
      name = base,
      line = list(color = colors[[base]], width = 1.2),
      hovertemplate = paste0(base, "<br>Base position: %{x:.2f}<br>Signal: %{y:.0f}<extra></extra>")
    )
  }

  # Base calls are separate text labels anchored to integer base positions. They
  # automatically move with the Plotly viewport and become readable as the user
  # zooms in. At a wide overview they may overlap by design; the range slider is
  # intended for navigation rather than detailed base reading.
  bases <- curated_raw_calls(result)
  # Reserve separate vertical bands for called bases and QC markers so they do
  # not collide with each other or with the legend.
  label_y <- rep(ymax * 1.055, n)
  p <- plotly::add_trace(
    p,
    x = seq_len(n),
    y = label_y,
    type = "scatter",
    mode = "text",
    text = bases,
    textfont = list(size = 10, color = "#374151"),
    hoverinfo = "skip",
    showlegend = FALSE,
    cliponaxis = FALSE
  )

  shapes <- list()
  annotations <- list()
  sm <- result$summary
  if (!is.na(sm$trim_start)) {
    shapes[[length(shapes) + 1L]] <- list(
      type = "line", x0 = sm$trim_start, x1 = sm$trim_start,
      y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#6b7280", dash = "dash", width = 1.5)
    )
  }
  if (!is.na(sm$trim_end)) {
    shapes[[length(shapes) + 1L]] <- list(
      type = "line", x0 = sm$trim_end, x1 = sm$trim_end,
      y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#6b7280", dash = "dash", width = 1.5)
    )
  }

  if (isTRUE(settings$enable_primer_mapping)) {
    pm <- primer_match_table(result, settings)
    if (nrow(pm)) {
      for (i in seq_len(nrow(pm))) {
        if (is.na(pm$Start[i]) || is.na(pm$End[i])) next
        col <- if (pm$Primer[i] == "Forward") "#d97706" else "#7e22ce"
        shapes[[length(shapes) + 1L]] <- list(
          type = "rect", x0 = pm$Start[i], x1 = pm$End[i],
          y0 = 1.045, y1 = 1.09, yref = "paper",
          fillcolor = col, opacity = 0.22,
          line = list(color = col, width = 1)
        )
        annotations[[length(annotations) + 1L]] <- list(
          x = mean(c(pm$Start[i], pm$End[i])), y = 1.095, yref = "paper",
          text = paste0(pm$Primer[i], " primer"), showarrow = FALSE,
          font = list(size = 10, color = col)
        )
      }
    }
  }

  # Ambiguous-peak flags are review markers only. The table below the plot is
  # the main inspection surface; these markers simply make suspicious positions
  # visible while navigating the chromatogram.
  if (isTRUE(show_flags) && !is.null(flags) && is.data.frame(flags) && nrow(flags)) {
    flag_colors <- ifelse(flags$Severity == "Strong", "#dc2626",
                          ifelse(flags$Severity == "Moderate", "#d97706", "#7c3aed"))
    hover <- paste0(
      "Base ", flags$Position,
      "<br>Call: ", flags$Call,
      "<br>Competitor: ", flags$Competing_channel,
      "<br>Peak ratio: ", flags$Peak_ratio,
      "<br>Competitor / called: ", flags$Competitor_percent, "%",
      "<br>", flags$Flag,
      ifelse(!is.null(flags$Auto_correct_candidate) & flags$Auto_correct_candidate, "<br>High-confidence auto-correction candidate", "")
    )
    p <- plotly::add_trace(
      p,
      x = flags$Position,
      y = rep(ymax * 1.155, nrow(flags)),
      type = "scatter",
      mode = "markers",
      name = "QC flags",
      marker = list(symbol = "triangle-down", size = 10, color = flag_colors,
                    line = list(color = "#ffffff", width = 0.5)),
      text = hover,
      hovertemplate = "%{text}<extra></extra>",
      showlegend = TRUE,
      cliponaxis = FALSE
    )
  }

  # Start at a useful detail range, but keep the complete read available in the
  # native Plotly range slider directly beneath the graph.
  initial_start <- if (!is.na(sm$trim_start)) max(1, sm$trim_start - 10) else 1
  initial_end <- min(n, initial_start + 79)

  p <- plotly::layout(
    p,
    # Keep the internal title to the sample ID only. The previous em dash was
    # rendered as byte markers on some Windows/R locale combinations.
    title = list(text = as.character(result$sample_id)[1], x = 0),
    dragmode = "pan",
    hovermode = "x unified",
    margin = list(l = 65, r = 125, t = 75, b = 55),
    legend = list(
      orientation = "v", x = 1.02, xanchor = "left",
      y = 1, yanchor = "top",
      bgcolor = "rgba(255,255,255,0.92)",
      bordercolor = "#e5e7eb", borderwidth = 1
    ),
    xaxis = list(
      title = "Base position",
      range = c(initial_start, initial_end),
      rangeslider = list(visible = TRUE, thickness = 0.10)
    ),
    yaxis = list(title = "Signal", range = c(0, ymax * 1.22), fixedrange = FALSE),
    shapes = shapes,
    annotations = annotations
  )

  plotly::config(
    p,
    scrollZoom = TRUE,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("select2d", "lasso2d")
  )
}


# ============================================================
# QC plots used both on-screen and in exported checkpoints
# ============================================================

draw_qc_metrics <- function(result, settings) {
  metrics <- result$metrics
  sm <- result$summary
  if (is.null(metrics) || !nrow(metrics)) {
    plot.new(); text(0.5, 0.5, "No QC metrics available")
    return(invisible(NULL))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 5, 3, 2))

  plot(
    metrics$index, metrics$called_signal, type = "l",
    xlab = "Base-call index", ylab = "Called-base signal",
    main = paste0(result$sample_id, " — called-base signal")
  )
  if (!is.na(sm$trim_start)) abline(v = sm$trim_start, lty = 2)
  if (!is.na(sm$trim_end)) abline(v = sm$trim_end, lty = 2)
  if (!is.null(settings$absolute_max_base_index)) abline(v = settings$absolute_max_base_index, lty = 3)
  if (!is.na(sm$collapse_index)) abline(v = sm$collapse_index, lty = 4)

  finite_ratio <- metrics$peak_ratio[is.finite(metrics$peak_ratio)]
  ymax <- if (length(finite_ratio)) max(5, min(10, max(finite_ratio, na.rm = TRUE))) else 5
  plot(
    metrics$index, metrics$peak_ratio, type = "l",
    xlab = "Base-call index", ylab = "Called peak / second peak",
    main = paste0(result$sample_id, " — peak ratio"), ylim = c(0, ymax)
  )
  if (!is.null(settings$min_peak_ratio)) abline(h = settings$min_peak_ratio, lty = 2)
  if (!is.na(sm$trim_start)) abline(v = sm$trim_start, lty = 2)
  if (!is.na(sm$trim_end)) abline(v = sm$trim_end, lty = 2)
  if (!is.null(settings$absolute_max_base_index)) abline(v = settings$absolute_max_base_index, lty = 3)
  if (!is.na(sm$collapse_index)) abline(v = sm$collapse_index, lty = 4)
  invisible(NULL)
}

write_qc_plot_png <- function(result, settings, file) {
  grDevices::png(file, width = 1500, height = 900, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_qc_metrics(result, settings)
  invisible(file)
}
