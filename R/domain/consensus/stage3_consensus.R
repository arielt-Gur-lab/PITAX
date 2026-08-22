# ============================================================
# PITAX v3.0.0-alpha.8
# Stage 3: auditable Forward / Reverse consensus review and curation
# ============================================================

stage3_scalar_text <- function(x, fallback = "") {
  if (is.null(x) || !length(x) || is.na(x[1])) return(fallback)
  value <- trimws(as.character(x[1]))
  if (nzchar(value)) value else fallback
}

stage3_sanitize_dna <- function(x) {
  x <- toupper(gsub("[^A-Za-z]", "", stage3_scalar_text(x)))
  gsub("[^ACGTRYSWKMBDHVN]", "N", x)
}

stage3_reverse_complement <- function(x) {
  chars <- strsplit(stage3_sanitize_dna(x), "", fixed = TRUE)[[1]]
  if (!length(chars)) return("")
  map <- c(
    A = "T", C = "G", G = "C", T = "A",
    R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
    B = "V", V = "B", D = "H", H = "D", N = "N"
  )
  paste0(unname(map[rev(chars)]), collapse = "")
}

stage3_clean_name_part <- function(x) {
  x <- gsub("\\s+", "-", stage3_scalar_text(x), perl = TRUE)
  x <- gsub("[^A-Za-z0-9_.-]", "-", x, perl = TRUE)
  x <- gsub("-+", "-", x, perl = TRUE)
  gsub("(^[-.]+|[-.]+$)", "", x, perl = TRUE)
}

stage3_compose_sequence_name <- function(isolate, locus) {
  isolate <- stage3_clean_name_part(isolate)
  locus <- stage3_clean_name_part(locus)
  if (!nzchar(isolate) || !nzchar(locus)) return("")
  paste(isolate, locus, sep = "_")
}

stage3_empty_summary <- function() {
  data.frame(
    Consensus_ID = character(), Final_Name = character(), Isolate = character(), Locus = character(),
    Status = character(), Source_Reads = integer(), Forward_Read = character(), Reverse_Read = character(),
    Length = integer(), Overlap = integer(), Identity_percent = numeric(), Mismatches = integer(),
    Indels = integer(), Review_positions = integer(), Revision = integer(),
    Algorithm = character(), Built_at = character(),
    stringsAsFactors = FALSE
  )
}

stage3_empty_consensus_audit <- function() {
  data.frame(
    Revision = integer(), Timestamp = character(), Action = character(),
    Alignment_Column = integer(), Previous_Call = character(), New_Call = character(),
    Method = character(), Note = character(), stringsAsFactors = FALSE
  )
}

stage3_record_snapshot <- function(record) {
  list(
    status = record$status,
    sequence = record$sequence,
    evidence = record$evidence,
    metrics = record$metrics
  )
}

stage3_restore_record_snapshot <- function(record, snapshot) {
  record$status <- snapshot$status
  record$sequence <- snapshot$sequence
  record$evidence <- snapshot$evidence
  record$metrics <- snapshot$metrics
  record
}

stage3_ensure_record_curation <- function(record) {
  if (!is.list(record)) return(record)
  evidence <- if (is.data.frame(record$evidence)) record$evidence else data.frame()
  n <- nrow(evidence)
  if (!"Automatic_Call" %in% names(evidence)) evidence$Automatic_Call <- if (n) as.character(evidence$Consensus_Call) else character()
  if (!"Original_Needs_Review" %in% names(evidence)) evidence$Original_Needs_Review <- if (n) as.logical(evidence$Needs_Review) else logical()
  if (!"Manual_Review" %in% names(evidence)) evidence$Manual_Review <- rep(FALSE, n)
  if (!"Review_Note" %in% names(evidence)) evidence$Review_Note <- rep("", n)
  if (!"Review_Revision" %in% names(evidence)) evidence$Review_Revision <- rep(NA_integer_, n)
  record$evidence <- evidence

  cur <- if (is.list(record$curation)) record$curation else list()
  revision <- suppressWarnings(as.integer(cur$revision[1]))
  if (!length(revision) || !is.finite(revision[1])) revision <- 0L
  audit <- if (is.data.frame(cur$audit_log)) cur$audit_log else stage3_empty_consensus_audit()
  for (nm in names(stage3_empty_consensus_audit())) {
    if (!nm %in% names(audit)) {
      audit[[nm]] <- if (nm %in% c("Revision", "Alignment_Column")) rep(NA_integer_, nrow(audit)) else rep("", nrow(audit))
    }
  }
  audit <- audit[, names(stage3_empty_consensus_audit()), drop = FALSE]
  record$curation <- list(
    revision = revision,
    audit_log = audit,
    undo_stack = if (is.list(cur$undo_stack)) cur$undo_stack else list(),
    redo_stack = if (is.list(cur$redo_stack)) cur$redo_stack else list(),
    automatic_sequence = stage3_scalar_text(cur$automatic_sequence, stage3_scalar_text(record$sequence)),
    automatic_status = stage3_scalar_text(cur$automatic_status, stage3_scalar_text(record$status))
  )
  record
}

stage3_refresh_record <- function(record) {
  record <- stage3_ensure_record_curation(record)
  if (!is.data.frame(record$evidence) || !nrow(record$evidence)) return(record)
  review_count <- sum(as.logical(record$evidence$Needs_Review), na.rm = TRUE)
  record$metrics$review_positions <- as.integer(review_count)
  if (!stage3_scalar_text(record$status) %in% c("NO_RELIABLE_OVERLAP", "SOURCE_MISSING")) {
    record$sequence <- paste0(as.character(record$evidence$Consensus_Call), collapse = "")
    if (!is.null(record$alignment)) {
      record$status <- if (review_count) "REVIEW_REQUIRED" else "READY"
    }
  }
  record
}

stage3_append_consensus_audit <- function(record, action, alignment_column = NA_integer_,
                                          previous_call = "", new_call = "", method = "", note = "") {
  row <- data.frame(
    Revision = as.integer(record$curation$revision), Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Action = as.character(action), Alignment_Column = as.integer(alignment_column),
    Previous_Call = as.character(previous_call), New_Call = as.character(new_call),
    Method = as.character(method), Note = as.character(note), stringsAsFactors = FALSE
  )
  record$curation$audit_log <- rbind(record$curation$audit_log, row)
  record
}

stage3_apply_consensus_call <- function(record, alignment_column, call, method = "Manual review", note = "") {
  record <- stage3_ensure_record_curation(record)
  if (stage3_scalar_text(record$status) %in% c("NO_RELIABLE_OVERLAP", "SOURCE_MISSING")) {
    stop("This record cannot be curated until usable source evidence and a reliable overlap exist.", call. = FALSE)
  }
  column <- suppressWarnings(as.integer(alignment_column[1]))
  call <- stage3_sanitize_dna(call)
  if (!length(column) || !is.finite(column[1]) || nchar(call) != 1L) stop("Choose one valid IUPAC base call for a review position.", call. = FALSE)
  idx <- match(column, record$evidence$Alignment_Column)
  if (is.na(idx)) stop("The selected alignment column is not present in this consensus.", call. = FALSE)
  reviewable <- isTRUE(record$evidence$Original_Needs_Review[idx]) || isTRUE(record$evidence$Needs_Review[idx]) || isTRUE(record$evidence$Manual_Review[idx])
  if (!reviewable) stop("Only an evidence position flagged for consensus review can be edited here.", call. = FALSE)
  allowed <- unique(c(record$evidence$Forward_Base[idx], record$evidence$Reverse_Base[idx], record$evidence$Automatic_Call[idx]))
  allowed <- allowed[allowed != "-" & nzchar(allowed)]
  if (!call %in% allowed) stop("Choose the Forward, Reverse or automatic IUPAC call shown for this evidence position.", call. = FALSE)

  before <- stage3_record_snapshot(record)
  previous <- as.character(record$evidence$Consensus_Call[idx])
  record$curation$undo_stack <- c(record$curation$undo_stack, list(before))
  record$curation$redo_stack <- list()
  record$curation$revision <- as.integer(record$curation$revision) + 1L
  record$evidence$Consensus_Call[idx] <- call
  record$evidence$Needs_Review[idx] <- FALSE
  record$evidence$Manual_Review[idx] <- TRUE
  record$evidence$Review_Note[idx] <- stage3_scalar_text(note)
  record$evidence$Review_Revision[idx] <- record$curation$revision
  record$evidence$Decision[idx] <- paste0("Manual review | ", stage3_scalar_text(method, "Selected call"))
  record <- stage3_refresh_record(record)
  stage3_append_consensus_audit(record, "APPLY", column, previous, call, method, note)
}

stage3_consensus_undo <- function(record) {
  record <- stage3_ensure_record_curation(record)
  if (!length(record$curation$undo_stack)) stop("There is no consensus edit to undo.", call. = FALSE)
  current <- stage3_record_snapshot(record)
  target <- record$curation$undo_stack[[length(record$curation$undo_stack)]]
  cur <- record$curation
  cur$undo_stack <- cur$undo_stack[-length(cur$undo_stack)]
  cur$redo_stack <- c(cur$redo_stack, list(current))
  cur$revision <- as.integer(cur$revision) + 1L
  record <- stage3_restore_record_snapshot(record, target)
  record$curation <- cur
  record <- stage3_refresh_record(record)
  stage3_append_consensus_audit(record, "UNDO", NA_integer_, "", "", "Undo", "Restored the previous consensus state")
}

stage3_consensus_redo <- function(record) {
  record <- stage3_ensure_record_curation(record)
  if (!length(record$curation$redo_stack)) stop("There is no consensus edit to redo.", call. = FALSE)
  current <- stage3_record_snapshot(record)
  target <- record$curation$redo_stack[[length(record$curation$redo_stack)]]
  cur <- record$curation
  cur$redo_stack <- cur$redo_stack[-length(cur$redo_stack)]
  cur$undo_stack <- c(cur$undo_stack, list(current))
  cur$revision <- as.integer(cur$revision) + 1L
  record <- stage3_restore_record_snapshot(record, target)
  record$curation <- cur
  record <- stage3_refresh_record(record)
  stage3_append_consensus_audit(record, "REDO", NA_integer_, "", "", "Redo", "Restored the next consensus state")
}

stage3_reviewable_evidence <- function(record) {
  record <- stage3_ensure_record_curation(record)
  if (!nrow(record$evidence)) return(record$evidence)
  keep <- as.logical(record$evidence$Original_Needs_Review) | as.logical(record$evidence$Needs_Review) | as.logical(record$evidence$Manual_Review)
  record$evidence[which(keep), , drop = FALSE]
}

stage3_empty_consensus_set <- function() {
  list(
    schema = "pitax-consensus-set-v1",
    algorithm = "pitax-overlap-consensus-v1",
    built_at = "",
    settings = list(),
    records = list(),
    summary = stage3_empty_summary()
  )
}

stage3_run_locus_error <- function(assignments) {
  if (!is.data.frame(assignments) || !nrow(assignments) || !"Locus" %in% names(assignments)) {
    return("No assigned reads are available.")
  }
  loci <- unique(trimws(as.character(assignments$Locus)))
  loci <- loci[nzchar(loci)]
  if (length(loci) != 1L) {
    return("A sequencing run must contain exactly one gene/locus. Process different loci as separate runs; Stage 4 will combine their isolate-level profiles.")
  }
  NULL
}

stage3_result_evidence <- function(result, direction) {
  sequence <- if (is.list(result)) stage3_sanitize_dna(result$seq) else ""
  n <- nchar(sequence)
  if (!n) {
    return(list(sequence = "", quality = numeric(), peak_ratio = numeric(), raw_position = integer(), oriented_position = integer()))
  }

  start <- NA_integer_
  if (is.list(result$curation)) start <- suppressWarnings(as.integer(result$curation$trim_start[1]))
  if (!is.finite(start) && is.data.frame(result$summary) && nrow(result$summary)) {
    start <- suppressWarnings(as.integer(result$summary$trim_start[1]))
  }
  if (!is.finite(start)) start <- 1L
  raw_position <- seq.int(start, length.out = n)

  quality <- rep(NA_real_, n)
  detail <- if (is.list(result$ab1_evidence)) result$ab1_evidence$detail else NULL
  if (is.data.frame(detail) && "Basecaller_quality" %in% names(detail)) {
    valid <- raw_position >= 1L & raw_position <= nrow(detail)
    quality[valid] <- suppressWarnings(as.numeric(detail$Basecaller_quality[raw_position[valid]]))
  }

  peak_ratio <- rep(NA_real_, n)
  if (is.data.frame(result$metrics) && "peak_ratio" %in% names(result$metrics)) {
    valid <- raw_position >= 1L & raw_position <= nrow(result$metrics)
    peak_ratio[valid] <- suppressWarnings(as.numeric(result$metrics$peak_ratio[raw_position[valid]]))
  }

  direction <- stage3_scalar_text(direction, "Forward")
  if (identical(direction, "Reverse")) {
    sequence <- stage3_reverse_complement(sequence)
    quality <- rev(quality)
    peak_ratio <- rev(peak_ratio)
    raw_position <- rev(raw_position)
  }

  list(
    sequence = sequence,
    quality = quality,
    peak_ratio = peak_ratio,
    raw_position = as.integer(raw_position),
    oriented_position = seq_len(n)
  )
}

# Semi-global overlap alignment. Leading and trailing overhangs are free, while
# the shared part is scored. This is intentionally a small, deterministic
# alignment routine for two Sanger reads, not a multiple-sequence aligner.
stage3_overlap_align <- function(forward, reverse_oriented, match = 2, mismatch = -3, gap = -4) {
  a <- strsplit(stage3_sanitize_dna(forward), "", fixed = TRUE)[[1]]
  b <- strsplit(stage3_sanitize_dna(reverse_oriented), "", fixed = TRUE)[[1]]
  n <- length(a); m <- length(b)
  if (!n || !m) stop("Both Forward and Reverse sequences are required for overlap alignment.", call. = FALSE)

  score <- matrix(0, nrow = n + 1L, ncol = m + 1L)
  trace <- matrix(0L, nrow = n + 1L, ncol = m + 1L)
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      diag_score <- score[i, j] + if (a[i] == b[j]) match else mismatch
      up_score <- score[i, j + 1L] + gap
      left_score <- score[i + 1L, j] + gap
      values <- c(diag_score, up_score, left_score)
      direction <- which.max(values)
      score[i + 1L, j + 1L] <- values[direction]
      trace[i + 1L, j + 1L] <- direction
    }
  }

  last_row <- score[n + 1L, seq_len(m + 1L)]
  last_col <- score[seq_len(n + 1L), m + 1L]
  row_j <- which.max(last_row) - 1L
  col_i <- which.max(last_col) - 1L
  if (last_row[row_j + 1L] >= last_col[col_i + 1L]) {
    end_i <- n; end_j <- row_j
  } else {
    end_i <- col_i; end_j <- m
  }

  i <- end_i; j <- end_j
  aa <- bb <- character()
  apos <- bpos <- integer()
  while (i > 0L && j > 0L) {
    direction <- trace[i + 1L, j + 1L]
    if (direction == 1L) {
      aa <- c(aa, a[i]); bb <- c(bb, b[j]); apos <- c(apos, i); bpos <- c(bpos, j)
      i <- i - 1L; j <- j - 1L
    } else if (direction == 2L) {
      aa <- c(aa, a[i]); bb <- c(bb, "-"); apos <- c(apos, i); bpos <- c(bpos, NA_integer_)
      i <- i - 1L
    } else {
      aa <- c(aa, "-"); bb <- c(bb, b[j]); apos <- c(apos, NA_integer_); bpos <- c(bpos, j)
      j <- j - 1L
    }
  }
  while (i > 0L) {
    aa <- c(aa, a[i]); bb <- c(bb, "-"); apos <- c(apos, i); bpos <- c(bpos, NA_integer_); i <- i - 1L
  }
  while (j > 0L) {
    aa <- c(aa, "-"); bb <- c(bb, b[j]); apos <- c(apos, NA_integer_); bpos <- c(bpos, j); j <- j - 1L
  }
  aa <- rev(aa); bb <- rev(bb); apos <- rev(apos); bpos <- rev(bpos)

  if (end_i < n) {
    idx <- seq.int(end_i + 1L, n)
    aa <- c(aa, a[idx]); bb <- c(bb, rep("-", length(idx)))
    apos <- c(apos, idx); bpos <- c(bpos, rep(NA_integer_, length(idx)))
  }
  if (end_j < m) {
    idx <- seq.int(end_j + 1L, m)
    aa <- c(aa, rep("-", length(idx))); bb <- c(bb, b[idx])
    apos <- c(apos, rep(NA_integer_, length(idx))); bpos <- c(bpos, idx)
  }

  shared <- aa != "-" & bb != "-"
  overlap <- sum(shared)
  matches <- sum(shared & aa == bb)
  list(
    forward_aligned = paste0(aa, collapse = ""),
    reverse_aligned = paste0(bb, collapse = ""),
    forward_base = aa,
    reverse_base = bb,
    forward_position = as.integer(apos),
    reverse_position = as.integer(bpos),
    score = max(c(last_row, last_col)),
    overlap = as.integer(overlap),
    matches = as.integer(matches),
    identity_percent = if (overlap) round(100 * matches / overlap, 2) else NA_real_
  )
}

stage3_iupac_pair <- function(a, b) {
  key <- paste(sort(unique(c(a, b))), collapse = "")
  map <- c(AG = "R", CT = "Y", CG = "S", AT = "W", GT = "K", AC = "M")
  if (key %in% names(map)) unname(map[[key]]) else "N"
}

stage3_consensus_from_alignment <- function(alignment, forward_evidence, reverse_evidence,
                                            min_overlap = 40L, min_identity = 85,
                                            quality_delta = 10, strong_quality = 20) {
  overlap_ok <- alignment$overlap >= as.integer(min_overlap) &&
    is.finite(alignment$identity_percent) && alignment$identity_percent >= as.numeric(min_identity)
  ncol_aln <- length(alignment$forward_base)
  if (!ncol_aln) stop("Alignment contains no columns.", call. = FALSE)

  shared_cols <- which(alignment$forward_base != "-" & alignment$reverse_base != "-")
  first_shared <- if (length(shared_cols)) min(shared_cols) else Inf
  last_shared <- if (length(shared_cols)) max(shared_cols) else -Inf
  rows <- vector("list", ncol_aln)
  consensus <- character(ncol_aln)
  consensus_position <- 0L

  for (k in seq_len(ncol_aln)) {
    fb <- alignment$forward_base[k]; rb <- alignment$reverse_base[k]
    fp <- alignment$forward_position[k]; rp <- alignment$reverse_position[k]
    fq <- if (is.finite(fp) && fp <= length(forward_evidence$quality)) forward_evidence$quality[fp] else NA_real_
    rq <- if (is.finite(rp) && rp <= length(reverse_evidence$quality)) reverse_evidence$quality[rp] else NA_real_
    fraw <- if (is.finite(fp) && fp <= length(forward_evidence$raw_position)) forward_evidence$raw_position[fp] else NA_integer_
    rraw <- if (is.finite(rp) && rp <= length(reverse_evidence$raw_position)) reverse_evidence$raw_position[rp] else NA_integer_
    within_overlap <- k >= first_shared && k <= last_shared
    call <- "N"; decision <- "Unresolved"; review <- FALSE

    if (fb == "-" && rb != "-") {
      call <- rb
      if (within_overlap) { decision <- "Indel | Reverse base retained"; review <- TRUE } else decision <- "Reverse-only overhang"
    } else if (rb == "-" && fb != "-") {
      call <- fb
      if (within_overlap) { decision <- "Indel | Forward base retained"; review <- TRUE } else decision <- "Forward-only overhang"
    } else if (fb == rb) {
      call <- fb
      decision <- if (fb == "N") "Shared N" else "Reads agree"
      review <- fb == "N"
    } else if (fb == "N" && rb != "N") {
      call <- rb; decision <- "Reverse resolves Forward N"
    } else if (rb == "N" && fb != "N") {
      call <- fb; decision <- "Forward resolves Reverse N"
    } else {
      forward_wins <- is.finite(fq) && is.finite(rq) && fq >= strong_quality && fq - rq >= quality_delta
      reverse_wins <- is.finite(fq) && is.finite(rq) && rq >= strong_quality && rq - fq >= quality_delta
      if (forward_wins && !reverse_wins) {
        call <- fb; decision <- "Forward quality dominates"
      } else if (reverse_wins && !forward_wins) {
        call <- rb; decision <- "Reverse quality dominates"
      } else {
        call <- stage3_iupac_pair(fb, rb)
        decision <- "Conflict | IUPAC review call"
        review <- TRUE
      }
    }

    consensus_position <- consensus_position + 1L
    consensus[k] <- call
    rows[[k]] <- data.frame(
      Alignment_Column = k, Consensus_Position = consensus_position,
      Forward_Base = fb, Reverse_Base = rb, Consensus_Call = call,
      Forward_Quality = if (is.finite(fq)) round(fq, 2) else NA_real_,
      Reverse_Quality = if (is.finite(rq)) round(rq, 2) else NA_real_,
      Forward_Raw_Position = fraw, Reverse_Raw_Position = rraw,
      Decision = decision, Needs_Review = review,
      stringsAsFactors = FALSE
    )
  }

  evidence <- do.call(rbind, rows)
  if (!overlap_ok) {
    evidence$Needs_Review <- TRUE
    status <- "NO_RELIABLE_OVERLAP"
    sequence <- ""
  } else {
    status <- if (any(evidence$Needs_Review)) "REVIEW_REQUIRED" else "READY"
    sequence <- paste0(consensus, collapse = "")
  }

  list(
    status = status,
    sequence = sequence,
    evidence = evidence,
    overlap = alignment$overlap,
    identity_percent = alignment$identity_percent,
    mismatches = sum(alignment$forward_base != "-" & alignment$reverse_base != "-" & alignment$forward_base != alignment$reverse_base),
    indels = sum(seq_len(ncol_aln) >= first_shared & seq_len(ncol_aln) <= last_shared &
                   xor(alignment$forward_base == "-", alignment$reverse_base == "-")),
    review_positions = sum(evidence$Needs_Review),
    alignment = alignment
  )
}

stage3_single_read_record <- function(assignment, result, consensus_id, built_at, independent = FALSE) {
  direction <- stage3_scalar_text(assignment$Direction)
  ev <- stage3_result_evidence(result, direction)
  seq_chars <- strsplit(ev$sequence, "", fixed = TRUE)[[1]]
  evidence <- if (length(seq_chars)) data.frame(
    Alignment_Column = seq_along(seq_chars), Consensus_Position = seq_along(seq_chars),
    Forward_Base = if (direction == "Forward") seq_chars else "-",
    Reverse_Base = if (direction == "Reverse") seq_chars else "-",
    Consensus_Call = seq_chars,
    Forward_Quality = if (direction == "Forward") ev$quality else NA_real_,
    Reverse_Quality = if (direction == "Reverse") ev$quality else NA_real_,
    Forward_Raw_Position = if (direction == "Forward") ev$raw_position else NA_integer_,
    Reverse_Raw_Position = if (direction == "Reverse") ev$raw_position else NA_integer_,
    Decision = paste0(direction, " single-read representative"), Needs_Review = FALSE,
    stringsAsFactors = FALSE
  ) else data.frame()

  list(
    schema = "pitax-consensus-record-v1", consensus_id = consensus_id,
    final_name = if (isTRUE(independent)) stage3_scalar_text(assignment$Final_Name, stage3_compose_sequence_name(assignment$Isolate, assignment$Locus)) else stage3_compose_sequence_name(assignment$Isolate, assignment$Locus),
    isolate = stage3_scalar_text(assignment$Isolate), locus = stage3_scalar_text(assignment$Locus),
    status = if (nzchar(ev$sequence)) if (isTRUE(independent)) "INDEPENDENT_READ" else "SINGLE_READ" else "SOURCE_MISSING",
    sequence = ev$sequence, source_read_ids = stage3_scalar_text(assignment$Source_ID),
    forward_read = if (direction == "Forward") stage3_scalar_text(assignment$Source_ID) else "",
    reverse_read = if (direction == "Reverse") stage3_scalar_text(assignment$Source_ID) else "",
    source_sequences = setNames(list(stage3_scalar_text(result$seq)), stage3_scalar_text(assignment$Source_ID)),
    source_revisions = setNames(as.list(if (is.list(result$curation)) as.integer(result$curation$revision) else 0L), stage3_scalar_text(assignment$Source_ID)),
    alignment = NULL, evidence = evidence,
    metrics = list(overlap = 0L, identity_percent = NA_real_, mismatches = 0L, indels = 0L, review_positions = 0L),
    algorithm = if (isTRUE(independent)) "independent-read-orientation-v1" else "single-read-orientation-normalization-v1", built_at = built_at
  )
}

stage3_pair_record <- function(assignments, results, consensus_id, settings, built_at) {
  frow <- assignments[assignments$Direction == "Forward", , drop = FALSE][1, , drop = FALSE]
  rrow <- assignments[assignments$Direction == "Reverse", , drop = FALSE][1, , drop = FALSE]
  fresult <- results[[stage3_scalar_text(frow$Source_ID)]]
  rresult <- results[[stage3_scalar_text(rrow$Source_ID)]]
  base_record <- list(
    schema = "pitax-consensus-record-v1", consensus_id = consensus_id,
    final_name = stage3_compose_sequence_name(frow$Isolate, frow$Locus),
    isolate = stage3_scalar_text(frow$Isolate), locus = stage3_scalar_text(frow$Locus),
    source_read_ids = c(stage3_scalar_text(frow$Source_ID), stage3_scalar_text(rrow$Source_ID)),
    forward_read = stage3_scalar_text(frow$Source_ID), reverse_read = stage3_scalar_text(rrow$Source_ID),
    algorithm = "pitax-overlap-consensus-v1", built_at = built_at
  )
  if (!is.list(fresult) || !is.list(rresult) || !nzchar(stage3_sanitize_dna(fresult$seq)) || !nzchar(stage3_sanitize_dna(rresult$seq))) {
    return(c(base_record, list(
      status = "SOURCE_MISSING", sequence = "", source_sequences = list(), source_revisions = list(),
      alignment = NULL, evidence = data.frame(),
      metrics = list(overlap = 0L, identity_percent = NA_real_, mismatches = 0L, indels = 0L, review_positions = 0L)
    )))
  }

  fe <- stage3_result_evidence(fresult, "Forward")
  re <- stage3_result_evidence(rresult, "Reverse")
  alignment <- stage3_overlap_align(
    fe$sequence, re$sequence,
    match = settings$alignment_match,
    mismatch = settings$alignment_mismatch,
    gap = settings$alignment_gap
  )
  call <- stage3_consensus_from_alignment(
    alignment, fe, re,
    min_overlap = settings$min_overlap,
    min_identity = settings$min_identity,
    quality_delta = settings$quality_delta,
    strong_quality = settings$strong_quality
  )
  c(base_record, list(
    status = call$status, sequence = call$sequence,
    source_sequences = setNames(list(stage3_scalar_text(fresult$seq), stage3_scalar_text(rresult$seq)), c(base_record$forward_read, base_record$reverse_read)),
    source_revisions = setNames(list(
      if (is.list(fresult$curation)) as.integer(fresult$curation$revision) else 0L,
      if (is.list(rresult$curation)) as.integer(rresult$curation$revision) else 0L
    ), c(base_record$forward_read, base_record$reverse_read)),
    alignment = call$alignment, evidence = call$evidence,
    metrics = list(overlap = call$overlap, identity_percent = call$identity_percent,
                   mismatches = call$mismatches, indels = call$indels, review_positions = call$review_positions)
  ))
}

stage3_summary_row <- function(record) {
  record <- stage3_ensure_record_curation(record)
  data.frame(
    Consensus_ID = stage3_scalar_text(record$consensus_id), Final_Name = stage3_scalar_text(record$final_name),
    Isolate = stage3_scalar_text(record$isolate), Locus = stage3_scalar_text(record$locus),
    Status = stage3_scalar_text(record$status), Source_Reads = length(record$source_read_ids),
    Forward_Read = stage3_scalar_text(record$forward_read), Reverse_Read = stage3_scalar_text(record$reverse_read),
    Length = nchar(stage3_scalar_text(record$sequence)), Overlap = as.integer(record$metrics$overlap),
    Identity_percent = suppressWarnings(as.numeric(record$metrics$identity_percent)),
    Mismatches = as.integer(record$metrics$mismatches), Indels = as.integer(record$metrics$indels),
    Review_positions = as.integer(record$metrics$review_positions), Revision = as.integer(record$curation$revision),
    Algorithm = stage3_scalar_text(record$algorithm),
    Built_at = stage3_scalar_text(record$built_at), stringsAsFactors = FALSE
  )
}

stage3_build_consensus_set <- function(assignments, results, min_overlap = 40L, min_identity = 85,
                                       quality_delta = 10, strong_quality = 20,
                                       project_mode = c("paired_consensus", "simple")) {
  project_mode <- match.arg(project_mode)
  locus_error <- stage3_run_locus_error(assignments)
  if (!is.null(locus_error)) stop(locus_error, call. = FALSE)
  if (!is.list(results) || !length(results)) stop("No processed read results are available.", call. = FALSE)

  settings <- list(
    project_mode = project_mode,
    min_overlap = as.integer(min_overlap), min_identity = as.numeric(min_identity),
    quality_delta = as.numeric(quality_delta), strong_quality = as.numeric(strong_quality),
    alignment_match = 2, alignment_mismatch = -3, alignment_gap = -4
  )
  if (!is.finite(settings$min_overlap) || settings$min_overlap < 1L) stop("Minimum overlap must be at least 1 bp.", call. = FALSE)
  if (!is.finite(settings$min_identity) || settings$min_identity < 0 || settings$min_identity > 100) stop("Minimum identity must be between 0 and 100%.", call. = FALSE)
  if (!is.finite(settings$quality_delta) || settings$quality_delta < 0) stop("Quality delta must be non-negative.", call. = FALSE)
  if (!is.finite(settings$strong_quality) || settings$strong_quality < 0) stop("Strong base quality must be non-negative.", call. = FALSE)

  keys <- if (identical(project_mode, "simple")) {
    paste0("read\r", seq_len(nrow(assignments)))
  } else {
    paste(trimws(as.character(assignments$Isolate)), trimws(as.character(assignments$Locus)), sep = "\r")
  }
  groups <- split(seq_len(nrow(assignments)), keys)
  built_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  records <- list()
  for (g in seq_along(groups)) {
    rows <- assignments[groups[[g]], , drop = FALSE]
    consensus_id <- sprintf("consensus_%03d", g)
    if (nrow(rows) == 1L) {
      source <- stage3_scalar_text(rows$Source_ID)
      record <- stage3_single_read_record(rows, results[[source]], consensus_id, built_at, independent = identical(project_mode, "simple"))
    } else if (nrow(rows) == 2L && identical(sort(as.character(rows$Direction)), c("Forward", "Reverse"))) {
      record <- stage3_pair_record(rows, results, consensus_id, settings, built_at)
    } else {
      stop(paste0("Isolate/locus ", rows$Isolate[1], "/", rows$Locus[1], " has multiple reads but not one explicit Forward/Reverse pair."), call. = FALSE)
    }
    records[[consensus_id]] <- stage3_ensure_record_curation(record)
  }
  summary <- if (length(records)) do.call(rbind, lapply(records, stage3_summary_row)) else stage3_empty_summary()
  rownames(summary) <- NULL
  list(
    schema = "pitax-consensus-set-v1", algorithm = "pitax-overlap-consensus-v1",
    built_at = built_at, settings = settings, records = records, summary = summary
  )
}

stage3_refresh_consensus_summary <- function(consensus_set) {
  if (!is.list(consensus_set)) return(stage3_empty_consensus_set())
  records <- if (is.list(consensus_set$records)) consensus_set$records else list()
  consensus_set$summary <- if (length(records)) do.call(rbind, lapply(records, stage3_summary_row)) else stage3_empty_summary()
  rownames(consensus_set$summary) <- NULL
  consensus_set
}

stage3_ensure_consensus_set <- function(consensus_set) {
  if (!is.list(consensus_set) || !identical(stage3_scalar_text(consensus_set$schema), "pitax-consensus-set-v1")) {
    return(stage3_empty_consensus_set())
  }
  if (!is.list(consensus_set$records)) consensus_set$records <- list()
  consensus_set$records <- lapply(consensus_set$records, stage3_ensure_record_curation)
  stage3_refresh_consensus_summary(consensus_set)
}

stage3_consensus_set_is_current <- function(consensus_set, results) {
  if (!is.list(consensus_set) || !length(consensus_set$records) || !is.list(results)) return(FALSE)
  for (record in consensus_set$records) {
    ids <- as.character(record$source_read_ids)
    if (!length(ids) || any(!ids %in% names(results))) return(FALSE)
    for (id in ids) {
      current <- stage3_scalar_text(results[[id]]$seq)
      stored <- stage3_scalar_text(record$source_sequences[[id]])
      if (!identical(current, stored)) return(FALSE)
    }
  }
  TRUE
}

stage3_consensus_gate_error <- function(consensus_set, results) {
  if (!is.list(consensus_set) || !length(consensus_set$records)) return("Build the analysis sequences before continuing.")
  if (!stage3_consensus_set_is_current(consensus_set, results)) return("One or more curated source reads changed. Rebuild the analysis sequences.")
  statuses <- vapply(consensus_set$records, function(x) stage3_scalar_text(x$status), character(1))
  if (any(statuses == "NO_RELIABLE_OVERLAP")) return("At least one Forward/Reverse pair has no reliable overlap.")
  if (any(statuses == "SOURCE_MISSING")) return("At least one isolate/locus is missing a usable processed source read.")
  if (any(statuses == "REVIEW_REQUIRED")) return("At least one consensus contains unresolved conflicts that require review.")
  NULL
}

stage3_analysis_records <- function(consensus_set) {
  if (!is.list(consensus_set) || !length(consensus_set$records)) return(list())
  out <- list()
  for (id in names(consensus_set$records)) {
    record <- consensus_set$records[[id]]
    if (!nzchar(stage3_scalar_text(record$sequence))) next
    record <- stage3_ensure_record_curation(record)
    analysis_reason <- switch(stage3_scalar_text(record$status),
      SINGLE_READ = "SINGLE_READ_REPRESENTATIVE",
      INDEPENDENT_READ = "INDEPENDENT_READ_ORIENTED",
      "FORWARD_REVERSE_CONSENSUS"
    )
    out[[id]] <- list(
      sample_id = id, seq = stage3_scalar_text(record$sequence), final_name = stage3_scalar_text(record$final_name),
      consensus = record,
      summary = data.frame(
        sample_id = id, target = stage3_scalar_text(record$locus), raw_length = NA_integer_,
        trimmed_length = nchar(stage3_scalar_text(record$sequence)), trim_start = NA_integer_, trim_end = NA_integer_,
        reason = analysis_reason, status = stage3_scalar_text(record$status),
        consensus_revision = as.integer(record$curation$revision), stringsAsFactors = FALSE
      )
    )
  }
  out
}

stage3_analysis_summary <- function(consensus_set) {
  records <- stage3_analysis_records(consensus_set)
  if (!length(records)) return(data.frame())
  out <- do.call(rbind, lapply(records, function(x) x$summary))
  out$final_name <- vapply(records, function(x) stage3_scalar_text(x$final_name), character(1))
  rownames(out) <- NULL
  out
}

stage3_migrate_v3_state <- function(state) {
  if (!is.list(state)) state <- list()
  state$consensus_set <- stage3_empty_consensus_set()
  old_log <- stage3_scalar_text(state$migration_log)
  note <- "Migrated project schema 3 to schema 4. Existing read/QC/BLAST/taxonomy evidence was preserved; Stage 3 isolate-level sequences must be built explicitly."
  state$migration_log <- if (nzchar(old_log)) paste(old_log, note, sep = " ") else note
  state
}
