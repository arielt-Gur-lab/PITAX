# One-read edits must not re-align all unrelated records or lose their
# manually reviewed consensus and BLAST stable IDs.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "domain", "consensus", "stage3_consensus.R"))

make_assignments <- function(ids, mode = "simple") {
  directions <- if (identical(mode, "simple")) rep("Forward", length(ids)) else
    rep(c("Forward", "Reverse"), length.out = length(ids))
  isolates <- if (identical(mode, "simple")) ids else
    rep(sprintf("ISO%03d", seq_len(length(ids) / 2L)), each = 2L)
  data.frame(
    Source_ID = ids, Isolate = isolates, Locus = "TEF1",
    Name_Label = "TEF1", Direction = directions,
    Final_Name = paste0(ids, "_TEF1_", ifelse(directions == "Forward", "F", "R")),
    stringsAsFactors = FALSE
  )
}
settings <- list(min_overlap = 10L, min_identity = 85,
                 quality_delta = 10, strong_quality = 20)
build <- function(assignments, results, mode) {
  stage3_build_consensus_set(
    assignments, results,
    min_overlap = settings$min_overlap,
    min_identity = settings$min_identity,
    quality_delta = settings$quality_delta,
    strong_quality = settings$strong_quality,
    project_mode = mode
  )
}
rebuild <- function(previous, assignments, results, id, mode, threshold = settings$min_overlap) {
  stage3_rebuild_changed_source(
    previous, assignments, results, source_id = id,
    min_overlap = threshold,
    min_identity = settings$min_identity,
    quality_delta = settings$quality_delta,
    strong_quality = settings$strong_quality,
    project_mode = mode
  )
}

# Realistic 50-read simple project: change only sample 25, keep all 49 others
# byte-for-byte including their prior curation revisions.
ids <- sprintf("FB%03d", seq_len(50))
assignments <- make_assignments(ids)
results <- stats::setNames(lapply(ids, function(id) list(seq = paste0(strrep("A", 20L), "N"),
  summary = data.frame(trim_start = 1L, trim_end = 21L))),
  ids)
original <- build(assignments, results, "simple")
stopifnot(length(original$records) == 50L)
original$records[[2]]$curation$revision <- 12L
target <- names(original$records)[which(vapply(original$records, function(x) "FB025" %in% x$source_read_ids, logical(1)))]
stopifnot(length(target) == 1L)
results[["FB025"]]$seq <- paste0(strrep("A", 20L), "V")
updated <- rebuild(original, assignments, results, "FB025", "simple")
stopifnot(is.list(updated))
stopifnot(identical(names(updated$records), names(original$records)))
stopifnot(identical(updated$records[[2]], original$records[[2]]))
stopifnot(updated$records[[target]]$source_sequences[["FB025"]] == results[["FB025"]]$seq)
stopifnot(updated$records[[target]]$sequence == results[["FB025"]]$seq)
stopifnot(isTRUE(stage3_consensus_set_is_current(updated, results)))
stopifnot(identical(updated$records[[2]]$curation$revision, 12L))
stopifnot(nrow(updated$summary) == 50L)

# Any renamed/new/removed read or altered consensus alignment settings must
# fall back to a complete rebuild: never produce a silent partial mismerge.
stopifnot(is.null(rebuild(original, assignments, results, "missing", "simple")))
renamed <- assignments
renamed$Source_ID[2] <- "NEW_SOURCE"
stopifnot(is.null(rebuild(original, renamed, results, "FB025", "simple")))
stopifnot(is.null(rebuild(original, assignments, results, "FB025", "simple", threshold = 11L)))
missing_record <- original
missing_record$records <- missing_record$records[-1]
stopifnot(is.null(rebuild(missing_record, assignments, results, "FB025", "simple")))

# Paired project: changing one Forward read should re-align only that pair,
# not discard previous reviews of an unrelated pair.
paired_ids <- c("ISO001_F", "ISO001_R", "ISO002_F", "ISO002_R")
paired_assign <- make_assignments(paired_ids, mode = "paired_consensus")
paired_results <- stats::setNames(lapply(seq_along(paired_ids), function(i) {
  list(seq = if (i %% 2L == 1L) strrep("A", 30L) else strrep("T", 30L),
       summary = data.frame(trim_start = 1L, trim_end = 30L))
}), paired_ids)
pair_old <- build(paired_assign, paired_results, "paired_consensus")
stopifnot(length(pair_old$records) == 2L)
other <- names(pair_old$records)[vapply(pair_old$records, function(x) "ISO002_F" %in% x$source_read_ids, logical(1))]
target_pair <- names(pair_old$records)[vapply(pair_old$records, function(x) "ISO001_F" %in% x$source_read_ids, logical(1))]
stopifnot(length(other) == 1L, length(target_pair) == 1L)
pair_old$records[[other]]$curation$revision <- 8L
paired_results[["ISO001_F"]]$seq <- paste0("V", strrep("A", 29L))
pair_new <- rebuild(pair_old, paired_assign, paired_results, "ISO001_F", "paired_consensus")
stopifnot(is.list(pair_new))
stopifnot(identical(pair_old$records[[other]], pair_new$records[[other]]))
stopifnot(pair_new$records[[target_pair]]$source_sequences[["ISO001_F"]] == paired_results[["ISO001_F"]]$seq)
stopifnot(isTRUE(stage3_consensus_set_is_current(pair_new, paired_results)))

cat("incremental_consensus_smoke: OK\n")
