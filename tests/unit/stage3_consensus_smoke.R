# PITAX v3.0.0-alpha.8 - Stage 3 consensus and curation tests.

source(file.path("R", "domain", "consensus", "stage3_consensus.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

rc <- stage3_reverse_complement
assert_true(rc("ARYKMBVDHN") == "NDHBVKMRYT", "IUPAC ambiguity codes were not reverse-complemented correctly.")

make_result <- function(sequence, quality = rep(30, nchar(sequence)), revision = 0L) {
  n <- nchar(sequence)
  data.frame_quality <- data.frame(Basecaller_quality = as.numeric(quality), stringsAsFactors = FALSE)
  list(
    seq = sequence,
    curation = list(trim_start = 1L, trim_end = n, revision = revision),
    summary = data.frame(trim_start = 1L, trim_end = n, stringsAsFactors = FALSE),
    ab1_evidence = list(detail = data.frame_quality),
    metrics = data.frame(peak_ratio = rep(3, n), stringsAsFactors = FALSE)
  )
}

assignment_pair <- data.frame(
  Source_ID = c("F001", "R001"), Isolate = c("ISO1", "ISO1"), Locus = c("ITS", "ITS"),
  Direction = c("Forward", "Reverse"), stringsAsFactors = FALSE
)

# 1. Reverse-complement orientation and suffix/prefix overlap reconstruct one
# isolate-level sequence without changing either stored source read.
truth <- "ACGTTGCAACCTGGAATC"
forward <- substr(truth, 1, 14)
reverse_oriented <- substr(truth, 7, nchar(truth))
reverse_raw <- rc(reverse_oriented)
results <- list(F001 = make_result(forward), R001 = make_result(reverse_raw))
built <- stage3_build_consensus_set(assignment_pair, results, min_overlap = 6, min_identity = 90)
record <- built$records$consensus_001
assert_true(record$status == "READY", "A clean Forward/Reverse overlap was not ready.")
assert_true(record$sequence == truth, "The clean pair did not reconstruct the expected isolate-level sequence.")
assert_true(record$forward_read == "F001" && record$reverse_read == "R001", "Source-read provenance was lost.")
assert_true(is.null(stage3_consensus_gate_error(built, results)), "A clean pair did not pass the Stage 3 gate.")

# 2. A single Reverse read remains valid but is explicitly labeled as a
# single-read representative, not as a two-read consensus.
single_assignment <- data.frame(Source_ID = "R002", Isolate = "ISO2", Locus = "ITS", Direction = "Reverse", stringsAsFactors = FALSE)
single_oriented <- "GATTACAGG"
single_results <- list(R002 = make_result(rc(single_oriented)))
single <- stage3_build_consensus_set(single_assignment, single_results, min_overlap = 4)
single_record <- single$records$consensus_001
assert_true(single_record$status == "SINGLE_READ", "A single read was not labeled as a single-read representative.")
assert_true(single_record$sequence == single_oriented, "A single Reverse read was not normalized to the Forward biological orientation.")
assert_true(is.null(stage3_consensus_gate_error(single, single_results)), "A valid single-read representative was incorrectly blocked.")

# 2b. Simple-project mode never groups reads by Isolate/Locus. It keeps the
# explicit read-level name and still orients Reverse reads for downstream use.
simple_assignment <- data.frame(
  Source_ID = c("F003", "R003"), Final_Name = c("ISO3_ITS_F", "ISO3_ITS_R"),
  Isolate = c("ISO3", "ISO3"), Locus = c("ITS", "ITS"),
  Direction = c("Forward", "Reverse"), stringsAsFactors = FALSE
)
simple_results <- list(F003 = make_result("ACGTACGT"), R003 = make_result(rc("GATTACAG")))
simple <- stage3_build_consensus_set(simple_assignment, simple_results, min_overlap = 4, project_mode = "simple")
assert_true(length(simple$records) == 2L, "Simple-project reads were incorrectly grouped into a pair.")
assert_true(all(simple$summary$Status == "INDEPENDENT_READ"), "Simple-project reads were not labeled as independent.")
assert_true(identical(as.character(simple$summary$Final_Name), c("ISO3_ITS_F", "ISO3_ITS_R")), "Simple-project read names were collapsed.")
assert_true(simple$records$consensus_002$sequence == "GATTACAG", "Simple-project Reverse read was not oriented.")

# 3. One high-quality base may resolve a low-quality contradiction when the
# explicit quality advantage is large enough.
f_conflict <- "ACGTACGT"
r_oriented <- "ACGTTCGT"
fq <- rep(30, 8); fq[5] <- 40
rq_oriented <- rep(30, 8); rq_oriented[5] <- 10
quality_results <- list(
  F001 = make_result(f_conflict, fq),
  R001 = make_result(rc(r_oriented), rev(rq_oriented))
)
quality_built <- stage3_build_consensus_set(assignment_pair, quality_results, min_overlap = 8, min_identity = 80, quality_delta = 10, strong_quality = 20)
quality_record <- quality_built$records$consensus_001
assert_true(quality_record$status == "READY", "A strongly asymmetric quality conflict was not resolved.")
assert_true(quality_record$sequence == f_conflict, "The higher-quality Forward base did not win the conflict.")
assert_true(any(quality_record$evidence$Decision == "Forward quality dominates"), "The quality-resolution provenance was not recorded.")

# 4. Strong contradictory evidence is not silently voted away. It produces an
# IUPAC ambiguity call and keeps the downstream gate closed for review.
strong_results <- list(
  F001 = make_result(f_conflict, rep(30, 8)),
  R001 = make_result(rc(r_oriented), rep(30, 8))
)
strong_built <- stage3_build_consensus_set(assignment_pair, strong_results, min_overlap = 8, min_identity = 80)
strong_record <- strong_built$records$consensus_001
assert_true(strong_record$status == "REVIEW_REQUIRED", "Strong contradictory evidence was not sent to review.")
assert_true(substr(strong_record$sequence, 5, 5) == "W", "The unresolved A/T conflict was not retained as an IUPAC W call.")
assert_true(grepl("review", stage3_consensus_gate_error(strong_built, strong_results), ignore.case = TRUE), "The review-required record did not close the downstream gate.")

# 4b. A reviewer can accept either source call or intentionally retain the
# automatic IUPAC call. Every action is revisioned and Undo/Redo are themselves
# new, auditable revisions rather than a rollback of history.
review_column <- strong_record$evidence$Alignment_Column[strong_record$evidence$Needs_Review][1]
forward_review <- stage3_apply_consensus_call(strong_record, review_column, "A", method = "Forward", note = "Forward trace is cleaner")
assert_true(forward_review$status == "READY", "A reviewed conflict did not open the consensus gate.")
assert_true(substr(forward_review$sequence, 5, 5) == "A", "The accepted Forward call was not materialized.")
assert_true(forward_review$curation$revision == 1L && nrow(forward_review$curation$audit_log) == 1L, "The consensus review was not revisioned and audited.")
reviewed_set <- strong_built
reviewed_set$records$consensus_001 <- forward_review
reviewed_set <- stage3_refresh_consensus_summary(reviewed_set)
assert_true(is.null(stage3_consensus_gate_error(reviewed_set, strong_results)), "A fully reviewed consensus set did not pass the downstream gate.")
assert_true(reviewed_set$summary$Revision[1] == 1L && reviewed_set$summary$Review_positions[1] == 0L, "The summary did not refresh after consensus curation.")

undone_review <- stage3_consensus_undo(forward_review)
assert_true(undone_review$status == "REVIEW_REQUIRED" && undone_review$curation$revision == 2L, "Undo did not restore the unresolved state as a new revision.")
assert_true(substr(undone_review$sequence, 5, 5) == "W", "Undo did not restore the automatic IUPAC call.")
redone_review <- stage3_consensus_redo(undone_review)
assert_true(redone_review$status == "READY" && redone_review$curation$revision == 3L, "Redo did not restore the reviewed state as a new revision.")
assert_true(identical(redone_review$source_revisions, strong_record$source_revisions), "Consensus curation changed its source-revision binding.")

iupac_review <- stage3_apply_consensus_call(strong_record, review_column, "W", method = "IUPAC", note = "Both traces remain plausible")
assert_true(iupac_review$status == "READY", "Intentional acceptance of an IUPAC call remained blocked.")
assert_true(isTRUE(iupac_review$evidence$Manual_Review[match(review_column, iupac_review$evidence$Alignment_Column)]), "Intentional IUPAC acceptance was not marked as manual review.")

# 5. A weak/unrelated overlap is blocked and emits no downstream sequence.
no_overlap_results <- list(F001 = make_result("AAAAAAAAAAAA"), R001 = make_result("GGGGGGGGGGGG"))
blocked <- stage3_build_consensus_set(assignment_pair, no_overlap_results, min_overlap = 8, min_identity = 75)
assert_true(blocked$records$consensus_001$status == "NO_RELIABLE_OVERLAP", "An unrelated pair was not blocked.")
assert_true(blocked$records$consensus_001$sequence == "", "A blocked pair emitted a downstream sequence.")

# 6. One upload/run is single-locus even though the project architecture can
# later combine separately processed loci in Stage 4.
mixed <- rbind(assignment_pair, data.frame(Source_ID = "F002", Isolate = "ISO2", Locus = "TEF1", Direction = "Forward", stringsAsFactors = FALSE))
assert_true(grepl("exactly one", stage3_run_locus_error(mixed), fixed = TRUE), "A mixed-locus run was not blocked.")

# 7. Consensus becomes stale when a curated source sequence changes.
changed <- results
changed$F001$seq <- paste0(changed$F001$seq, "A")
assert_true(!stage3_consensus_set_is_current(built, changed), "A source-sequence change did not stale the consensus.")

# 8. Schema migration is additive and preserves established evidence.
legacy <- list(results = list(marker = "keep"), blast_hits = data.frame(accession = "NR_1"), migration_log = "Stage 2 closed.")
migrated <- stage3_migrate_v3_state(legacy)
assert_true(identical(migrated$results$marker, "keep"), "Schema migration changed read evidence.")
assert_true(identical(migrated$blast_hits$accession, "NR_1"), "Schema migration changed BLAST evidence.")
assert_true(identical(migrated$consensus_set$schema, "pitax-consensus-set-v1"), "Schema migration did not initialize the Stage 3 store.")

cat("v3.0.0-alpha.8 Stage 3 consensus tests passed.\n")
