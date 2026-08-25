# PITAX 3.0.4 - reflection follow-up correctness smokes.

source(file.path("R", "domain", "assay", "assay_profiles.R"))
source(file.path("R", "domain", "assignment", "stage2_architecture.R"))
source(file.path("tests", "helpers", "test_paths.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

# --- Assay_ID belongs in the assignment identity signature ---
assignment_state_signature <- function(assignments) {
  assignments <- stage2_coerce_assignments(assignments)
  if (!nrow(assignments)) return("")
  cols <- c("Source_ID", "Isolate", "Assay_ID", "Locus", "Direction", "Final_Name")
  paste(apply(assignments[, cols, drop = FALSE], 1, paste, collapse = "\r"), collapse = "\n")
}

profiles <- rbind(
  assay_profile_from_legacy_settings(list(target = "TEF1", forward_primer = "EF1", reverse_primer = "EF2"), assay_id = "assay-tef1-a"),
  assay_profile_from_legacy_settings(list(target = "TEF1", forward_primer = "EF3", reverse_primer = "EF4"), assay_id = "assay-tef1-b")
)
a1 <- stage2_coerce_assignments(data.frame(
  Source_ID = "S1", Isolate = "ISO1", Locus = "TEF1", Direction = "Forward",
  Primer = "EF1", Notes = "", Assay_ID = "assay-tef1-a",
  File = "S1.ab1", Read_ID = "r1", Final_Name = "ISO1_TEF1_F",
  stringsAsFactors = FALSE
))
a1 <- stage2_sync_generated_names(a1, assay_profiles = profiles)
a2 <- a1
a2$Assay_ID <- "assay-tef1-b"
a2 <- stage2_sync_generated_names(a2, assay_profiles = profiles)
assert_true(identical(a1$Final_Name, a2$Final_Name), "Same-locus assay swap should keep Final_Name.")
assert_true(!identical(assignment_state_signature(a1), assignment_state_signature(a2)),
            "Assay_ID changes must alter the assignment identity signature.")

# --- Taxonomy RID selection mirrors READY-only policy ---
latest_blast_rid_for_sample <- function(jobs, original_name) {
  jobs <- jobs[jobs$original_name == original_name, , drop = FALSE]
  if (!nrow(jobs)) return("")
  ready <- jobs[as.character(jobs$status) == "READY", , drop = FALSE]
  if (!nrow(ready)) return("")
  as.character(ready$rid[nrow(ready)])
}
jobs <- data.frame(
  original_name = c("ISO1_ITS", "ISO1_ITS", "ISO1_ITS"),
  rid = c("RID_OLD", "RID_READY", "RID_STALE"),
  status = c("READY", "READY", "STALE"),
  stringsAsFactors = FALSE
)
assert_true(identical(latest_blast_rid_for_sample(jobs, "ISO1_ITS"), "RID_READY"),
            "Taxonomy must prefer the latest READY RID, not a trailing STALE row.")
assert_true(identical(latest_blast_rid_for_sample(jobs[3, , drop = FALSE], "ISO1_ITS"), ""),
            "Taxonomy must return empty when no READY RID exists.")

# --- Curation invalidation breadth: clear-all consensus invalidates every prior ID ---
prior_ids <- c("ISO1_ITS", "ISO2_ITS")
sample_name <- "readA_F"
invalidated <- unique(c(prior_ids, sample_name))
assert_true(all(c("ISO1_ITS", "ISO2_ITS", "readA_F") %in% invalidated),
            "When the consensus set is cleared, every prior consensus ID must be invalidated.")

# --- Contract markers in server sources ---
upload_text <- pitax_read_text("R", "server", "stages", "20_upload.R")
project_text <- pitax_read_text("R", "server", "stages", "10_project.R")
blast_text <- pitax_read_text("R", "server", "stages", "90_blast.R")
taxonomy_text <- pitax_read_text("R", "server", "stages", "100_taxonomy.R")
assert_true(grepl('cols <- c("Source_ID", "Isolate", "Assay_ID", "Locus", "Direction", "Final_Name")', upload_text, fixed = TRUE),
            "Upload signature must include Assay_ID.")
assert_true(grepl("prior_consensus_ids", project_text, fixed = TRUE),
            "Project curation must invalidate all prior consensus IDs after a sequence change.")
assert_true(grepl("job_revision", blast_text, fixed = TRUE) && grepl("current_revision", blast_text, fixed = TRUE),
            "BLAST retrieve must compare stored and current consensus revisions.")
assert_true(grepl('jobs\\[as\\.character\\(jobs\\$status\\) == "READY"', taxonomy_text),
            "Taxonomy RID helper must filter READY jobs.")

cat("Reflection follow-up correctness smokes passed.\n")
