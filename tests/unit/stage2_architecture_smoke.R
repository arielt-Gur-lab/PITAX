# PITAX v3.0.0-alpha.6 - Stage 2 explicit identity and migration tests.

source(file.path("R", "domain", "assignment", "stage2_architecture.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

# Upload filenames/barcodes are immutable source identifiers, never biological identity.
raw <- stage2_initialize_read_assignment(
  "ISO7_ITS_F.ab1", default_locus = "LSU", default_direction = "Reverse",
  forward_primer = "ITS1F", reverse_primer = "ITS4"
)
assert_true(raw$Source_ID == "ISO7_ITS_F", "Upload source barcode was not preserved.")
assert_true(raw$Isolate == "", "PITAX incorrectly inferred isolate identity from the filename.")
assert_true(raw$Locus == "LSU" && raw$Direction == "Reverse", "Visible assay defaults were not used only as initialization values.")
assert_true(raw$Final_Name == "", "A final name was generated before explicit isolate assignment.")

# Explicit fields are the source of truth; PITAX generates the output name from them.
raw$Isolate <- "FB_120"
raw$Locus <- "TEF1 / EF1-alpha"
raw$Direction <- "Forward"
assigned <- stage2_sync_generated_names(raw, forward_primer = "EF1-728F")
assert_true(assigned$Source_ID == "ISO7_ITS_F", "Assignment changed the immutable upload barcode.")
assert_true(assigned$Isolate == "FB_120", "Explicit isolate field was not preserved.")
assert_true(assigned$Final_Name == "FB_120_TEF1-EF1-alpha_F", "Final name was not generated from explicit fields.")
assert_true(assigned$Primer == "EF1-728F", "Primer metadata was not synchronized by explicit direction.")
assert_true(is.null(stage2_identity_error(assigned)), "Complete explicit identity was rejected.")

missing <- assigned
missing$Isolate <- ""
missing <- stage2_sync_generated_names(missing)
assert_true(grepl("isolate", stage2_identity_error(missing), ignore.case = TRUE), "Missing explicit isolate was not blocked.")

# Two distinct source reads form one pair only when both explicit assignments share isolate/locus.
paired <- stage2_make_read_assignments(
  c("barcode_001.ab1", "barcode_002.ab1"),
  default_locus = "ITS", default_direction = "Unknown"
)
paired$Isolate <- c("ISO7", "ISO7")
paired$Locus <- c("ITS", "ITS")
paired$Direction <- c("Forward", "Reverse")
paired <- stage2_sync_generated_names(paired, forward_primer = "ITS1F", reverse_primer = "ITS4")
architecture <- stage2_build_architecture(paired, project_id = "pair_test")
summary <- stage2_architecture_summary(architecture)
assert_true(summary$Reads == 2L, "Both source reads must remain distinct.")
assert_true(summary$Paired_loci == 1L, "Explicit Forward/Reverse pair was not counted.")
assert_true(nrow(architecture$isolates) == 1L && nrow(architecture$loci) == 1L, "Pair did not share one isolate and locus.")

# A single read remains valid and can never be counted as a pair.
single <- stage2_make_read_assignments("barcode_003.ab1", default_locus = "LSU", default_direction = "Forward")
single$Isolate <- "FB121"
single <- stage2_sync_generated_names(single)
single_summary <- stage2_architecture_summary(stage2_build_architecture(single))
assert_true(single_summary$Single_read_loci == 1L && single_summary$Reads == 1L, "Single-read locus is not fully valid.")
assert_true(single_summary$Paired_loci == 0L, "One read was incorrectly reported as a Forward/Reverse pair.")

# Duplicate generated FASTA names are invalid even though source barcodes differ.
duplicate_identity <- paired
duplicate_identity$Direction <- "Forward"
duplicate_identity <- stage2_sync_generated_names(duplicate_identity)
assert_true(grepl("unique", stage2_identity_error(duplicate_identity), ignore.case = TRUE), "Duplicate generated names were not blocked.")

# Schema-1 migration preserves established evidence and legacy output names without parsing them.
legacy_state <- list(
  results = list(raw_001 = list(seq = "ACGT", evidence_marker = "keep-me")),
  summary = data.frame(sample_id = "raw_001", stringsAsFactors = FALSE),
  rename = data.frame(Original_name = "raw_001", New_name = "FB001", stringsAsFactors = FALSE),
  settings = list(target = "ITS", sequencing_primer = "Forward", forward_primer = "ITS1F"),
  blast_hits = data.frame(accession = "NR_000001.1", stringsAsFactors = FALSE)
)
migrated <- stage2_migrate_v1_state(legacy_state)
assert_true(identical(migrated$results$raw_001$evidence_marker, "keep-me"), "Migration changed established read evidence.")
assert_true(identical(migrated$blast_hits$accession, "NR_000001.1"), "Migration changed BLAST evidence.")
assert_true(migrated$read_assignments$Source_ID == "raw_001", "Migration lost original sample identity.")
assert_true(migrated$read_assignments$Final_Name == "FB001", "Migration changed an existing legacy output name.")

# Duplicate upload basenames remain unsafe because the established store is keyed by Source_ID.
duplicate_source <- stage2_make_read_assignments(c("folder_a/barcode_001.ab1", "folder_b/barcode_001.ab1"), default_locus = "ITS")
assert_true(grepl("duplicate AB1 basenames", stage2_validate_assignments(duplicate_source), fixed = TRUE), "Duplicate source IDs were not blocked.")

cat("v3.0.0-alpha.6 Stage 2 architecture tests passed.\n")
