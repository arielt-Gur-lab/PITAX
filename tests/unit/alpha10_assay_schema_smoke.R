# PITAX - schema 6 assay/locus foundation tests.

source(file.path("R", "domain", "assay", "assay_profiles.R"))
source(file.path("R", "domain", "assignment", "stage2_architecture.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

vocabulary <- pitax_locus_vocabulary()
assert_true(all(c("Locus_ID", "Display_Name") %in% names(vocabulary)), "Controlled locus vocabulary columns are missing.")
assert_true(!"Other" %in% vocabulary$Locus_ID, "Free-text Other locus must not be available in schema 6.")
assert_true(identical(pitax_normalize_locus_id("TEF1 / EF1-alpha"), "TEF1"), "Legacy TEF1 label was not normalized.")
assert_true(identical(pitax_normalize_locus_id("Beta-tubulin"), "TUB2"), "Legacy beta-tubulin label was not normalized.")

its <- assay_profile_from_legacy_settings(list(
  target = "ITS", forward_primer = "ITS1F", reverse_primer = "ITS4",
  expected_amplicon_len = 650L, absolute_max_base_index = 680L
), assay_id = "assay-its")
tef <- assay_profile_from_legacy_settings(list(
  target = "TEF1 / EF1-alpha", forward_primer = "EF1-728F", reverse_primer = "EF1-986R",
  expected_amplicon_len = 350L, absolute_max_base_index = 500L
), assay_id = "assay-tef1")
profiles <- rbind(its, tef)
assert_true(is.null(assay_validate_profiles(profiles)), "Valid multi-assay profiles were rejected.")

reads <- stage2_make_read_assignments(
  c("FB120_1.ab1", "FB120_2.ab1", "FB121_1.ab1", "FB121_2.ab1"),
  default_assay_id = "assay-its"
)
reads$Isolate <- c("FB120", "FB120", "FB121", "FB121")
reads$Assay_ID <- c("assay-its", "assay-its", "assay-tef1", "assay-tef1")
reads$Direction <- c("Forward", "Reverse", "Forward", "Reverse")
reads <- stage2_sync_generated_names(reads, assay_profiles = profiles)
assert_true(identical(reads$Locus, c("ITS", "ITS", "TEF1", "TEF1")), "Read loci were not inherited from assay profiles.")
assert_true(identical(reads$Primer, c("ITS1F", "ITS4", "EF1-728F", "EF1-986R")), "Direction-specific primer provenance was not inherited from assay profiles.")
assert_true(is.null(stage2_validate_assignments(reads, assay_profiles = profiles)), "Valid assay-linked read assignments were rejected.")
architecture <- stage2_build_architecture(reads, project_id = "alpha10-test", assay_profiles = profiles)
summary <- stage2_architecture_summary(architecture)
assert_true(identical(architecture$schema, "pitax-project-architecture-v2"), "Schema 6 architecture marker is missing.")
assert_true(nrow(architecture$assays) == 2L && summary$Isolates == 2L && summary$Loci == 2L, "Multi-assay architecture counts are incorrect.")

bad_reads <- reads
bad_reads$Assay_ID[1] <- "missing-assay"
assert_true(grepl("does not exist", stage2_validate_assignments(bad_reads, assay_profiles = profiles), fixed = TRUE), "Dangling assay references were not blocked.")

legacy_assignments <- reads[1:2, , drop = FALSE]
legacy_assignments$Assay_ID <- NULL
legacy_assignments$Locus <- "ITS"
legacy_state <- list(
  settings = list(
    target = "ITS", forward_primer = "ITS1F", reverse_primer = "ITS4",
    expected_amplicon_len = 650L, absolute_max_base_index = 680L,
    window = 25L, min_peak_ratio = 3
  ),
  read_assignments = legacy_assignments,
  results = list(FB120_1 = list(seq = "ACGT", evidence_marker = "preserve-me")),
  blast_hits = data.frame(accession = "NR_000001.1", stringsAsFactors = FALSE)
)
migrated <- assay_migrate_schema5_state(legacy_state)
assert_true(nrow(migrated$assay_profiles) == 1L, "Schema-5 settings did not become exactly one assay profile.")
assert_true(all(migrated$read_assignments$Assay_ID == migrated$assay_profiles$Assay_ID[1]), "Schema-5 reads were not linked to the migrated assay.")
assert_true(identical(migrated$results$FB120_1$evidence_marker, "preserve-me"), "Schema-5 migration changed read evidence.")
assert_true(identical(migrated$blast_hits$accession, "NR_000001.1"), "Schema-5 migration changed BLAST evidence.")

cat("Assay/schema-6 foundation tests passed.\n")
