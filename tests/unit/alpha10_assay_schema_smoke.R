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

other_state <- legacy_state
other_state$settings$target <- "Other"
other_error <- tryCatch(assay_migrate_schema5_state(other_state), error = function(e) conditionMessage(e))
assert_true(is.character(other_error) && grepl("controlled PITAX locus vocabulary", other_error, fixed = TRUE),
            "Schema-5 migration must reject free-text Other instead of remapping to ITS.")

unknown_state <- legacy_state
unknown_state$settings$target <- "SomethingWeird"
unknown_error <- tryCatch(assay_migrate_schema5_state(unknown_state), error = function(e) conditionMessage(e))
assert_true(is.character(unknown_error) && grepl("controlled PITAX locus vocabulary", unknown_error, fixed = TRUE),
            "Schema-5 migration must reject unknown loci instead of remapping to ITS.")

editor_from_profile <- function(row, amplicon = row$Expected_Amplicon_Length[1], max_pos = row$Maximum_Sequence_Position[1]) {
  list(
    assay_name = row$Assay_Name[1],
    target = row$Locus_ID[1],
    forward_primer = row$Forward_Primer_Name[1],
    reverse_primer = row$Reverse_Primer_Name[1],
    forward_primer_seq = row$Forward_Primer_Sequence[1],
    reverse_primer_seq = row$Reverse_Primer_Sequence[1],
    expected_amplicon_len = amplicon,
    absolute_max_base_index = max_pos
  )
}

simulate_commit_active_assay <- function(profiles, assay_id, editor, direction = "Forward") {
  result <- assay_try_apply_editor_inputs(profiles, assay_id, editor)
  if (!isTRUE(result$committed)) return(result)
  idx <- match(result$assay_id, result$profiles$Assay_ID)
  result$settings <- assay_resolve_read_settings(result$profiles[idx, , drop = FALSE], NULL, direction)
  result
}

amplicon_change <- tryCatch(
  simulate_commit_active_assay(its, "assay-its", editor_from_profile(its, amplicon = 720L, max_pos = 680L)),
  error = function(e) e
)
assert_true(!inherits(amplicon_change, "error"), "Changing amplicon length must not crash assay commit/resolve.")
assert_true(isTRUE(amplicon_change$committed), "A complete amplicon-length edit was not committed.")
assert_true(identical(amplicon_change$profiles$Expected_Amplicon_Length[1], 720L), "Committed amplicon length was not stored.")
assert_true(identical(amplicon_change$settings$expected_amplicon_len, 720L), "Resolved settings did not pick up the new amplicon length.")

for (draft_max in list(NULL, NA_real_, NA_integer_, "", "  ", numeric(0))) {
  skipped <- tryCatch(
    simulate_commit_active_assay(its, "assay-its", editor_from_profile(its, amplicon = 720L, max_pos = draft_max)),
    error = function(e) e
  )
  assert_true(!inherits(skipped, "error"), "Transient missing max position must not throw.")
  assert_true(identical(skipped$committed, FALSE), "Transient missing max position must not commit.")
  assert_true(identical(skipped$reason, "incomplete"), "Transient missing max position must be treated as an incomplete draft.")
  assert_true(identical(skipped$profiles$Expected_Amplicon_Length[1], 650L), "Incomplete drafts must leave the stored amplicon length unchanged.")
  assert_true(identical(skipped$profiles$Maximum_Sequence_Position[1], 680L), "Incomplete drafts must leave the stored max position unchanged.")
}

rejected_49 <- tryCatch(
  simulate_commit_active_assay(its, "assay-its", editor_from_profile(its, amplicon = 650L, max_pos = 49L)),
  error = function(e) e
)
assert_true(!inherits(rejected_49, "error"), "A real max position below 50 must be rejected without crashing Shiny commit.")
assert_true(identical(rejected_49$committed, FALSE) && identical(rejected_49$reason, "invalid"), "A max position of 49 must not be committed.")
assert_true(is.character(rejected_49$error) && grepl("at least 50 bp", rejected_49$error, fixed = TRUE), "A max position of 49 must keep the strict validation message.")
assert_true(identical(rejected_49$profiles$Maximum_Sequence_Position[1], 680L), "Rejected max position 49 must not overwrite the stored profile.")
direct_49 <- its
direct_49$Maximum_Sequence_Position[1] <- 49L
resolve_49 <- tryCatch(assay_resolve_read_settings(direct_49, NULL, "Forward"), error = function(e) conditionMessage(e))
assert_true(is.character(resolve_49) && grepl("at least 50 bp", resolve_49, fixed = TRUE),
            "assay_resolve_read_settings must still reject a real committed max position below 50.")

accepted_50 <- tryCatch(
  simulate_commit_active_assay(its, "assay-its", editor_from_profile(its, amplicon = 650L, max_pos = 50L)),
  error = function(e) e
)
assert_true(!inherits(accepted_50, "error"), "A max position of 50 must not crash assay commit/resolve.")
assert_true(isTRUE(accepted_50$committed), "A max position of 50 must be accepted.")
assert_true(identical(accepted_50$profiles$Maximum_Sequence_Position[1], 50L), "Committed max position 50 was not stored.")
assert_true(identical(accepted_50$settings$absolute_max_base_index, 50L), "Resolved settings did not pick up max position 50.")

upload_text <- paste(readLines(file.path("R", "server", "stages", "20_upload.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
assert_true(grepl("assay_try_apply_editor_inputs(", upload_text, fixed = TRUE),
            "commit_active_assay_from_inputs must apply editor inputs through assay_try_apply_editor_inputs.")
assert_true(!grepl("Maximum_Sequence_Position[idx] <- 680L", upload_text, fixed = TRUE),
            "Assay commit must not silently replace incomplete max position drafts with 680.")

cat("Assay/schema-6 foundation tests passed.\n")
