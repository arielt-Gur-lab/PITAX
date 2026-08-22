# PITAX v3.0.0-alpha.9.2 - Stage 4 multi-locus profile tests.

source(file.path("R", "domain", "consensus", "stage3_consensus.R"))
source(file.path("R", "export", "export_tools.R"))
source(file.path("R", "domain", "multilocus", "stage4_multilocus.R"))

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)

make_result <- function(sequence) {
  n <- nchar(sequence)
  list(
    seq = sequence,
    curation = list(trim_start = 1L, trim_end = n, revision = 0L),
    summary = data.frame(trim_start = 1L, trim_end = n, stringsAsFactors = FALSE),
    ab1_evidence = list(detail = data.frame(Basecaller_quality = rep(30, n))),
    metrics = data.frame(peak_ratio = rep(3, n), stringsAsFactors = FALSE)
  )
}

make_project <- function(locus, isolate = "ISO1", identification = "", level = "", genus = "", accession = "") {
  source_id <- paste0(isolate, "_", locus, "_F")
  assignments <- data.frame(
    Source_ID = source_id, Final_Name = source_id, Isolate = isolate,
    Locus = locus, Direction = "Forward", stringsAsFactors = FALSE
  )
  results <- setNames(list(make_result("ACGTACGTACGT")), source_id)
  consensus <- stage3_build_consensus_set(assignments, results, min_overlap = 4)
  taxonomy <- data.frame()
  if (nzchar(identification)) {
    taxonomy <- data.frame(
      original_name = "consensus_001", final_name = paste(isolate, locus, sep = "_"),
      recommended_identification = identification, recommended_level = level,
      confidence = "High", best_match_genus = genus,
      best_molecular_match = identification, best_match_accession = accession,
      best_match_identity_percent = 99.5, best_match_query_coverage_percent = 100,
      reference_support = "Curated reference context", locus_discrimination = "Good",
      rid = paste0("RID_", locus), analyzed_at = "2026-08-20 10:00:00",
      stringsAsFactors = FALSE
    )
  }
  list(
    format = "SangerSequencePipelineProject", schema_version = 5L,
    app_version = "3.0.0-alpha.9", saved_at = "2026-08-20 10:00:00",
    state = list(results = results, consensus_set = consensus, taxonomy_summary = taxonomy)
  )
}

its <- make_project("ITS", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "NR_001")
tef1 <- make_project("TEF1", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "ON_002")
profile <- stage4_build_profile(list(its, tef1), c("ITS project", "TEF1 project"), c("md5-a", "md5-b"))

assert_true(is.null(stage4_profile_gate_error(profile)), "A valid two-locus profile did not pass the structural gate.")
assert_true(nrow(profile$evidence) == 2L && nrow(profile$profiles) == 1L, "The two locus rows were not combined into one isolate profile.")
assert_true(profile$profiles$Profile_Status[1] == "CONCORDANT_SPECIES", "Concordant species evidence was not recognized.")
assert_true(profile$profiles$Supported_Species[1] == "Fusarium oxysporum", "The concordant species call was not retained.")
assert_true(all(c("Source", "Source_MD5", "Consensus_Revision", "Best_Match_Accession", "RID") %in% names(profile$evidence)), "Required per-locus provenance fields are missing.")

selected_profile <- stage4_isolate_profile(profile, "iso1")
selected_evidence <- stage4_isolate_evidence(profile, "ISO1")
overview <- stage4_profile_overview(profile)
assert_true(nrow(selected_profile) == 1L && selected_profile$Isolate[1] == "ISO1", "The isolate profile selector did not resolve its profile row.")
assert_true(nrow(selected_evidence) == 2L && all(c("ITS", "TEF1") %in% selected_evidence$Locus), "The isolate evidence selector did not retain all loci.")
assert_true(overview$Isolates[1] == 1L && overview$Loci[1] == 2L && overview$Concordant[1] == 1L, "The visual overview counts are incorrect.")

# The visual selector must keep isolates separate while showing every locus for each one.
its_iso2 <- make_project("ITS", isolate = "ISO2", identification = "Aspergillus niger", level = "species", genus = "Aspergillus", accession = "NR_010")
tef1_iso2 <- make_project("TEF1", isolate = "ISO2", identification = "Aspergillus niger", level = "species", genus = "Aspergillus", accession = "ON_011")
two_isolates <- stage4_build_profile(
  list(its, tef1, its_iso2, tef1_iso2),
  c("ISO1 ITS", "ISO1 TEF1", "ISO2 ITS", "ISO2 TEF1"),
  c("md5-j", "md5-k", "md5-l", "md5-m")
)
two_overview <- stage4_profile_overview(two_isolates)
assert_true(two_overview$Isolates[1] == 2L && two_overview$Loci[1] == 2L, "The overview confused isolate count with locus count.")
assert_true(nrow(stage4_isolate_evidence(two_isolates, "ISO2")) == 2L, "The second isolate did not expose both of its loci.")
assert_true(all(stage4_isolate_evidence(two_isolates, "ISO2")$Isolate == "ISO2"), "Evidence from different isolates was mixed in the visual selector.")

# A 2:1 numerical majority must not erase a cross-genus conflict.
rpb2_same <- make_project("RPB2", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "ON_003")
bt_conflict <- make_project("BT2", identification = "Aspergillus niger", level = "species", genus = "Aspergillus", accession = "ON_004")
conflict <- stage4_build_profile(
  list(its, rpb2_same, bt_conflict), c("ITS", "RPB2", "BT2"), c("md5-c", "md5-d", "md5-e")
)
assert_true(conflict$profiles$Profile_Status[1] == "GENUS_CONFLICT", "A 2:1 locus majority incorrectly voted away a genus conflict.")
assert_true(!nzchar(conflict$profiles$Supported_Genus[1]), "A combined genus was reported despite conflicting loci.")

# The same isolate/locus cannot be imported twice, even with different source names.
duplicate_blocked <- tryCatch({
  stage4_build_profile(list(its, its), c("ITS A", "ITS B"), c("md5-f", "md5-g"))
  FALSE
}, error = function(e) grepl("Duplicate Isolate/Locus", conditionMessage(e), fixed = TRUE))
assert_true(duplicate_blocked, "Duplicate Isolate/Locus evidence was not blocked.")

# Stage 4 can retain sequences before taxonomy, but it reports the missing evidence.
raw_tef1 <- make_project("TEF1")
partial <- stage4_build_profile(list(its, raw_tef1), c("ITS", "TEF1 raw"), c("md5-h", "md5-i"))
assert_true(partial$profiles$Profile_Status[1] == "PARTIAL_EVIDENCE", "Missing per-locus taxonomy was not reported as partial evidence.")
assert_true(partial$profiles$Taxonomy_Complete[1] == "1/2", "Taxonomy completeness was not counted correctly.")

# Current-session binding detects later sequence or taxonomy changes.
assert_true(stage4_current_project_matches(profile, its, "ITS project"), "An unchanged source project was marked stale.")
changed <- its
changed$state$results[[1]]$seq <- "ACGTACGTACGA"
assert_true(!stage4_current_project_matches(profile, changed, "ITS project"), "A changed source sequence did not stale its profile snapshot.")

fasta <- stage4_make_fasta(profile)
assert_true(length(gregexpr(">", fasta, fixed = TRUE)[[1]]) == 2L, "The multi-locus FASTA did not contain one record per locus.")

legacy <- list(results = list(keep = TRUE), migration_log = "Stage 3 retained.")
migrated <- stage4_migrate_v4_state(legacy)
assert_true(isTRUE(migrated$results$keep), "Stage 4 migration changed existing project evidence.")
assert_true(identical(migrated$multilocus_profile$schema, "pitax-multilocus-profile-v1"), "Stage 4 migration did not initialize its store.")

cat("v3.0.0-alpha.9.2 Stage 4 multi-locus tests passed.\n")
