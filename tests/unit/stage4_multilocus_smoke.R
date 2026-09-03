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

# Current-session TEF1 must still join when taxonomy was stored under the source read id
# or final name instead of consensus_001. This is the FB83 / Current session miss.
fb83_lsu <- make_project("LSU", isolate = "FB83", identification = "Pleurotus", level = "genus", genus = "Pleurotus", accession = "KX787096.1")
fb83_tef <- make_project("TEF1", isolate = "FB83", identification = "Pleurotus eryngii", level = "species", genus = "Pleurotus", accession = "TEF_001")
source_id <- names(fb83_tef$state$results)[1]
fb83_tef_by_read <- fb83_tef
fb83_tef_by_read$state$taxonomy_summary$original_name <- source_id
fb83_tef_by_read$state$taxonomy_summary$final_name <- paste0("FB83_TEF1_F")
joined_by_read <- stage4_extract_project_evidence(fb83_tef_by_read, "Current session")
assert_true(identical(joined_by_read$Taxonomy_Status[1], "Analyzed"),
            "TEF1 taxonomy keyed by source read id was not attached to the current-session locus card.")
assert_true(identical(joined_by_read$Recommended_Identification[1], "Pleurotus eryngii"),
            "TEF1 identification was not displayed after joining by source read id.")

fb83_tef_by_name <- fb83_tef
fb83_tef_by_name$state$taxonomy_summary$original_name <- "stale_consensus_999"
fb83_tef_by_name$state$taxonomy_summary$final_name <- "FB83_TEF1_Forward"
fb83_tef_by_name$state$taxonomy_summary$target <- "TEF1"
joined_by_isolate <- stage4_extract_project_evidence(fb83_tef_by_name, "Current session")
assert_true(identical(joined_by_isolate$Taxonomy_Status[1], "Analyzed"),
            "TEF1 taxonomy keyed by isolate+locus was not attached to the current-session locus card.")

fb83_profile <- stage4_build_profile(
  list(fb83_tef_by_read, fb83_lsu),
  c("Current session", "LSU-2.9.26.sangerproject"),
  c("", "md5-lsu")
)
fb83_evidence <- stage4_isolate_evidence(fb83_profile, "FB83")
assert_true(nrow(fb83_evidence) == 2L && all(c("TEF1", "LSU") %in% fb83_evidence$Locus),
            "FB83 did not keep both current-session TEF1 and imported LSU evidence.")
assert_true(identical(as.character(fb83_evidence$Taxonomy_Status[fb83_evidence$Locus == "TEF1"])[1], "Analyzed"),
            "FB83 TEF1 from the current session was shown as not analyzed despite taxonomy evidence.")

# BLAST-only TEF1 should still plot identity/coverage even before taxonomy is run.
blast_only <- make_project("TEF1", isolate = "FB83")
blast_only$state$blast_hits <- data.frame(
  original_name = source_id, final_name = "FB83_TEF1", rid = "RIDBLAST1",
  organism = "Pleurotus ostreatus", accession = "BLASTACC", rank = 1,
  identity_percent = 99.1, query_coverage_percent = 98.4,
  stringsAsFactors = FALSE
)
blast_ev <- stage4_extract_project_evidence(blast_only, "Current session")
assert_true(identical(blast_ev$Taxonomy_Status[1], "Not analyzed"),
            "BLAST-only TEF1 must not be marked as taxonomically analyzed.")
assert_true(isTRUE(abs(blast_ev$Best_Match_Identity[1] - 99.1) < 0.01),
            "BLAST-only TEF1 identity was not displayed on the locus card.")
assert_true(identical(blast_ev$RID[1], "RIDBLAST1"),
            "BLAST-only TEF1 RID was not displayed on the locus card.")

# An unresolved taxonomic call is still analyzed. The card must not say "Not analyzed".
fb50_tef <- make_project("TEF1", isolate = "FB50", identification = "Pleurotus", level = "genus", genus = "Pleurotus", accession = "OZ415662.1")
fb50_lsu <- make_project("LSU", isolate = "FB50", identification = "Unresolved", level = "unresolved", genus = "Pleurotus", accession = "PQ652238.1")
fb50_lsu$state$taxonomy_summary$confidence <- "Low / review"
fb50_lsu$state$taxonomy_summary$locus_discrimination <- "Poor at genus and species level"
fb50_lsu$state$taxonomy_summary$best_molecular_match <- "Pleurotus pulmonarius"
fb50_unresolved <- stage4_build_profile(
  list(fb50_tef, fb50_lsu),
  c("Current session", "LSU-2.9.26.sangerproject"),
  c("", "md5-fb50-lsu")
)
fb50_ev <- stage4_isolate_evidence(fb50_unresolved, "FB50")
lsu_row <- fb50_ev[fb50_ev$Locus == "LSU", , drop = FALSE]
tef_row <- fb50_ev[fb50_ev$Locus == "TEF1", , drop = FALSE]
assert_true(identical(lsu_row$Taxonomy_Status[1], "Analyzed"), "Unresolved LSU must remain taxonomically analyzed.")
assert_true(identical(lsu_row$Recommended_Identification[1], "Unresolved"),
            "Unresolved LSU identification was blanked instead of being retained.")
assert_true(identical(stage4_display_call(lsu_row$Taxonomy_Status[1], lsu_row$Recommended_Identification[1]), "Unresolved"),
            "Analyzed + Unresolved must display Unresolved, not Not analyzed.")
assert_true(identical(stage4_display_call(tef_row$Taxonomy_Status[1], tef_row$Recommended_Identification[1]), "Pleurotus"),
            "A supported genus call must still display the taxon name.")
assert_true(!nzchar(lsu_row$Supported_Genus[1]),
            "An unresolved LSU call must not donate a supported genus from the BLAST best match.")
fb50_prof <- stage4_isolate_profile(fb50_unresolved, "FB50")
assert_true(identical(fb50_prof$Taxonomy_Complete[1], "2/2"), "Unresolved LSU must still count as interpreted.")
assert_true(identical(fb50_prof$Profile_Status[1], "PARTIAL_EVIDENCE"),
            "Genus + unresolved LSU must not be treated as a combined call.")
assert_true(!grepl("Complete missing locus-level BLAST", fb50_prof$Next_Action[1], fixed = TRUE),
            "2/2 interpreted profiles must not ask the user to complete missing BLAST/taxonomy.")
assert_true(grepl("unresolved", fb50_prof$Next_Action[1], ignore.case = TRUE),
            "The next action must tell the user to review the unresolved locus.")
assert_true(identical(stage4_display_call("Not analyzed", ""), "Not analyzed"),
            "A locus with no taxonomy row must still display Not analyzed.")

# --- N loci (3+): no two-source ceiling ---------------------------------
lsu3 <- make_project("LSU", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "LSU_001")
three_sources <- stage4_build_profile(
  list(its, tef1, lsu3),
  c("ITS project", "TEF1 project", "LSU project"),
  c("md5-n1", "md5-n2", "md5-n3")
)
assert_true(is.null(stage4_profile_gate_error(three_sources)), "A three-locus profile must pass the structural gate.")
assert_true(nrow(three_sources$evidence) == 3L && nrow(three_sources$profiles) == 1L,
            "Three imported locus sources must become three evidence rows and one isolate profile.")
assert_true(all(c("ITS", "TEF1", "LSU") %in% three_sources$evidence$Locus),
            "ITS, TEF1 and LSU must all be retained in a three-source profile.")
assert_true(identical(three_sources$profiles$Profile_Status[1], "CONCORDANT_SPECIES"),
            "Three concordant species loci must still report CONCORDANT_SPECIES.")
three_overview <- stage4_profile_overview(three_sources)
assert_true(identical(as.integer(three_overview$Loci[1]), 3L), "Overview must count three distinct loci.")

make_multi_locus_session <- function(isolate, locus_specs) {
  assignments <- do.call(rbind, lapply(locus_specs, function(spec) {
    source_id <- paste0(isolate, "_", spec$locus, "_F")
    data.frame(
      Source_ID = source_id, Final_Name = source_id, Isolate = isolate,
      Locus = spec$locus, Direction = "Forward", stringsAsFactors = FALSE
    )
  }))
  results <- list()
  taxonomy_rows <- list()
  for (spec in locus_specs) {
    source_id <- paste0(isolate, "_", spec$locus, "_F")
    results[[source_id]] <- make_result("ACGTACGTACGT")
    if (!is.null(spec$identification) && nzchar(spec$identification)) {
      taxonomy_rows[[length(taxonomy_rows) + 1L]] <- data.frame(
        original_name = NA_character_, final_name = paste(isolate, spec$locus, sep = "_"),
        recommended_identification = spec$identification,
        recommended_level = if (!is.null(spec$level)) spec$level else "species",
        confidence = if (!is.null(spec$confidence)) spec$confidence else "High",
        best_match_genus = if (!is.null(spec$genus)) spec$genus else "",
        best_molecular_match = spec$identification,
        best_match_accession = if (!is.null(spec$accession)) spec$accession else "ACC",
        best_match_identity_percent = 99.5, best_match_query_coverage_percent = 100,
        reference_support = "Curated reference context", locus_discrimination = "Good",
        rid = paste0("RID_", spec$locus), analyzed_at = "2026-08-20 10:00:00",
        target = spec$locus, stringsAsFactors = FALSE
      )
    }
  }
  consensus <- stage3_build_consensus_set(assignments, results, min_overlap = 4)
  taxonomy <- if (length(taxonomy_rows)) do.call(rbind, taxonomy_rows) else data.frame()
  if (nrow(taxonomy)) {
    for (i in seq_len(nrow(taxonomy))) {
      match_id <- names(consensus$records)[
        vapply(consensus$records, function(r) identical(stage3_scalar_text(r$locus), taxonomy$target[i]), logical(1))
      ]
      if (length(match_id)) taxonomy$original_name[i] <- match_id[1]
    }
  }
  list(
    format = "SangerSequencePipelineProject", schema_version = 6L,
    app_version = "3.3.0", saved_at = "2026-09-03 15:00:00",
    state = list(results = results, consensus_set = consensus, taxonomy_summary = taxonomy)
  )
}

session3 <- make_multi_locus_session("ISO1", list(
  list(locus = "ITS", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "S_ITS"),
  list(locus = "TEF1", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "S_TEF"),
  list(locus = "LSU", identification = "Fusarium oxysporum", level = "species", genus = "Fusarium", accession = "S_LSU")
))
session_profile <- stage4_build_profile(list(session3), "Current session", "")
assert_true(is.null(stage4_profile_gate_error(session_profile)), "A single session with three loci must pass the gate.")
assert_true(nrow(session_profile$evidence) == 3L,
            "A single multi-locus session must emit one evidence row per locus without imports.")
assert_true(identical(session_profile$profiles$Profile_Status[1], "CONCORDANT_SPECIES"),
            "Three same-species loci in one session must be concordant at species.")

single_only <- stage4_build_profile(list(its), "ITS only", "md5-one")
assert_true(identical(single_only$profiles$Profile_Status[1], "SINGLE_LOCUS"), "One locus remains SINGLE_LOCUS.")
assert_true(grepl("another locus", single_only$profiles$Next_Action[1], fixed = TRUE),
            "SINGLE_LOCUS next action must ask for another locus, not imply a two-locus ceiling.")

mixed3 <- make_multi_locus_session("FB50", list(
  list(locus = "TEF1", identification = "Pleurotus", level = "genus", genus = "Pleurotus", accession = "M_TEF"),
  list(locus = "ITS", identification = "Pleurotus", level = "genus", genus = "Pleurotus", accession = "M_ITS"),
  list(locus = "LSU", identification = "Unresolved", level = "unresolved", genus = "Pleurotus", accession = "M_LSU", confidence = "Low / review")
))
mixed_profile <- stage4_build_profile(list(mixed3), "Current session", "")
mixed_row <- stage4_isolate_profile(mixed_profile, "FB50")
assert_true(identical(mixed_row$Taxonomy_Complete[1], "3/3"), "Two genus + one unresolved must still count 3/3 interpreted.")
assert_true(identical(mixed_row$Profile_Status[1], "CONCORDANT_GENUS"),
            "Two supported genus calls plus unresolved LSU must remain concordant at genus (unresolved does not donate a genus).")
assert_true(!grepl("Complete missing locus-level BLAST", mixed_row$Next_Action[1], fixed = TRUE),
            "Fully interpreted 3-locus profiles must not ask to complete missing BLAST.")

legacy <- list(results = list(keep = TRUE), migration_log = "Stage 3 retained.")
migrated <- stage4_migrate_v4_state(legacy)
assert_true(isTRUE(migrated$results$keep), "Stage 4 migration changed existing project evidence.")
assert_true(identical(migrated$multilocus_profile$schema, "pitax-multilocus-profile-v1"), "Stage 4 migration did not initialize its store.")

cat("v3.0.0-alpha.9.2 Stage 4 multi-locus tests passed.\n")
