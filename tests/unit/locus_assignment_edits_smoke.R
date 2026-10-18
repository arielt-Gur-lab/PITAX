# Project-local locus catalog and literal assignment text transformations.
source(file.path("R", "domain", "assay", "assay_profiles.R"))
source(file.path("R", "domain", "assignment", "stage2_architecture.R"))

check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
rejects <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  check(is.character(error) && grepl(pattern, error, fixed = TRUE), paste("Expected rejection:", pattern))
}

defaults <- pitax_locus_vocabulary()
loci <- pitax_save_locus(defaults, "sdhc", "SDH subunit C")
loci <- pitax_save_locus(loci, "SDHD", "SDH subunit D")
check(all(c("SDHC", "SDHD") %in% loci$Locus_ID), "Custom loci were not registered.")
check(identical(pitax_locus_vocabulary(), defaults), "A project edit leaked into built-in/global state.")
check(unname(pitax_locus_choices(loci)["SDH subunit C"]) == "SDHC", "Dropdown choices lost the custom ID.")
rejects(pitax_save_locus(loci, "sdhc", "Duplicate"), "already exists")
rejects(pitax_save_locus(loci, "EF1-alpha", "Duplicate alias"), "already exists")
rejects(pitax_save_locus(loci, "SDH C", "Malformed ID"), "Locus IDs")
rejects(pitax_save_locus(loci, "Other", "Other"), "Locus IDs")
rejects(pitax_save_locus(loci, "NEWGENE", ""), "display name")
rejects(pitax_save_locus(loci, "NEWGENE", "SDH subunit C"), "display names must be unique")
rejects(pitax_save_locus(loci, "NEWID", "Renamed", existing_id = "SDHC"), "ID is stable")

loci <- pitax_save_locus(loci, "SDHC", "Succinate dehydrogenase C", existing_id = "SDHC")
loci <- pitax_save_locus(loci, "ITS", "Internal transcribed spacer", existing_id = "ITS")
profile <- assay_apply_locus_labels(assay_profile_from_legacy_settings(list(target = "SDHC")), loci)
check(is.null(assay_validate_profiles(profile, loci)), "Registered custom profile was rejected.")
check(grepl("Succinate", assay_coerce_profiles(profile)$Locus_Display_Name), "Coercion erased the edited label.")
check(assay_resolve_read_settings(profile)$target == "SDHC", "Processing settings remapped a custom locus.")
check(!is.null(assay_validate_profiles(profile, defaults)), "An unregistered profile was accepted against a catalog.")
edited <- assay_try_apply_editor_inputs(profile, profile$Assay_ID, list(
  target = "SDHD", assay_name = "Test assay", expected_amplicon_len = 600L,
  absolute_max_base_index = 700L
), loci)
check(isTRUE(edited$committed) && edited$profiles$Locus_ID == "SDHD", "Editor lost a custom locus.")
invalid <- assay_try_apply_editor_inputs(profile, profile$Assay_ID, list(
  target = "NOT_REGISTERED", expected_amplicon_len = 600L, absolute_max_base_index = 700L
), loci)
check(!isTRUE(invalid$committed) && identical(invalid$profiles, profile), "Invalid editor input silently changed the assay.")

reads <- stage2_make_read_assignments(c("barcode.1.ab1", "barcode.2.ab1", "barcode.3.ab1"),
                                     default_assay_id = profile$Assay_ID, default_direction = "Forward")
reads$Isolate <- c("FB.1", "FB.2", "FB.3")
original <- reads
batch <- stage2_apply_batch_edits(reads, c(1L, 3L), ".", "-", prefix = "X", suffix = "Y")
check(identical(batch$assignments$Isolate, c("XFB-1Y", "FB.2", "XFB-3Y")), "Literal replacement or row scope failed.")
check(batch$matched == 2L, "Match count is incorrect.")
check(identical(batch$assignments$Source_ID, original$Source_ID), "Batch edits changed immutable source IDs.")
check(identical(batch$assignments$Read_ID, original$Read_ID), "Batch edits changed read keys.")
removed <- stage2_apply_batch_edits(reads, 1L, ".", "")
check(removed$assignments$Isolate[1] == "FB1", "An empty replacement must delete the literal match.")
reads$Isolate[1] <- "FB 1"
spaces <- stage2_apply_batch_edits(reads, 1L, " ", "")
check(spaces$assignments$Isolate[1] == "FB1" && spaces$matched == 1L, "Literal whitespace searches were discarded.")
reads <- original
miss <- stage2_apply_batch_edits(reads, 1:3, "MISSING", "X")
check(miss$matched == 0L && identical(miss$assignments, reads), "A missing match was reported/treated as a change.")
reads$Isolate <- rep("", 3)
empty <- stage2_apply_batch_edits(reads, 1:3, "barcode.", "FB")
check(empty$matched == 0L && all(empty$assignments$Isolate == ""), "Barcode identity was inferred without explicit selection.")
copied <- stage2_apply_batch_edits(reads, 1:3, "barcode.", "FB", text_source = "source")
check(identical(copied$assignments$Isolate, c("FB1", "FB2", "FB3")), "Explicit barcode-to-isolate replacement failed.")
synced <- stage2_sync_generated_names(copied$assignments, assay_profiles = profile)
check(all(synced$Locus == "SDHC") && synced$Final_Name[1] == "FB1_SDHC_F", "Custom locus/name propagation failed.")
check(is.null(stage2_validate_assignments(synced, assay_profiles = profile)), "Custom read assignments failed validation.")
architecture <- stage2_build_architecture(synced, assay_profiles = profile)
check(all(architecture$loci$Locus == "SDHC"), "Architecture dropped the custom locus.")
rejects(stage2_apply_batch_edits(reads, 0L), "Invalid assignment")
rejects(stage2_apply_batch_edits(reads, 1L, text_source = "final_name"), "Choose isolate codes")

path <- tempfile(fileext = ".sangerproject")
saveRDS(list(locus_vocabulary = loci, assay_profiles = profile, read_assignments = synced), path)
restored <- readRDS(path)
unlink(path)
check(identical(pitax_coerce_locus_vocabulary(restored$locus_vocabulary), loci), "Saved catalog changed on readback.")
check(identical(pitax_coerce_locus_vocabulary(NULL), defaults), "Legacy schema-6 catalog default changed.")
check("SDHD" %in% restored$locus_vocabulary$Locus_ID, "An unused custom locus was not preserved.")
# Profile-specific naming changes labels, never the canonical biological locus.
source(file.path("R", "domain", "consensus", "stage3_consensus.R"))
source(file.path("R", "domain", "multilocus", "stage4_multilocus.R"))
source(file.path("R", "export", "export_tools.R"))
source(file.path("R", "export", "availability.R"))
first <- assay_profile_from_legacy_settings(list(target = "ITS", assay_name = "Assay One"), "assay-one")
second <- assay_profile_from_legacy_settings(list(target = "LSU", assay_name = "Assay Two"), "assay-two")
first$Name_Source <- "assay"
profiles <- rbind(first, second)
named <- stage2_make_read_assignments(c("F.ab1", "R.ab1"), default_assay_id = "assay-one")
named$Isolate <- "FB1"; named$Direction <- c("Forward", "Reverse")
named <- stage2_sync_generated_names(named, assay_profiles = profiles)
check(identical(named$Final_Name, c("FB1_Assay-One_F", "FB1_Assay-One_R")), "Assay display name was not used.")
check(all(named$Locus == "ITS"), "Naming choice changed the biological locus.")
check(is.null(stage2_identity_error(named)) && is.null(stage2_validate_assignments(named, assay_profiles = profiles)), "Assay-based names failed validation.")
legacy_profiles <- profiles; legacy_profiles$Name_Source <- NULL
check(all(assay_coerce_profiles(legacy_profiles)$Name_Source == "locus"), "Old profiles did not default to locus names.")
invalid_profile <- profiles; invalid_profile$Name_Source[1] <- "invalid"
check(!is.null(assay_validate_profiles(invalid_profile)), "Invalid naming mode was accepted.")
check(assay_resolve_read_settings(first)$assay_name_source == "assay", "Processing metadata lost the naming choice.")

result <- function() list(seq = "ACGTACGT", curation = list(trim_start = 1L, trim_end = 8L, revision = 0L),
  summary = data.frame(trim_start = 1L, trim_end = 8L),
  ab1_evidence = list(detail = data.frame(Basecaller_quality = rep(30, 8))),
  metrics = data.frame(peak_ratio = rep(3, 8)))
results <- list(F = result(), R = result())
paired <- stage3_build_consensus_set(named, results, min_overlap = 4L)
check(paired$records[[1]]$final_name == "FB1_Assay-One" && paired$records[[1]]$locus == "ITS", "Paired naming lost the selected label or biological locus.")
simple <- stage3_build_consensus_set(named, results, project_mode = "simple", min_overlap = 4L)
check(simple$records[[1]]$final_name == "FB1_Assay-One_F", "Independent analysis ignored the naming choice.")
single <- stage3_build_consensus_set(named[1, , drop = FALSE], results, min_overlap = 4L)
check(single$records[[1]]$final_name == "FB1_Assay-One", "Single-read paired representative ignored the label.")
check(grepl(">FB1_Assay-One", make_fasta(stage3_analysis_records(paired)), fixed = TRUE), "Analysis FASTA ignored the chosen name.")
ambiguous <- named; ambiguous$Name_Label[2] <- "Different assay"
rejects(stage3_build_consensus_set(ambiguous, results, min_overlap = 4L), "different naming labels")

mixed <- named; mixed$Isolate <- c("FB1", "FB2"); mixed$Assay_ID[2] <- "assay-two"
mixed <- stage2_sync_generated_names(mixed, assay_profiles = profiles)
check(identical(mixed$Final_Name, c("FB1_Assay-One_F", "FB2_LSU_R")), "Naming mode leaked between profiles.")
profiles$Name_Source[1] <- "locus"
locus_named <- stage2_sync_generated_names(named, assay_profiles = profiles)
check(identical(locus_named$Final_Name, c("FB1_ITS_F", "FB1_ITS_R")), "Switching back did not restore locus naming.")
profiles$Name_Source <- "assay"; profiles$Assay_Name <- "Same Name"
collision <- stage2_sync_generated_names(mixed, assay_profiles = profiles)
collision$Isolate <- "FB1"; collision$Direction <- "Forward"
collision <- stage2_sync_generated_names(collision, assay_profiles = profiles)
check(grepl("unique", stage2_identity_error(collision)), "Cleaned naming collisions were accepted.")

project <- list(format = "SangerSequencePipelineProject", schema_version = 6L, app_version = "3.5.1",
  state = list(results = results, consensus_set = paired))
evidence <- stage4_extract_project_evidence(project, "Test project")
check(evidence$Locus == "ITS" && evidence$Name_Label == "Assay One", "Multi-locus evidence mixed the label with biological identity.")
profile <- stage4_empty_profile(); profile$evidence <- evidence
check(grepl(">FB1_Assay-One", stage4_make_fasta(profile), fixed = TRUE), "Multi-locus FASTA ignored the chosen label.")
snap <- list(multilocus_profile = profile, read_assignments = named)
exported <- pitax_export_extract_multilocus_sequences(snap)
check(exported[[1]]$final_name == "FB1_Assay-One" && exported[[1]]$locus == "ITS", "Shared export lost name choice or locus identity.")

cat("Project locus and assignment batch behavior tests passed.\n")

