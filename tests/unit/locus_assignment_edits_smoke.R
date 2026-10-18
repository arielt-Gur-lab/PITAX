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
cat("Project locus and assignment batch behavior tests passed.\n")

