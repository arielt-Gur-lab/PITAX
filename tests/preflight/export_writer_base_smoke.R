# Base-R export smoke regression: no CRAN packages required.
# Intended for GitHub Actions startup/CI preflight and local development.
stopifnot(file.exists("VERSION.txt"))
for (path in c("R/export/export_tools.R", "R/export/catalog.R",
               "R/export/availability.R",
               "R/export/plan.R", "R/export/writers.R", "R/export/execute.R",
               "R/export/header_templates.R")) source(path, local = FALSE)

empty18 <- as.data.frame(setNames(rep(list(character()), 18L),
                                  paste0("Col", seq_len(18L))),
                         stringsAsFactors = FALSE)
ready <- pitax_export_prepare_sheet_data(
  empty18, NULL, list(row_mode = "all", column_mode = "display", column_ids = character())
)
stopifnot(nrow(ready) == 1L, identical(names(ready), "Message"))
stopifnot(identical(pitax_export_prepare_sheet_data(
  data.frame(), NULL, list())$Message, "No rows available for this selection"))

config <- pitax_header_default_config()
config$draft$ordered_fields <- c("organism", "isolate", "locus", "direction",
                                 "molecule", "experiment", "seq_date")
config <- pitax_header_refresh_derived(config)
records <- list()
for (i in seq_len(18L)) {
  rid <- paste0("S", i)
  records[[rid]] <- list(final_name = paste0("ISO", i, "_TEF_R"),
                         sequence = "ACTGCATG", direction = "Reverse",
                         isolate = paste0("ISO", i), locus = "TEF", missing = FALSE)
  key <- pitax_header_record_key("qc_trimmed_sequences", rid)
  for (pair in list(c("organism", "En"), c("experiment", "Resistance_X"),
                    c("seq_date", "2026-11-05"))) {
    config <- pitax_header_set_value(config, key, pair[1], pair[2])
  }
}
rendered <- pitax_header_apply_to_records(records, "qc_trimmed_sequences", config)
stopifnot(length(rendered) == 18L)
out <- tempfile(fileext = ".fasta")
pitax_export_write_fasta_records(out, rendered, header_mode = "custom")
lines <- readLines(out, warn = FALSE)
stopifnot(sum(startsWith(lines, ">")) == 18L,
          identical(lines[1], paste0(
            '>ISO1_TEF_R organism=En isolate=ISO1 locus=TEF',
            ' direction=R molecule=DNA experiment=Resistance_X seq_date=2026-11-05')))

# Removed fields must not require stored values; experiment stays out unless selected.
config <- pitax_header_profile_remove_field(config, "experiment")
config <- pitax_header_profile_remove_field(config, "seq_date")
config <- pitax_header_profile_remove_field(config, "isolate")
config <- pitax_header_profile_remove_field(config, "locus")
config <- pitax_header_profile_remove_field(config, "direction")
config <- pitax_header_profile_remove_field(config, "molecule")
rendered <- pitax_header_apply_to_records(records, "qc_trimmed_sequences", config)
stopifnot(identical(rendered[["S18"]]$custom_header, "ISO18_TEF_R organism=En"))
asgs <- data.frame(Source_ID = paste0("S", seq_len(18L)),
                   Isolate = paste0("ISO", seq_len(18L)),
                   Locus = rep("TEF", 18L),
                   Direction = rep("Reverse", 18L),
                   Final_Name = paste0("ISO", seq_len(18L), "_TEF_R"),
                   stringsAsFactors = FALSE)
snapshot <- list(read_assignments = asgs, results = records,
                 consensus_set = list(records = list()),
                 multilocus_profile = list())
display <- pitax_header_active_edit_table(snapshot, config)
stopifnot(nrow(display) == 18L,
          !("experiment" %in% names(display)),
          !("seq_date" %in% names(display)),
          "organism" %in% names(display),
          identical(pitax_header_active_metadata_fields(config), "organism"))
config <- pitax_header_save_profile(config, "Without Experiment")
restored <- pitax_header_coerce_config(unserialize(serialize(config, NULL)))
stopifnot(identical(pitax_header_active_metadata_fields(restored), "organism"))

# Metadata lifecycle under the structured profile model.
lifecycle <- pitax_header_default_config()
lifecycle$draft$ordered_fields <- c("organism")
lifecycle <- pitax_header_save_profile(lifecycle, "Minimal")
step1 <- pitax_header_register_field(lifecycle, "Collection_Site")
stopifnot(isTRUE(step1$created),
          "collection_site" %in% pitax_header_active_metadata_fields(step1$config),
          grepl("collection_site={collection_site}", step1$config$draft_template, fixed = TRUE),
          !grepl("experiment=", step1$config$draft_template, fixed = TRUE),
          !grepl("collection_site", step1$config$templates[["Minimal"]], fixed = TRUE))
visible <- pitax_header_active_edit_table(snapshot, step1$config)
stopifnot("collection_site" %in% names(visible),
          !"experiment" %in% names(visible))
registered <- pitax_header_registered_field_table(step1$config)
stopifnot("Collection_Site" %in% registered$Field,
          identical(registered$Status[match("Collection_Site", registered$Field)],
                    "In active header"))
step2 <- pitax_header_register_field(step1$config, "collection_site")
stopifnot(!step2$created,
          sum(names(step2$config$custom_fields) == "collection_site") == 1L)
step2 <- pitax_header_profile_remove_field(step2$config, "collection_site")
stopifnot(!"collection_site" %in% names(pitax_header_active_edit_table(snapshot, step2)),
          "Collection_Site" %in% pitax_header_registered_field_table(step2)$Field)
step3 <- pitax_header_activate_field(step2, "collection_site")
stopifnot("collection_site" %in% names(pitax_header_active_edit_table(snapshot, step3)))
step3 <- pitax_header_activate_field(step3, "Collection_Site")
stopifnot(sum(step3$draft$ordered_fields == "collection_site") == 1L)
stopifnot(inherits(try(pitax_header_register_field(step3, "Final_Name"), silent = TRUE), "try-error"),
          inherits(try(pitax_header_register_field(step3, "1invalid"), silent = TRUE), "try-error"))
step3 <- pitax_header_set_value(step3, "read::S1", "collection_site", "Golan")
step3 <- pitax_header_save_profile(step3, "Site")
reloaded <- pitax_header_coerce_config(unserialize(serialize(step3, NULL)))
stopifnot(identical(reloaded$values$collection_site[
            match("read::S1", reloaded$values$Record_Key)], "Golan"),
          "collection_site" %in% pitax_header_active_metadata_fields(reloaded))

# Removing from active draft must not erase data; builtin experiment is not a deletable custom field.
removal <- pitax_header_profile_add_field(pitax_header_default_config(), "organism")
removal <- pitax_header_profile_add_field(removal, "experiment")
removal <- pitax_header_profile_add_field(removal, "location")
removal <- pitax_header_set_value(removal, "read::S1", "experiment", "keep")
no_experiment <- pitax_header_remove_from_draft(removal, "experiment")
stopifnot(!grepl("experiment=", no_experiment$draft_template, fixed = TRUE),
          "experiment" %in% names(no_experiment$values))
protected <- try(pitax_header_delete_custom_field(removal, "experiment"), silent = TRUE)
stopifnot(inherits(protected, "try-error"))
# location is builtin notes-folded, not custom — delete_custom_field rejects it.
deleted_try <- try(pitax_header_delete_custom_field(no_experiment, "location"), silent = TRUE)
stopifnot(inherits(deleted_try, "try-error"))
no_location <- pitax_header_remove_from_draft(no_experiment, "location")
stopifnot(!"location" %in% no_location$draft$ordered_fields,
          "organism" %in% no_location$draft$ordered_fields)

# The former plan builder assessed EVERY component, so a broken unrelated
# taxonomy table could abort even a FASTA-only export.
fake_assessor <- function(id, snap, selection) {
  if (identical(id, "taxonomy_team_summary")) {
    broken <- data.frame()
    broken$unexpected <- seq_len(18L)
  }
  pitax_export_status("available", record_count = 18L)
}
component_status <- pitax_export_assess_all(
  snapshot, list(), component_ids = c("qc_trimmed_sequences", "taxonomy_team_summary"),
  assessor = fake_assessor
)
stopifnot(identical(component_status$qc_trimmed_sequences$status, "available"),
          identical(component_status$taxonomy_team_summary$status, "blocked"),
          grepl("replacement has 18 rows, data has 0",
                component_status$taxonomy_team_summary$reason, fixed = TRUE),
          grepl("taxonomy_team_summary",
                component_status$taxonomy_team_summary$reason, fixed = TRUE))

# Exercise the actual export plan with 18 source reads and deliberately
# incomplete/unavailable non-FASTA project components.
snapshot$results <- lapply(records, function(x) list(seq = x$sequence))
snapshot$summary <- data.frame()
snapshot$upload_files <- data.frame()
snapshot$settings <- list()
snapshot$assay_profiles <- data.frame()
snapshot$project_defaults <- list()
snapshot$rename <- data.frame()
snapshot$architecture <- list()
snapshot$blast_jobs <- data.frame()
snapshot$blast_hits <- data.frame()
snapshot$blast_raw <- list()
snapshot$taxonomy_summary <- data.frame()
snapshot$taxonomy_hits <- data.frame()
snapshot$taxonomy_counts <- data.frame()
snapshot$project_mode <- "simple"
snapshot$snapshot_id <- "18_read_smoke_snapshot"
selection <- pitax_export_default_selection("qc", "stage")
selection$component_ids <- "qc_trimmed_sequences"
selection$include_plots <- FALSE
selection$sequence_header <- "custom"
selection$header_config <- config
plan <- pitax_export_build_plan(snapshot, selection)
stopifnot(identical(plan$sequence_components, "qc_trimmed_sequences"),
          length(Filter(function(f) identical(f$kind, "fasta"), plan$files)) == 18L,
          identical(plan$availability$qc_trimmed_sequences$status, "available"))

# End-to-end: the same 18 values that previously crashed modifyList must
# survive build -> validate/rebuild -> writer -> ZIP manifest.
snapshot$app_version <- "3.6.4"
snapshot$schema_version <- 6L
snapshot$captured_at <- Sys.time()
snapshot$timezone <- "UTC"
snapshot$source_stage <- "qc"
full_zip <- tempfile(fileext = ".zip")
pitax_export_execute(snapshot, plan, full_zip)
stopifnot(file.exists(full_zip), file.info(full_zip)$size > 0)
# The download handler must serve the exact ZIP bytes prepared on the server.
delivered <- tempfile(fileext = ".zip")
pitax_export_transfer_prepared_zip(full_zip, delivered)
stopifnot(file.exists(delivered),
          file.info(delivered)$size == file.info(full_zip)$size,
          identical(unname(tools::md5sum(delivered)), unname(tools::md5sum(full_zip))))
unpacked <- tempfile("pitax_18_record_export_")
dir.create(unpacked)
utils::unzip(full_zip, exdir = unpacked)
fasta_files <- list.files(unpacked, pattern = "[.]fasta$", recursive = TRUE, full.names = TRUE)
stopifnot(length(fasta_files) == 18L)
first_fasta <- file.path(unpacked, "sequences", "qc_trimmed_sequences",
                         "ISO1_TEF_R.fasta")
if (!file.exists(first_fasta)) {
  stop("ZIP archive layout wrong: ", paste(utils::unzip(full_zip, list = TRUE)$Name,
                                          collapse = ", "))
}
first_header <- readLines(first_fasta, n = 1L)
stopifnot(identical(first_header, ">ISO1_TEF_R organism=En"),
          !grepl("experiment=", first_header, fixed = TRUE))
manifest <- utils::read.csv(file.path(unpacked, "sequences", "FASTA_header_manifest.csv"),
                            stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 18L,
          !any(grepl("experiment=", manifest$Header, fixed = TRUE)))

cat("End-to-end 18-read FASTA ZIP, plan isolation, optional template fields and empty-table smoke passed (18 records).\n")
