# Base-R export smoke regression: no CRAN packages required.
# Intended for GitHub Actions startup/CI preflight and local development.
stopifnot(file.exists("VERSION.txt"))
for (path in c("R/export/export_tools.R", "R/export/availability.R",
               "R/export/plan.R", "R/export/writers.R",
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
config$draft_template <- paste0(
  ">{Final_Name} organism={Organism} isolate={Isolate} locus={Locus}",
  ' direction={Direction} molecule=DNA experiment="{Experiment}" date={Seq_Date}'
)
records <- list()
for (i in seq_len(18L)) {
  rid <- paste0("S", i)
  records[[rid]] <- list(final_name = paste0("ISO", i, "_TEF_R"),
                         sequence = "ACTGCATG", direction = "Reverse",
                         isolate = paste0("ISO", i), locus = "TEF", missing = FALSE)
  key <- pitax_header_record_key("qc_trimmed_sequences", rid)
  for (pair in list(c("Organism", "En"), c("Experiment", "Resistance_X"),
                    c("Seq_Date", "2026-11-05"))) {
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
            ' direction=Reverse molecule=DNA experiment="Resistance_X" date=2026-11-05')))

# Removed tokens must not require stored values.
config$draft_template <- ">{Final_Name} organism={Organism}"
config$values$Experiment[] <- ""
config$values$Seq_Date[] <- ""
rendered <- pitax_header_apply_to_records(records, "qc_trimmed_sequences", config)
stopifnot(identical(rendered[["S18"]]$custom_header, "ISO18_TEF_R organism=En"))
# Active template controls the visible table fields; old saved templates do not
# force Experiment into newly saved profiles that omit it.
config$draft_template <- ">{Final_Name} organism={Organism}"
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
          !("Experiment" %in% names(display)),
          !("Seq_Date" %in% names(display)),
          "Organism" %in% names(display),
          identical(pitax_header_active_metadata_fields(config), "Organism"))
config$templates[["Without Experiment"]] <- config$draft_template
config$active_template <- "Without Experiment"
restored <- pitax_header_coerce_config(unserialize(serialize(config, NULL)))
stopifnot(identical(pitax_header_active_metadata_fields(restored), "Organism"))

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

cat("Base-R FASTA plan isolation, optional template fields and empty-table smoke passed (18 records).\n")
