# Shared paths and source-corpus helpers for the PITAX test suite.

pitax_project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

pitax_source <- function(...) {
  source(file.path(pitax_project_root, ...))
}

pitax_read_text <- function(...) {
  path <- file.path(pitax_project_root, ...)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

pitax_app_contract_files <- c(
  "app.R",
  file.path("R", "app", "bootstrap.R"),
  file.path("R", "ui", "components.R"),
  file.path("R", "ui", "app_ui.R"),
  file.path("R", "server", "app_server.R"),
  file.path("R", "server", "stages", c(
    "00_state.R", "10_project.R", "20_upload.R", "30_trimming.R",
    "40_qc_summary.R", "50_evidence_review.R", "60_assignment.R",
    "70_consensus.R", "80_export.R", "90_blast.R", "100_taxonomy.R",
    "110_multilocus.R", "120_reset.R"
  )),
  file.path("www", "pitax.js"),
  file.path("www", "pitax.css")
)

pitax_read_app_contract <- function() {
  paste(vapply(pitax_app_contract_files, function(path) {
    pitax_read_text(path)
  }, character(1)), collapse = "\n")
}
