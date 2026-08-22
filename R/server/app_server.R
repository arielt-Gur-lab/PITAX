# ============================================================
# Server composition
# ============================================================

PITAX_SERVER_MODULES <- file.path(
  "R", "server", "stages",
  c(
    "00_state.R",
    "10_project.R",
    "20_upload.R",
    "30_trimming.R",
    "40_qc_summary.R",
    "50_evidence_review.R",
    "60_assignment.R",
    "70_consensus.R",
    "80_export.R",
    "90_blast.R",
    "100_taxonomy.R",
    "110_multilocus.R",
    "120_reset.R"
  )
)

server <- function(input, output, session) {
  server_environment <- environment()
  for (module_path in PITAX_SERVER_MODULES) {
    source(module_path, local = server_environment, encoding = "UTF-8")
  }
}
