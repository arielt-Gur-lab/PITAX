# PITAX v3.0.0-alpha.9.2 - retained Stage 3 Shiny integration and startup tests.

source(file.path("tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
consensus_text <- pitax_read_text("R", "domain", "consensus", "stage3_consensus.R")
export_text <- pitax_read_text("R", "export", "export_tools.R")

must_contain <- function(text, marker) {
  if (!grepl(marker, text, fixed = TRUE)) stop(paste("Missing Stage 3 integration marker:", marker), call. = FALSE)
}

must_contain(app_text, 'pitax_source(file.path("R", "domain", "consensus", "stage3_consensus.R"), local = TRUE)')
must_contain(app_text, "PROJECT_SCHEMA_VERSION <- 6L")
must_contain(app_text, 'tabPanel("5 \\u00B7 Analysis Sequence", value = "consensus"')
must_contain(app_text, 'actionButton("to_consensus", "Continue to Consensus"')
must_contain(app_text, 'stage3_build_consensus_set(')
must_contain(app_text, 'stage3_consensus_gate_error(rv$consensus_set, rv$results)')
must_contain(app_text, 'rv$consensus_set <- stage3_empty_consensus_set()')
must_contain(app_text, 'if (source_schema == 5L) st <- assay_migrate_schema5_state(st)')
must_contain(app_text, 'consensus_set = rv$consensus_set')
must_contain(app_text, 'stage3_analysis_records(rv$consensus_set)')
must_contain(app_text, 'write_consensus_checkpoint_zip(file, rv$consensus_set, rv$read_assignments, rv$settings)')
must_contain(app_text, 'locus_error <- stage3_run_locus_error(rv$read_assignments)')
must_contain(app_text, 'project_mode = rv$project_mode')
must_contain(app_text, 'project_mode = rv$project_mode,')
must_contain(app_text, 'output$consensus_mode_note <- renderUI({')
must_contain(app_text, 'hideTab(inputId = "pipeline_step", target = "consensus"')
must_contain(app_text, 'sync_project_mode_navigation <- function(mode)')
must_contain(app_text, 'session$onFlushed(function() sync_project_mode_navigation(isolate(rv$project_mode)), once = TRUE)')
must_contain(app_text, 'if (identical(rv$project_mode, "simple")) {')
must_contain(app_text, 'build_analysis_sequences(notify = TRUE)')
must_contain(app_text, 'actionButton("apply_consensus_review", "Apply decision"')
must_contain(app_text, 'output$consensus_forward_conflict_plot <- plotly::renderPlotly')
must_contain(app_text, 'output$consensus_reverse_conflict_plot <- plotly::renderPlotly')
must_contain(app_text, 'plotly::plotlyOutput("consensus_forward_conflict_plot"')
must_contain(app_text, 'plotly::plotlyOutput("consensus_reverse_conflict_plot"')
if (grepl("(?<!:)plotlyOutput\\(", app_text, perl = TRUE)) stop("An unqualified plotlyOutput() call can prevent application startup.", call. = FALSE)
if (grepl("sync_project_mode_navigation\\(\\)", app_text, perl = TRUE)) stop("Project-mode navigation was called without an explicit or isolated mode.", call. = FALSE)
nav_start <- regexpr("sync_project_mode_navigation <- function(mode)", app_text, fixed = TRUE)[1]
nav_end <- regexpr("session$onFlushed", app_text, fixed = TRUE)[1]
if (nav_start < 1L || nav_end <= nav_start) stop("Could not inspect project-mode navigation helper.", call. = FALSE)
nav_helper <- substr(app_text, nav_start, nav_end - 1L)
if (grepl("rv$", nav_helper, fixed = TRUE) || grepl("input$", nav_helper, fixed = TRUE)) stop("Navigation helper must remain reactive-free and receive mode explicitly.", call. = FALSE)
must_contain(app_text, 'stage3_consensus_undo(old)')
must_contain(app_text, 'stage3_consensus_redo(old)')
must_contain(app_text, 'autoWidth = TRUE')
must_contain(app_text, 'consensus_revision=if (is.list(r$consensus$curation))')

must_contain(consensus_text, 'algorithm = "pitax-overlap-consensus-v1"')
must_contain(consensus_text, 'stage3_reverse_complement')
must_contain(consensus_text, 'stage3_overlap_align')
must_contain(consensus_text, '"Conflict \\u00B7 IUPAC review call"')
must_contain(consensus_text, 'status <- "NO_RELIABLE_OVERLAP"')
must_contain(consensus_text, '"SINGLE_READ"')
must_contain(consensus_text, '"INDEPENDENT_READ"')
must_contain(consensus_text, 'source_sequences')
must_contain(consensus_text, 'Forward_Raw_Position')
must_contain(consensus_text, 'Reverse_Raw_Position')
must_contain(consensus_text, 'stage3_apply_consensus_call')
must_contain(consensus_text, 'stage3_empty_consensus_audit')
must_contain(consensus_text, 'Review_Revision')

must_contain(export_text, 'write_consensus_artifacts')
must_contain(export_text, '"consensus_summary.csv"')
must_contain(export_text, '"isolate_level_sequences.fasta"')
must_contain(export_text, '"_column_evidence.csv"')

# Source the complete application and create a mock Shiny session. This catches
# UI-construction failures and illegal reactive-value reads in one-time
# callbacks, which text-presence contracts alone cannot detect.
app_env <- new.env(parent = globalenv())
sys.source("app.R", envir = app_env)
stopifnot(is.function(app_env$server))
fixture_path <- normalizePath(file.path("tests", "fixtures", "Stage3_synthetic_pair", "TESTPAIR001_forward.ab1"), winslash = "/", mustWork = TRUE)
fixture_upload <- data.frame(
  name = basename(fixture_path), size = file.info(fixture_path)$size,
  type = "application/octet-stream", datapath = fixture_path,
  stringsAsFactors = FALSE
)
shiny::testServer(app_env$server, {
  session$flushReact()
  session$setInputs(project_mode = "paired_consensus")
  session$flushReact()
  session$setInputs(project_mode = "simple")
  session$flushReact()
  session$setInputs(
    target = "ITS", sequencing_primer = "Forward",
    forward_primer = "", reverse_primer = "",
    ab1_files = fixture_upload
  )
  session$flushReact()
})

cat("v3.0.0-alpha.9.2 retained Stage 3 app integration and startup tests passed.\n")
