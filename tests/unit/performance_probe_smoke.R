# Opt-in session-only performance metrics and browser probe wiring.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "services", "perf_diagnostics.R"))

store <- pitax_perf_new(max_events = 10L)
stopifnot(isFALSE(store$enabled))
stopifnot(nrow(store$events) == 0L)
pitax_perf_append(store, "Server", "hidden", 30)
stopifnot(nrow(store$events) == 0L)
stopifnot(identical(pitax_perf_measure(store, "disabled", 123L), 123L))
stopifnot(nrow(store$events) == 0L)

store$enabled <- TRUE
stopifnot(identical(pitax_perf_measure(store, "noop", { 7L * 6L }), 42L))
stopifnot(nrow(store$events) == 1L)
stopifnot(store$events$Source[1] == "Server")
stopifnot(store$events$Operation[1] == "noop")
stopifnot(is.finite(store$events$Duration_ms[1]))
stopifnot(pitax_perf_append(store, "Browser", "guided.row.to_shiny_idle", 83))
stopifnot(nrow(store$events) == 2L)
stopifnot(!pitax_perf_append(store, "Invalid", "test", 2))
stopifnot(!pitax_perf_append(store, "Server", "test", -2))

# Old entries expire; one session cannot grow without bound.
for (i in seq_len(50L)) {
  pitax_perf_append(store, "Server", paste0("event_", i), i)
}
stopifnot(nrow(store$events) == 10L)
stopifnot(identical(tail(store$events$Operation, 1L), "event_50"))
stopifnot(identical(store$events$Operation[1L], "event_41"))
summary <- pitax_perf_summary(store$events)
stopifnot(nrow(summary) == 10L)
stopifnot(all(summary$N == 1L))
stopifnot(all(is.finite(summary$Median_ms)))
store$enabled <- FALSE
stopifnot(!pitax_perf_append(store, "Server", "inactive", 15))
stopifnot(nrow(store$events) == 10L)
store$events <- pitax_perf_empty()
stopifnot(nrow(pitax_perf_summary(store$events)) == 0L)

# The UI is opt-in and browser probes are silent until explicitly checked.
app_ui <- paste(readLines(file.path(root, "R", "ui", "app_ui.R"), warn = FALSE), collapse = "\n")
js <- paste(readLines(file.path(root, "www", "pitax.js"), warn = FALSE), collapse = "\n")
server <- paste(readLines(file.path(root, "R", "server", "stages", "00_perf.R"), warn = FALSE), collapse = "\n")
qc <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"), warn = FALSE), collapse = "\n")
project <- paste(readLines(file.path(root, "R", "server", "stages", "10_project.R"), warn = FALSE), collapse = "\n")
consensus <- paste(readLines(file.path(root, "R", "server", "stages", "70_consensus.R"), warn = FALSE), collapse = "\n")
stopifnot(all(vapply(c(
  'checkboxInput("pitax_perf_enabled", "Enable diagnostics", value = FALSE)',
  'tags$details(class = "panel-box pitax-dev-performance-panel"',
  'tableOutput("pitax_perf_summary")',
  'tableOutput("pitax_perf_recent")',
  'downloadButton("pitax_perf_download"'
), grepl, logical(1), x = app_ui, fixed = TRUE)))
stopifnot(all(vapply(c(
  "perf_store <- pitax_perf_new", "pitax_perf_browser_event",
  "output$pitax_perf_download", "pitax_perf_summary"
), grepl, logical(1), x = server, fixed = TRUE)))
stopifnot(all(vapply(c(
  "perf_measure", "qc.guided_plot_render", "guided.open_server",
  "qc.peak_flags", "qc.full_chromatogram_build",
  "qc.ambiguous_peak_table_render", "qc.guided_side_ui_render",
  "qc.ambiguous_peak_summary_render"
), grepl, logical(1), x = qc, fixed = TRUE)))
stopifnot(grepl("qc.commit_curated_result", project, fixed = TRUE))
stopifnot(grepl("analysis.rebuild_after_edit", project, fixed = TRUE))
stopifnot(grepl("analysis.single_group_rebuild", consensus, fixed = TRUE))
stopifnot(grepl("analysis.full_project_rebuild", consensus, fixed = TRUE))
stopifnot(grepl("session$onFlushed", server, fixed = TRUE))
stopifnot(grepl("qc.perf_summary_render", server, fixed = TRUE))
stopifnot(grepl("publish = FALSE", server, fixed = TRUE))
stopifnot(grepl("to_first_message", server, fixed = TRUE))
stopifnot(grepl("performance.now()", js, fixed = TRUE))
stopifnot(grepl("'shiny:idle.pitaxDevPerf'", js, fixed = TRUE))
stopifnot(grepl("'shiny:message.pitaxDevPerf'", js, fixed = TRUE))
stopifnot(grepl("'shiny:value.pitaxDevPerf'", js, fixed = TRUE))
stopifnot(grepl("__pitaxPerfLog", js, fixed = TRUE))
stopifnot(grepl("__PITAX_PERF_CONSOLE__", js, fixed = TRUE))

# Do not silently persist diagnostics into .sangerproject payload.
a <- regexpr("make_project_bundle <- function", project, fixed = TRUE)[1]
b <- regexpr("restore_settings_inputs <- function", project, fixed = TRUE)[1]
stopifnot(a > 0L, b > a)
bundle <- substr(project, a, b - 1L)
stopifnot(!grepl("perf_store =", bundle, fixed = TRUE))
stopifnot(!grepl("perf_events =", bundle, fixed = TRUE))

cat("performance_probe_smoke: OK\n")
