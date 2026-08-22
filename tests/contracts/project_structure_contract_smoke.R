# PITAX - organized project structure contract.

required_files <- c(
  "app.R",
  file.path("R", "app", "bootstrap.R"),
  file.path("R", "ui", "components.R"),
  file.path("R", "ui", "app_ui.R"),
  file.path("R", "server", "app_server.R"),
  file.path("R", "domain", "assay", "assay_profiles.R"),
  file.path("R", "domain", "assignment", "stage2_architecture.R"),
  file.path("R", "domain", "consensus", "stage3_consensus.R"),
  file.path("R", "domain", "multilocus", "stage4_multilocus.R"),
  file.path("R", "domain", "sanger", "ab1_evidence.R"),
  file.path("R", "domain", "sanger", "core_sanger.R"),
  file.path("R", "domain", "sanger", "sequence_tools.R"),
  file.path("R", "services", "taxonomy_tools.R"),
  file.path("R", "export", "export_tools.R"),
  file.path("www", "pitax.css"),
  file.path("www", "pitax.js"),
  file.path("www", "logo.png"),
  file.path("docs", "PITAX_MASTER.md"),
  file.path("templates", "assignment_key_template.csv")
)

missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop("Missing organized project file(s): ", paste(missing, collapse = ", "), call. = FALSE)

runtime_source_files <- c(
  "app.R",
  list.files("R", pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("www", pattern = "[.](js|css)$", recursive = TRUE, full.names = TRUE)
)
contains_non_ascii <- vapply(runtime_source_files, function(path) {
  size <- file.info(path)$size
  if (!is.finite(size) || size < 1L) return(FALSE)
  any(as.integer(readBin(path, what = "raw", n = size)) > 127L)
}, logical(1))
if (any(contains_non_ascii)) {
  stop(
    "Runtime source must remain ASCII-safe: ",
    paste(runtime_source_files[contains_non_ascii], collapse = ", "),
    call. = FALSE
  )
}
runtime_text <- vapply(runtime_source_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
has_unicode_escape <- grepl("\\\\u[0-9A-Fa-f]{4}", runtime_text, perl = TRUE)
if (any(has_unicode_escape)) {
  stop(
    "Runtime source contains Unicode escape sequences that can leak as <U+XXXX> on locale-sensitive hosts: ",
    paste(runtime_source_files[has_unicode_escape], collapse = ", "),
    call. = FALSE
  )
}

app_lines <- readLines("app.R", warn = FALSE, encoding = "UTF-8")
app_text <- paste(app_lines, collapse = "\n")
if (length(app_lines) > 30L) stop("app.R must remain a small composition root.", call. = FALSE)
if (!grepl('pitax_source(file.path("R", "app", "bootstrap.R"), local = TRUE)', app_text, fixed = TRUE)) stop("app.R does not load bootstrap.R through the platform-safe source helper.", call. = FALSE)
if (!grepl('pitax_source(file.path("R", "server", "app_server.R"), local = TRUE)', app_text, fixed = TRUE)) stop("app.R does not load app_server.R through the platform-safe source helper.", call. = FALSE)
if (grepl('encoding = "UTF-8"', app_text, fixed = TRUE)) stop("Runtime modules must not require locale-dependent source conversion.", call. = FALSE)
if (grepl('(^|\n)[[:space:]]*source\\(file\\.path\\(', app_text, perl = TRUE)) stop("Application modules bypass the platform-safe source helper.", call. = FALSE)
if (grepl("reactiveValues(", app_text, fixed = TRUE) || grepl("fluidPage(", app_text, fixed = TRUE)) stop("UI or server implementation leaked back into app.R.", call. = FALSE)

ui_text <- paste(readLines(file.path("R", "ui", "app_ui.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl('href = "pitax.css"', ui_text, fixed = TRUE) || !grepl('src = "pitax.js"', ui_text, fixed = TRUE)) stop("UI does not load the external browser assets.", call. = FALSE)
if (grepl("tags$style(HTML(", ui_text, fixed = TRUE)) stop("Large CSS must remain in www/pitax.css.", call. = FALSE)

server_env <- new.env(parent = baseenv())
sys.source(file.path("R", "server", "app_server.R"), envir = server_env)
server_composition_text <- paste(readLines(file.path("R", "server", "app_server.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl('pitax_source(module_path, local = server_environment)', server_composition_text, fixed = TRUE)) stop("Server stages do not use the platform-safe source helper.", call. = FALSE)
bootstrap_text <- paste(readLines(file.path("R", "app", "bootstrap.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (grepl('encoding = "UTF-8"', bootstrap_text, fixed = TRUE) || grepl('encoding = "UTF-8"', server_composition_text, fixed = TRUE)) stop("Linux/Connect modules must not force locale-dependent UTF-8 conversion.", call. = FALSE)
expected_modules <- file.path(
  "R", "server", "stages",
  c("00_state.R", "10_project.R", "20_upload.R", "30_trimming.R", "40_qc_summary.R",
    "50_evidence_review.R", "60_assignment.R", "70_consensus.R", "80_export.R", "90_blast.R",
    "100_taxonomy.R", "110_multilocus.R", "120_reset.R")
)
if (!identical(server_env$PITAX_SERVER_MODULES, expected_modules)) stop("Server stage order changed unexpectedly.", call. = FALSE)
if (!all(file.exists(expected_modules))) stop("One or more server stage modules are missing.", call. = FALSE)

current_version <- trimws(readLines("VERSION.txt", warn = FALSE, n = 1L))
if (length(current_version) != 1L || !grepl("^[0-9]+[.][0-9]+[.][0-9]+$", current_version)) {
  stop("VERSION.txt must contain one X.Y.Z version string.", call. = FALSE)
}

canonical_logo_md5 <- "f1e67ad4a15d181f7998cf93f81025ed"
actual_logo_md5 <- unname(tools::md5sum(file.path("www", "logo.png")))
if (!identical(actual_logo_md5, canonical_logo_md5)) {
  stop("www/logo.png is not the approved canonical PITAX logo asset.", call. = FALSE)
}

cat("PITAX v", current_version, " project structure contract passed.\n", sep = "")
