# PITAX v3.0.0-alpha.6 - user-facing workflow order regression tests.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)
app_text <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
sequence_text <- paste(readLines(file.path(app_dir, "R", "sequence_tools.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

must_contain <- function(text, value) {
  if (!grepl(value, text, fixed = TRUE)) stop("Expected workflow marker was not found: ", value)
}

must_not_contain <- function(text, value) {
  if (grepl(value, text, fixed = TRUE)) stop("Unexpected legacy workflow marker remains: ", value)
}

must_contain(app_text, paste(
  '"Upload",',
  '    "Assay settings",',
  '    "Rename & assign",',
  '    "Trim & QC",',
  sep = "\n"
))
must_contain(app_text, 'tabPanel("2 · Assay", value = "settings"')
must_contain(app_text, 'tabPanel("3 · Rename & Assign", value = "rename"')
must_contain(app_text, 'tabPanel("4 · Trim & QC", value = "qc"')
must_contain(app_text, 'DTOutput("assignment_upload_table")')
must_contain(app_text, 'actionButton("to_rename", "Continue to Rename"')
must_contain(app_text, 'actionButton("run_trimming", "Start trimming"')
must_contain(app_text, 'updateTabsetPanel(session,"pipeline_step",selected="qc")')
must_contain(app_text, 'actionButton("to_export", "Continue to Export"')
must_contain(app_text, 'renameTab.insertBefore(qcTab)')
must_contain(app_text, '_checkpoint_A_rename.zip')
must_contain(app_text, '_checkpoint_B_qc.zip')
must_contain(app_text, 'df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))')
must_contain(sequence_text, 'result$display_name')
must_contain(sequence_text, 'yaxis = list(rangemode = "fixed", range = c(0, ymax))')

must_not_contain(app_text, 'actionButton("to_qc"')
must_not_contain(app_text, 'actionButton("run_trimming", "Run trimming"')
must_not_contain(app_text, 'Loading Rename workspace')

cat("v3.0.0-alpha.6 workflow order tests passed.\n")
