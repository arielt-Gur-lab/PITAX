# PITAX v3.0.0-alpha.4 - user-facing workflow order smoke tests.

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

must_contain(app_text, paste(
  '"Upload",',
  '    "Assay & Trim",',
  '    "Rename",',
  '    "QC",',
  sep = "\n"
))
must_contain(app_text, 'tabPanel("3 · Rename", value = "rename"')
must_contain(app_text, 'tabPanel("4 · QC & Chromatogram", value = "qc"')
must_contain(app_text, 'actionButton("to_qc", "Continue to QC"')
must_contain(app_text, 'actionButton("to_export", "Continue to Export"')
must_contain(app_text, 'Loading Rename workspace')
must_contain(app_text, 'Loading renamed QC workspace')
must_contain(app_text, 'renameTab.insertBefore(qcTab)')
must_contain(app_text, '_checkpoint_A_rename.zip')
must_contain(app_text, '_checkpoint_B_qc.zip')
must_contain(app_text, 'df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))')
must_contain(sequence_text, 'result$display_name')
must_contain(sequence_text, 'yaxis = list(rangemode = "fixed", range = c(0, ymax))')

if (grepl('actionButton("to_rename"', app_text, fixed = TRUE)) {
  stop("Legacy QC-to-Rename navigation button is still present.")
}

cat("v3.0.0-alpha.4 workflow order tests passed.\n")
