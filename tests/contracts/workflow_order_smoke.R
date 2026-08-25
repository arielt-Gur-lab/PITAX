# PITAX workflow order and gated Steps UI regression tests.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
source(file.path(app_dir, "tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
sequence_text <- pitax_read_text("R", "domain", "sanger", "sequence_tools.R")

must_contain <- function(text, value) {
  if (!grepl(value, text, fixed = TRUE)) stop("Expected workflow marker was not found: ", value)
}

must_not_contain <- function(text, value) {
  if (grepl(value, text, fixed = TRUE)) stop("Unexpected legacy workflow marker remains: ", value)
}

must_contain(app_text, 'pitax_workflow_stepper_ui')
must_contain(app_text, 'workflow_mark_unlocked')
must_contain(app_text, 'workflow-nav-grid')
must_contain(app_text, 'workflow-col-label')
must_contain(app_text, 'workflow_open_help')
must_contain(app_text, 'workflow_stage_actions')
must_contain(app_text, 'free_nav = TRUE')
app_ui_text <- pitax_read_text("R", "ui", "app_ui.R")
if (grepl("stage_topbar\\(", app_ui_text)) stop("stage_topbar() must not remain inside app_ui.R tab panels")
must_contain(app_text, 'tabsetPanel(id = "pipeline_step", type = "hidden"')
must_contain(app_text, 'uiOutput("workflow_stepper")')
must_not_contain(app_text, 'tabsetPanel(id = "pipeline_step", type = "tabs"')
must_contain(app_text, 'tabPanel("Assay", value = "settings"')
must_contain(app_text, 'tabPanel("Assign", value = "rename"')
must_contain(app_text, 'tabPanel("Trim & QC", value = "qc"')
# Assign must appear before Trim & QC in the tabset source order.
assign_pos <- regexpr('tabPanel("Assign", value = "rename"', app_text, fixed = TRUE)[1]
qc_pos <- regexpr('tabPanel("Trim & QC", value = "qc"', app_text, fixed = TRUE)[1]
if (assign_pos < 1 || qc_pos < 1 || assign_pos > qc_pos) {
  stop("Assign tabPanel must be declared before Trim & QC in app_ui.R")
}
must_contain(app_text, 'uiOutput("assignment_editor")')
must_not_contain(app_text, 'DTOutput("assignment_upload_table")')
must_contain(app_text, 'actionButton("to_rename", "Continue to Assign"')
must_contain(app_text, 'actionButton("run_trimming", "Start trimming"')
must_contain(app_text, 'updateTabsetPanel(session,"pipeline_step",selected="qc")')
must_contain(app_text, 'actionButton("to_blast", "Continue to NCBI BLAST"')
must_contain(app_text, 'actionButton("open_export_output"')
must_contain(app_text, 'assay_add_profile')
must_not_contain(app_text, 'renameTab.insertBefore(qcTab)')
must_contain(app_text, '_checkpoint_A_rename.zip')
must_contain(app_text, '_checkpoint_B_qc.zip')
must_contain(app_text, 'df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))')
must_contain(sequence_text, 'result$display_name')
must_contain(sequence_text, 'yaxis = list(rangemode = "fixed", range = c(0, ymax))')

must_not_contain(app_text, 'actionButton("to_qc"')
must_not_contain(app_text, 'actionButton("run_trimming", "Run trimming"')
must_not_contain(app_text, 'Loading Rename workspace')
must_not_contain(app_text, 'Continue to Export')

cat("Workflow order and gated Steps UI tests passed.\n")
