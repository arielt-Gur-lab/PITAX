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
must_contain(app_text, 'workflow_mark_completed')
must_contain(app_text, 'workspace-chrome')
must_contain(app_text, 'workflow_stage_heading')
must_contain(app_text, 'workflow_stage_actions')
must_contain(app_text, 'workflow_open_help')
must_contain(app_text, 'free_nav = TRUE')
must_contain(app_text, 'workflow_completed_steps')
must_contain(app_text, 'workflow_hide_loader')
must_contain(app_text, 'inputId = paste0("workflow_goto_"')
must_not_contain(app_text, 'step_idx < cur_idx')
must_not_contain(app_text, "Opening Help...")
must_not_contain(app_text, 'workflow_nav_step')
must_not_contain(app_text, "Loading step...")
must_not_contain(app_text, 'has_settings <- is.list(rv$settings)')
# Default assay settings alone must never paint Assay green.
must_not_contain(app_text, 'if (has_settings || has_results) done <- c(done, "settings")')
app_ui_text <- pitax_read_text("R", "ui", "app_ui.R")
if (grepl("stage_heading\\(", app_ui_text)) stop("stage_heading() must not remain inside app_ui.R tab panels")
if (grepl("stage_topbar\\(", app_ui_text)) stop("stage_topbar() must not remain inside app_ui.R tab panels")
must_contain(app_text, 'tabsetPanel(id = "pipeline_step", type = "hidden"')
must_contain(app_text, 'uiOutput("workflow_stepper")')
must_not_contain(app_text, 'tabsetPanel(id = "pipeline_step", type = "tabs"')
# Chip path must be Shiny actionButtons, not raw tags$button + JS setInputValue.
components_text <- pitax_read_text("R", "ui", "components.R")
if (!grepl('inputId = paste0\\("workflow_goto_"', components_text)) {
  stop("workflow chips must use Shiny actionButton(inputId = paste0(\"workflow_goto_\"...))")
}
if (grepl('tags\\$button', components_text) && grepl('workflow_goto_', components_text) &&
    grepl('tags\\$button\\([\\s\\S]{0,200}workflow_goto_', components_text, perl = TRUE)) {
  stop("workflow chips must not use tags$button for workflow_goto_ navigation")
}
js_text <- pitax_read_text("www", "pitax.js")
if (grepl("workflow_nav_step", js_text, fixed = TRUE)) {
  stop("pitax.js must not set workflow_nav_step for chip navigation")
}
if (grepl("Loading step...", js_text, fixed = TRUE)) {
  stop("pitax.js must not show a blocking loader on step chip clicks")
}
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
