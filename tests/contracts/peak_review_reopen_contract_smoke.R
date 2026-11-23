# Regression contract: peak review must reopen when the SAME flagged row is clicked
# after any Close path. DT single selection is not a per-click event.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
stage <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"), warn = FALSE), collapse = "\n")
js <- paste(readLines(file.path(root, "www", "pitax.js"), warn = FALSE), collapse = "\n")
state <- paste(readLines(file.path(root, "R", "server", "stages", "00_state.R"), warn = FALSE), collapse = "\n")

table_start <- regexpr("output\\$ambiguous_peak_table <- renderDT", stage)
table_end <- regexpr("  flag_evidence_text <- function", stage)
stopifnot(table_start[1] > 0L, table_end[1] > table_start[1])
table <- substring(stage, table_start[1], table_end[1] - 1L)
# A row click must send Position as a fresh event every time. Do not depend
# on rows_selected, since clicking an already selected row may do nothing.
stopifnot(grepl('selection = "none"', table, fixed = TRUE))
stopifnot(grepl("table.on('click', 'tbody tr'", table, fixed = TRUE))
stopifnot(grepl('table.row(this).data()', table, fixed = TRUE))
stopifnot(grepl('peak_flag_row_click', table, fixed = TRUE))
stopifnot(grepl("priority: 'event'", table, fixed = TRUE))
stopifnot(!grepl("ambiguous_peak_table_rows_selected", table, fixed = TRUE))

# Modal reopen is forced by row/context clicks, while Continue uses the
# existing shell (not a flickering series of newly mounted dialogs).
stopifnot(grepl("observeEvent(input$peak_flag_row_click", stage, fixed = TRUE))
stopifnot(grepl("open_peak_review_at(pos, reopen_shell = TRUE)", stage, fixed = TRUE))
stopifnot(grepl("show_peak_curation_menu(row[1, , drop = FALSE], reopen_shell = reopen_shell)",
                stage, fixed = TRUE))
stopifnot(grepl("if (isTRUE(reopen_shell) || !isTRUE(rv$peak_review_shell_open))", stage, fixed = TRUE))
stopifnot(grepl('easyClose = TRUE', stage, fixed = TRUE))

# ESC/backdrop/X closure (which does not trigger the Close actionButton) is
# reported with a generation token, so delayed close cannot reset a new modal.
stopifnot(grepl('pitax-peak-review-shell-marker', stage, fixed = TRUE))
stopifnot(grepl('peak_review_shell_generation', state, fixed = TRUE))
stopifnot(grepl("hide.bs.modal.pitaxPeakReview", js, fixed = TRUE))
stopifnot(grepl("hidden.bs.modal.pitaxPeakReview", js, fixed = TRUE))
stopifnot(grepl("peak_review_modal_dismissed", js, fixed = TRUE))
stopifnot(grepl("observeEvent(input$peak_review_modal_dismissed", stage, fixed = TRUE))

cat("peak_review_reopen_contract_smoke: OK\n")
