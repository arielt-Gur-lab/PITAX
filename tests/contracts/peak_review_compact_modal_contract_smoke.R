# Compact guided-review contract: preserve the quick-action workflow above
# the fold and make verbose diagnostics optional rather than required scroll.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
server <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"), warn = FALSE), collapse = "\n")
css <- paste(readLines(file.path(root, "www", "pitax.css"), warn = FALSE), collapse = "\n")

required <- c(
  "plotOutput(\"peak_review_local_plot\", height = \"340px\")",
  'class = "peak-review-metrics peak-review-metrics-main"',
  'class = "peak-review-evidence-details"',
  'tags$summary("Full signal evidence")',
  'class = "peak-review-dye-legend"',
  'class = "peak-review-yscale-toggle',
  '"peak_review_yscale_focus"',
  '"peak_review_yscale_full"',
  'class = "peak-review-trim-inline"',
  'uiOutput("peak_review_advanced_ui")',
  'class = "peak-review-primary-actions"',
  '"peak_review_continue"',
  '"peak_review_change"',
  '"peak_review_flag"',
  '"manual_base_choice"',
  '"manual_iupac_choice"',
  '"peak_review_apply_iupac"',
  '"peak_review_close"'
)
for (token in required) stopifnot(grepl(token, server, fixed = TRUE))

# The details panel must start collapsed and retain full QC diagnostics when
# expanded; no evidence/curation decision is silently discarded.
a <- regexpr("peak_review_metrics_ui <- function", server, fixed = TRUE)[1]
b <- regexpr("show_peak_review_shell <- function", server, fixed = TRUE)[1]
stopifnot(a > 0L, b > a)
metrics <- substr(server, a, b - 1L)
stopifnot(!grepl('open = TRUE', metrics, fixed = TRUE))
stopifnot(all(vapply(c(
  "Original basecaller Q", "Peak prominence", "Peak separation",
  "Called signal", "Competitor signal", "Peak offset", "Flag", "Status"
), function(x) grepl(x, metrics, fixed = TRUE), logical(1))))

# The four dye labels must sit in HTML ABOVE the image, not overlay calls.
# Manual trim actions now stay in the visible left column instead of hidden.
stopifnot(grepl('.peak-review-dye-legend', css, fixed = TRUE))
stopifnot(grepl('.peak-review-yscale-toggle', css, fixed = TRUE))
stopifnot(grepl('.peak-review-trim-controls', css, fixed = TRUE))
stopifnot(!grepl('tags$summary("Advanced trim options")', server, fixed = TRUE))
stopifnot(grepl('"context_trim_left"', server, fixed = TRUE))
stopifnot(grepl('"context_trim_right"', server, fixed = TRUE))
seq_tools <- paste(readLines(file.path(root, "R", "domain", "sanger", "sequence_tools.R"),
                             warn = FALSE), collapse = "\n")
stopifnot(grepl("pitax_chromatogram_focus_yscale <- function", seq_tools, fixed = TRUE))
start <- regexpr("draw_chromatogram_focus_fast <- function", seq_tools, fixed = TRUE)[1]
end <- regexpr("make_chromatogram_focus_plot <- function", seq_tools, fixed = TRUE)[1]
stopifnot(start > 0L, end > start)
local_plot <- substr(seq_tools, start, end - 1L)
stopifnot(!grepl('graphics::legend(', local_plot, fixed = TRUE))
stopifnot(grepl("yscale = c(\"focus\", \"full\")", local_plot, fixed = TRUE) ||
            grepl('yscale = c("focus", "full")', local_plot, fixed = TRUE))
stopifnot(grepl("pitax_chromatogram_focus_yscale", local_plot, fixed = TRUE))
stopifnot(grepl("Clipped peak", local_plot, fixed = TRUE))

# The old tall fixed layout caused routine scrolling. The CSS caps the modal
# height for tiny viewports, but desktop opens with 2 columns and compact
# key evidence; advanced material requires an explicit expansion.
stopifnot(grepl('#shiny-modal:has(.pitax-peak-review-shell-marker) .modal-content', css, fixed = TRUE))
stopifnot(grepl('max-height:calc(100dvh - 24px)', css, fixed = TRUE))
stopifnot(grepl('.peak-review-metrics-main', css, fixed = TRUE))
stopifnot(grepl('grid-template-columns:repeat(2,minmax(0,1fr))', css, fixed = TRUE))
stopifnot(grepl('max-height:140px', css, fixed = TRUE))
stopifnot(grepl('overflow-y:auto', css, fixed = TRUE))
stopifnot(grepl('@media (max-width:780px)', css, fixed = TRUE))

cat("peak_review_compact_modal_contract_smoke: OK\n")
