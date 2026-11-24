# On-demand QC signal-metrics plots: cache fingerprint, UI contract,
# export Prepare freshness, and Continue-path decoupling.
get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
root <- normalizePath(file.path(get_this_script_dir(), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "domain", "sanger", "sequence_tools.R"))

# --- Fingerprint: revision / trim / settings change marks outdated ---
mk_result <- function(revision = 1L, trim_start = 10L, trim_end = 100L, sample_id = "S1") {
  list(
    sample_id = sample_id,
    curation = list(revision = as.integer(revision)),
    summary = data.frame(
      trim_start = as.integer(trim_start),
      trim_end = as.integer(trim_end),
      collapse_index = NA_integer_,
      stringsAsFactors = FALSE
    ),
    metrics = data.frame(index = 1:3, called_signal = c(1, 2, 3), peak_ratio = c(2, 3, 4))
  )
}
settings_a <- list(min_peak_ratio = 1.5, absolute_max_base_index = 800L)
settings_b <- list(min_peak_ratio = 2.0, absolute_max_base_index = 800L)

fp1 <- qc_plot_fingerprint(mk_result(1L), settings_a)
fp2 <- qc_plot_fingerprint(mk_result(2L), settings_a)
fp3 <- qc_plot_fingerprint(mk_result(1L, trim_start = 20L), settings_a)
fp4 <- qc_plot_fingerprint(mk_result(1L), settings_b)
stopifnot(identical(fp1, qc_plot_fingerprint(mk_result(1L), settings_a)))
stopifnot(!identical(fp1, fp2))
stopifnot(!identical(fp1, fp3))
stopifnot(!identical(fp1, fp4))

# write_qc_plot_png still draws from live result (export path).
tmp <- tempfile(fileext = ".png")
on.exit(unlink(tmp), add = TRUE)
write_qc_plot_png(mk_result(), settings_a, tmp)
stopifnot(file.exists(tmp), file.info(tmp)$size > 0)

ui <- paste(readLines(file.path(root, "R", "ui", "app_ui.R"), warn = FALSE), collapse = "\n")
qc <- paste(readLines(file.path(root, "R", "server", "stages", "50_evidence_review.R"), warn = FALSE), collapse = "\n")
export <- paste(readLines(file.path(root, "R", "server", "stages", "80_export.R"), warn = FALSE), collapse = "\n")
css <- paste(readLines(file.path(root, "www", "pitax.css"), warn = FALSE), collapse = "\n")

# Generate button + imageOutput; no auto plotOutput("qc_plot").
stopifnot(grepl('"qc_plots_generate"', ui, fixed = TRUE))
stopifnot(grepl("Generate / Refresh QC Plots", ui, fixed = TRUE))
stopifnot(grepl('imageOutput("qc_plot"', ui, fixed = TRUE))
stopifnot(!grepl('plotOutput("qc_plot"', ui, fixed = TRUE))
stopifnot(grepl('uiOutput("qc_plot_status"', ui, fixed = TRUE))

# Server: on-demand Generate measures plot; renderImage must not call
# draw_qc_metrics(selected_result()) reactively.
stopifnot(grepl('observeEvent(input$qc_plots_generate', qc, fixed = TRUE))
stopifnot(grepl('perf_measure("qc.signal_metrics_plot"', qc, fixed = TRUE))
stopifnot(grepl("write_qc_plot_png", qc, fixed = TRUE))
stopifnot(grepl("qc_plot_fingerprint", qc, fixed = TRUE))
stopifnot(grepl("Outdated", qc, fixed = TRUE))
stopifnot(grepl("Not generated", qc, fixed = TRUE))
stopifnot(grepl("Current", qc, fixed = TRUE))
stopifnot(grepl("output$qc_plot <- renderImage", qc, fixed = TRUE))
stopifnot(!grepl("output$qc_plot <- renderPlot", qc, fixed = TRUE))

# renderImage body must not pull selected_result() (would auto-redraw).
img_start <- regexpr("output$qc_plot <- renderImage", qc, fixed = TRUE)[1]
img_end <- regexpr("session$onSessionEnded(function()", qc, fixed = TRUE)[1]
stopifnot(img_start > 0L, img_end > img_start)
img_block <- substr(qc, img_start, img_end - 1L)
stopifnot(!grepl("selected_result()", img_block, fixed = TRUE))
stopifnot(!grepl("draw_qc_metrics", img_block, fixed = TRUE))

# Guided Review local chromatogram stays live.
stopifnot(grepl("output$peak_review_local_plot <- renderPlot", qc, fixed = TRUE))
stopifnot(grepl("draw_chromatogram_focus_fast", qc, fixed = TRUE))
local_start <- regexpr("output$peak_review_local_plot <- renderPlot", qc, fixed = TRUE)[1]
local_end <- regexpr("observeEvent(input$peak_review_yscale_focus", qc, fixed = TRUE)[1]
stopifnot(local_start > 0L, local_end > local_start)
local_block <- substr(qc, local_start, local_end - 1L)
stopifnot(grepl("selected_result()", local_block, fixed = TRUE))

# Collapse contracts: details without default open; no nested Open evidence audit.
stopifnot(grepl('class = "panel-box saved-retrim-card saved-retrim-details"', ui, fixed = TRUE) ||
            grepl("saved-retrim-details", ui, fixed = TRUE))
stopifnot(grepl("saved-retrim-toggle", ui, fixed = TRUE))
stopifnot(grepl("stage1-audit-details", ui, fixed = TRUE))
stopifnot(grepl("Stage 1 | AB1 evidence audit", ui, fixed = TRUE))
stopifnot(!grepl("Open evidence audit", ui, fixed = TRUE))
# Collapsed by default: no open=TRUE on those details blocks.
retrim_a <- regexpr("saved-retrim-details", ui, fixed = TRUE)[1]
retrim_b <- regexpr("Trimming results", ui, fixed = TRUE)[1]
stopifnot(retrim_a > 0L, retrim_b > retrim_a)
retrim_block <- substr(ui, retrim_a - 40L, retrim_b)
stopifnot(!grepl("open\\s*=\\s*TRUE", retrim_block))
stage1_a <- regexpr("stage1-evidence-card", ui, fixed = TRUE)[1]
stage1_b <- regexpr("Ambiguous peak review", ui, fixed = TRUE)[1]
stopifnot(stage1_a > 0L, stage1_b > stage1_a)
stage1_block <- substr(ui, stage1_a - 40L, stage1_b)
stopifnot(!grepl("open\\s*=\\s*TRUE", stage1_block))
stopifnot(!grepl("Open evidence audit", stage1_block, fixed = TRUE))

# Guided Review labels stacked above controls.
stopifnot(grepl("Single base / unresolved", qc, fixed = TRUE))
stopifnot(grepl("flex-direction:column", css, fixed = TRUE))
stopifnot(grepl(".peak-review-modal-side .shiny-input-radiogroup > .control-label", css, fixed = TRUE) ||
            grepl(".peak-review-modal-side .shiny-input-radiogroup > .control-label,", css, fixed = TRUE))
stopifnot(grepl("float:none !important", css, fixed = TRUE))

# Prepare Export refreshes snapshot from live rv.
prep_a <- regexpr("observeEvent(input$export_prepare_download", export, fixed = TRUE)[1]
stopifnot(prep_a > 0L)
prep_block <- substr(export, prep_a, prep_a + 1800L)
stopifnot(grepl("pitax_export_snapshot_from_rv", prep_block, fixed = TRUE))
stopifnot(grepl("export_dialog_state$snap <- snap", prep_block, fixed = TRUE))

cat("qc_plot_ondemand_smoke: OK\n")
