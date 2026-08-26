# PITAX 3.1.1 - Trim/QC empty-state and amplicon_overview regression smokes.

source(file.path("R", "domain", "sanger", "sequence_tools.R"), local = TRUE)
source(file.path("tests", "helpers", "test_paths.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

# --- draw_amplicon_overview must not throw on empty / missing data ---
pdf(NULL)
ok_empty <- tryCatch({
  draw_amplicon_overview_empty("Trim & QC has nothing to show, make sure there is at least one AB1 file uploaded.")
  TRUE
}, error = function(e) FALSE)
assert_true(ok_empty, "Empty amplicon overview helper must not throw.")

ok_null <- tryCatch({
  draw_amplicon_overview(NULL, list())
  TRUE
}, error = function(e) FALSE)
assert_true(ok_null, "draw_amplicon_overview(NULL) must not throw.")

ok_blank <- tryCatch({
  draw_amplicon_overview(list(sample_id = "x", raw_seq = "", summary = data.frame()), list())
  TRUE
}, error = function(e) FALSE)
assert_true(ok_blank, "draw_amplicon_overview with empty raw_seq must not throw.")

# Valid minimal QC result still renders.
ok_valid <- tryCatch({
  draw_amplicon_overview(
    list(
      sample_id = "S1",
      raw_seq = paste(rep("ACGT", 25), collapse = ""),
      summary = data.frame(trim_start = 10L, trim_end = 80L, stringsAsFactors = FALSE)
    ),
    list(enable_primer_mapping = FALSE)
  )
  TRUE
}, error = function(e) FALSE)
assert_true(ok_valid, "Valid QC result must still render amplicon_overview.")
dev.off()

# --- selected_sample_key / amplicon_overview quiet empty contracts ---
qc_text <- pitax_read_text("R", "server", "stages", "40_qc_summary.R")
ev_text <- pitax_read_text("R", "server", "stages", "50_evidence_review.R")
assert_true(!grepl("validate\\(need\\(FALSE, qc_workspace_empty_message", qc_text),
            "selected_sample_key must not validate() the empty message onto plots.")
assert_true(grepl("req(FALSE)", qc_text, fixed = TRUE),
            "selected_sample_key must use quiet req(FALSE) for missing selection.")
assert_true(grepl("draw_amplicon_overview_empty", ev_text, fixed = TRUE),
            "amplicon_overview must use an explicit empty-state drawer.")
assert_true(grepl("do not fall back to another sample", ev_text, fixed = TRUE),
            "Stale/missing inspect_sample must not silently select another sample in amplicon_overview.")

# sync_qc_sample_choices: stale current falls back only when keys exist; empty clears.
trim_text <- pitax_read_text("R", "server", "stages", "30_trimming.R")
assert_true(grepl("current %in% keys", trim_text, fixed = TRUE),
            "QC sample sync must require the selected key to exist in results.")
assert_true(grepl('updateSelectInput\\(session, "inspect_sample", choices = character\\(\\), selected = character\\(\\)\\)', trim_text),
            "Empty results must clear inspect_sample rather than keep a stale key.")

# --- Plotly empty placeholders set type/mode ---
assert_true(grepl("pitax_empty_plotly", pitax_read_text("R", "domain", "sanger", "sequence_tools.R"), fixed = TRUE),
            "sequence_tools must define pitax_empty_plotly.")
tools_text <- pitax_read_text("R", "domain", "sanger", "sequence_tools.R")
assert_true(grepl('type = "scatter"', tools_text, fixed = TRUE) && grepl('mode = "markers"', tools_text, fixed = TRUE),
            "pitax_empty_plotly must set scatter/markers explicitly.")
assert_true(!grepl("plotly::plot_ly\\(\\)\\s*\\|>\\s*plotly::layout\\(title", tools_text),
            "Chromatogram empty states must not use bare plot_ly() placeholders.")

tax_text <- pitax_read_text("R", "server", "stages", "100_taxonomy.R")
ml_text <- pitax_read_text("R", "server", "stages", "110_multilocus.R")
assert_true(grepl("pitax_empty_plotly", tax_text, fixed = TRUE),
            "Taxonomy empty score plot must use pitax_empty_plotly.")
assert_true(grepl("pitax_empty_plotly", ml_text, fixed = TRUE),
            "Multi-locus empty evidence plot must use pitax_empty_plotly.")

# Runtime: empty plotly helper returns without relying on missing trace type.
if (requireNamespace("plotly", quietly = TRUE)) {
  p <- pitax_empty_plotly("placeholder")
  assert_true(inherits(p, "plotly"), "pitax_empty_plotly must return a plotly object.")
}

cat("Trim/QC empty-state and amplicon_overview smokes passed.\n")
