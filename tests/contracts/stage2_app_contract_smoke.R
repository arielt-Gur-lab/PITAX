# PITAX v3.0.0-alpha.7.1 - Stage 2 closed-gate regression contracts.

source(file.path("tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
architecture_text <- pitax_read_text("R", "domain", "assignment", "stage2_architecture.R")
export_text <- pitax_read_text("R", "export", "export_tools.R")

must_contain <- function(text, marker) {
  if (!grepl(marker, text, fixed = TRUE)) stop(paste("Missing Stage 2 integration marker:", marker), call. = FALSE)
}

must_not_contain <- function(text, marker) {
  if (grepl(marker, text, fixed = TRUE)) stop(paste("Unexpected marker in corrected Stage 2 alpha:", marker), call. = FALSE)
}

must_contain(app_text, 'pitax_source(file.path("R", "domain", "assignment", "stage2_architecture.R"), local = TRUE)')
must_contain(app_text, "PROJECT_SCHEMA_VERSION <- 6L")
must_contain(app_text, 'actionButton("to_rename", "Continue to Rename"')
must_contain(app_text, 'actionButton("run_trimming", "Start trimming"')
must_contain(app_text, 'uiOutput("assignment_editor")')
must_contain(app_text, 'actionButton("save_assignment_edits"')
must_contain(app_text, 'actionButton("apply_assignment_batch"')
must_contain(app_text, 'tags$option(value = "Forward"')
must_contain(app_text, 'tags$option(value = "Reverse"')
must_contain(app_text, 'collect_assignment_editor <- function()')
must_contain(app_text, 'radioButtons(')
must_contain(app_text, '"Simple reads \\u2014 no Forward/Reverse matching" = "simple"')
must_contain(app_text, '"Paired reads \\u2014 build Forward/Reverse consensus" = "paired_consensus"')
must_contain(app_text, 'Assignment key must contain old_id, isolate, direction, and assay_id or locus (gene).')
must_contain(app_text, 'name_error <- stage2_identity_error(rv$read_assignments)')
must_contain(app_text, 'sync_assignment_state()')
must_contain(app_text, 'stage2_sync_generated_names(')
must_contain(app_text, 'if (!length(rv$results) || is.null(rv$rename)) {')
must_contain(app_text, 'updateSelectInput(session, "blast_sample", choices = character(), selected = character())')
must_contain(app_text, 'if (!length(keys)) {')
must_contain(app_text, 'updateSelectInput(session, "inspect_sample", choices = character(), selected = character())')
must_contain(app_text, 'result$processing_settings <- read_settings')
must_contain(app_text, 'result$read_assignment <- as.list(assignment[1, , drop = FALSE])')
must_contain(app_text, 'if (source_schema == 5L) st <- assay_migrate_schema5_state(st)')
must_contain(app_text, 'build_taxonomic_consensus(enriched, target=sample_target, top_n=top_n)')
must_contain(app_text, 'write_assignment_checkpoint_zip(file, rv$rename, rv$read_assignments, rv$architecture, rv$settings)')
must_contain(export_text, 'write_assignment_checkpoint_zip <- function')
must_contain(architecture_text, 'assignments$Final_Name[i] <- stage2_compose_read_name(assignments$Isolate[i], assignments$Locus[i], assignments$Direction[i])')
must_contain(app_text, 'uiOutput("batch_assay_control")')
must_contain(app_text, 'assignment_input_id("assign_assay", read_id)')

blast_observer_pos <- regexpr('# Keep the BLAST sequence selector synchronized', app_text, fixed = TRUE)[1]
blast_guard_pos <- regexpr('if (!length(rv$results) || is.null(rv$rename)) {', app_text, fixed = TRUE)[1]
blast_set_names_tail <- substr(app_text, blast_observer_pos, nchar(app_text))
blast_set_names_pos <- regexpr('choices <- setNames(names(rec), final_names)', blast_set_names_tail, fixed = TRUE)[1] + blast_observer_pos - 1L
if (blast_observer_pos < 1L || blast_guard_pos < blast_observer_pos || blast_set_names_pos <= blast_guard_pos) {
  stop("The BLAST selector must guard the empty pre-trim state before setNames().", call. = FALSE)
}

qc_sync_pos <- regexpr('sync_qc_sample_choices <- function', app_text, fixed = TRUE)[1]
qc_guard_tail <- substr(app_text, qc_sync_pos, nchar(app_text))
qc_guard_pos <- regexpr('if (!length(keys)) {', qc_guard_tail, fixed = TRUE)[1]
qc_set_names_pos <- regexpr('choices <- stats::setNames(keys, labels)', qc_guard_tail, fixed = TRUE)[1]
if (qc_sync_pos < 1L || qc_guard_pos < 1L || qc_set_names_pos <= qc_guard_pos) {
  stop("QC selector synchronization must guard an all-failed/empty result set before setNames().", call. = FALSE)
}

to_rename_pos <- regexpr('observeEvent(input$to_rename, {', app_text, fixed = TRUE)[1]
run_trimming_pos <- regexpr('observeEvent(input$run_trimming, {', app_text, fixed = TRUE)[1]
if (to_rename_pos < 1L || run_trimming_pos <= to_rename_pos) {
  stop("Rename transition must be registered before the trimming action.", call. = FALSE)
}
to_rename_contract <- substr(app_text, to_rename_pos, run_trimming_pos - 1L)
must_not_contain(to_rename_contract, 'trim_one_ab1(')
must_contain(substr(app_text, run_trimming_pos, nchar(app_text)), 'trim_one_ab1(')

final_names_pos <- regexpr('card_title("Final read / FASTA names and biological identity"', app_text, fixed = TRUE)[1]
architecture_pos <- regexpr('card_title("Stage 2 \\u00B7 Architecture preview"', app_text, fixed = TRUE)[1]
if (final_names_pos < 1L || architecture_pos < 1L || final_names_pos >= architecture_pos) {
  stop("Generated final names must appear before the architecture preview.", call. = FALSE)
}

must_not_contain(architecture_text, "stage2_parse_structured_read_name")
must_not_contain(architecture_text, "stage2_detect_direction_suffix")
must_not_contain(architecture_text, "stage2_known_locus")
must_not_contain(app_text, 'DTOutput("read_assignment_table")')
must_not_contain(app_text, 'DTOutput("assignment_upload_table")')
must_not_contain(app_text, 'DTOutput("assignment_review_table")')
must_not_contain(app_text, 'actionButton("to_qc"')
must_not_contain(app_text, 'Loading Rename workspace')

# Stage 2 never performs a hidden merge. Stage 3 is a separate, explicit layer.
must_not_contain(app_text, "build_forward_reverse_consensus")
must_not_contain(app_text, "auto_merge_forward_reverse")

cat("v3.0.0-alpha.7.1 Stage 2 closed-gate integration tests passed.\n")
