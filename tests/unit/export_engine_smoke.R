# PITAX shared export engine behavior tests.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, "..", ".."), mustWork = TRUE)
setwd(app_dir)
source(file.path("tests", "helpers", "test_paths.R"))

sys.source(file.path("R", "domain", "sanger", "ab1_evidence.R"), envir = environment())
sys.source(file.path("R", "domain", "assay", "assay_profiles.R"), envir = environment())
sys.source(file.path("R", "domain", "assignment", "stage2_architecture.R"), envir = environment())
sys.source(file.path("R", "domain", "consensus", "stage3_consensus.R"), envir = environment())
sys.source(file.path("R", "domain", "multilocus", "stage4_multilocus.R"), envir = environment())
sys.source(file.path("R", "domain", "sanger", "core_sanger.R"), envir = environment())
sys.source(file.path("R", "domain", "sanger", "sequence_tools.R"), envir = environment())
sys.source(file.path("R", "export", "export_tools.R"), envir = environment())
sys.source(file.path("R", "export", "catalog.R"), envir = environment())
sys.source(file.path("R", "export", "snapshot.R"), envir = environment())
sys.source(file.path("R", "export", "availability.R"), envir = environment())
sys.source(file.path("R", "export", "header_templates.R"), envir = environment())
sys.source(file.path("R", "export", "plan.R"), envir = environment())
sys.source(file.path("R", "export", "writers.R"), envir = environment())
sys.source(file.path("R", "export", "execute.R"), envir = environment())

fail <- function(...) stop(paste0(...), call. = FALSE)
check <- function(cond, msg) if (!isTRUE(cond)) fail(msg)

make_min_result <- function(sid, seq = "ACGTACGT", settings = list(target = "ITS", window = 11, min_peak_ratio = 1.2)) {
  list(
    sample_id = sid,
    seq = seq,
    raw_seq = seq,
    peak_pos = seq_len(nchar(seq)),
    trace = matrix(1, nrow = nchar(seq), ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T"))),
    channel_map = c(A = 1L, C = 2L, G = 3L, T = 4L),
    metrics = data.frame(
      index = seq_len(nchar(seq)),
      called_signal = rep(100, nchar(seq)),
      peak_ratio = rep(5, nchar(seq)),
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      sample_id = sid, target = "ITS", raw_length = nchar(seq), trimmed_length = nchar(seq),
      trim_start = 1L, trim_end = nchar(seq), collapse_index = NA_integer_, reason = "ok",
      status = "OK", stringsAsFactors = FALSE
    ),
    processing_settings = settings,
    curation = list(
      revision = 1L,
      audit_log = data.frame(
        Timestamp = "2026-01-01 00:00:00", Transaction_ID = "t1", Revision = 1L,
        Action = "Base edit", Position = 2L, Before = "C", After = "G",
        Method = "Manual", Evidence = "test", Details = "active",
        stringsAsFactors = FALSE
      )
    )
  )
}

# --- Assign before trim ------------------------------------------------------
asg <- data.frame(
  Read_ID = c("r1", "r2"),
  Source_ID = c("S1", "S2"),
  File = c("S1_F.ab1", "S2_R.ab1"),
  Final_Name = c("ISO1_ITS_F", "ISO1_ITS_R"),
  Isolate = c("ISO1", "ISO1"),
  Assay_ID = c("assay-its", "assay-its"),
  Locus = c("ITS", "ITS"),
  Direction = c("Forward", "Reverse"),
  Primer = c("ITS1", "ITS4"),
  Inference = c("", ""),
  Notes = c("", ""),
  stringsAsFactors = FALSE
)
state_assign <- list(
  results = list(), summary = data.frame(), rename = data.frame(Original_name = asg$Source_ID, New_name = asg$Final_Name, stringsAsFactors = FALSE),
  settings = list(target = "ITS"), assay_profiles = assay_default_profiles(),
  project_defaults = assay_project_defaults_from_legacy_settings(),
  read_assignments = asg,
  architecture = stage2_build_architecture(asg, assay_profiles = assay_default_profiles()),
  consensus_set = stage3_empty_consensus_set(),
  blast_jobs = data.frame(), blast_hits = data.frame(), blast_ids = data.frame(), blast_raw = list(),
  taxonomy_summary = data.frame(), taxonomy_hits = data.frame(), taxonomy_counts = data.frame(),
  multilocus_profile = stage4_empty_profile()
)
snap_assign <- pitax_export_capture_snapshot(
  state_assign,
  upload_files = data.frame(name = asg$File, size = c(1000, 1100), datapath = c(tempfile(), tempfile()),
                            Source_ID = asg$Source_ID, Bytes_available = FALSE, stringsAsFactors = FALSE),
  pending_curation = list(type = "base", position = 4L, value = "A"),
  source_stage = "rename", project_mode = "paired_consensus", app_version = "3.4.0"
)
sel_assign <- pitax_export_default_selection("rename", "stage")
sel_assign$component_ids <- c("assign_identity_table", "qc_summary")
plan_assign <- pitax_export_build_plan(snap_assign, sel_assign)
check(identical(plan_assign$availability$assign_identity_table$status, "available"), "Assign table should be exportable before trim")
check(identical(plan_assign$availability$qc_summary$status, "not_created"), "QC should be not_created before trim")
check("assign_identity_table" %in% plan_assign$included_components, "Assign component should be included")
check(!"qc_summary" %in% plan_assign$included_components, "QC should be omitted before trim")
check(any(grepl("will NOT be applied", plan_assign$warnings)), "Pending curation must be warned and not applied")

zip_assign <- tempfile(fileext = ".zip")
res_assign <- pitax_export_execute(snap_assign, plan_assign, zip_assign)
check(file.exists(zip_assign), "Assign-before-trim package should be written")
tmpdir <- tempfile("exp_unzip_")
dir.create(tmpdir)
utils::unzip(zip_assign, exdir = tmpdir)
run_info <- paste(readLines(file.path(tmpdir, "RUN_INFO.txt"), warn = FALSE), collapse = "\n")
check(grepl("snapshot_id:", run_info, fixed = TRUE), "RUN_INFO must include snapshot_id")
check(grepl("pending|will NOT be applied|Unconfirmed", run_info), "RUN_INFO must document unapplied edits")
xlsx <- list.files(tmpdir, pattern = "[.]xlsx$", recursive = TRUE, full.names = TRUE)
check(length(xlsx) == 1L, "Assign export should produce one workbook")
wb <- openxlsx::loadWorkbook(xlsx[1])
sheets <- openxlsx::getSheetNames(xlsx[1])
check("Provenance" %in% sheets, "Workbook must keep a provenance sheet with source links")

# --- Row/column selection ----------------------------------------------------
df_id <- pitax_export_extract_assign_identity_table(snap_assign)
sel_cols <- sel_assign
sel_cols$component_ids <- "assign_identity_table"
sel_cols$column_mode <- "custom"
sel_cols$column_ids <- c("Source_ID", "Isolate", "Final_Name")
sel_cols$row_mode <- "selected"
sel_cols$selected_rows <- 1L
comp <- pitax_export_get_component("assign_identity_table")
filtered <- pitax_export_apply_row_column_filters(df_id, comp, sel_cols)
check(nrow(filtered) == 1L, "Selected-row export must keep one row")
check(identical(names(filtered), c("Source_ID", "Isolate", "Final_Name")), "Custom columns must be honored")

# --- Multi-locus sequences and name collisions -------------------------------
seq_recs <- list(
  a = list(id = "a", final_name = "ISO1 ITS", sequence = "AAA", sequence_role = "multilocus_locus",
           direction = "oriented", isolate = "ISO1", locus = "ITS", missing = FALSE),
  b = list(id = "b", final_name = "ISO1 ITS", sequence = "CCC", sequence_role = "multilocus_locus",
           direction = "oriented", isolate = "ISO1", locus = "TEF1", missing = FALSE),
  c = list(id = "c", final_name = "ISO1_LSU", sequence = "", sequence_role = "multilocus_locus",
           direction = "oriented", isolate = "ISO1", locus = "LSU", missing = TRUE)
)
fasta_path <- tempfile(fileext = ".fasta")
written <- pitax_export_write_fasta_records(fasta_path, seq_recs, "rich")
fasta_txt <- paste(readLines(fasta_path, warn = FALSE), collapse = "\n")
check(grepl("sequence_role=multilocus_locus", fasta_txt, fixed = TRUE), "FASTA must label sequence role")
check(length(unique(tolower(written$headers))) == 2L, "Name collisions must be disambiguated")
check(length(written$missing_lines) == 1L, "Missing sequences must be documented")
check(!grepl("AAACCC", gsub("\\s", "", fasta_txt)), "Loci must not be concatenated into one sequence")

# --- Sequence role distinction ----------------------------------------------
r1 <- make_min_result("S1", "ACGTACGT", list(target = "ITS", window = 9, min_peak_ratio = 1.5, assay_id = "assay_ITS"))
r2 <- make_min_result("S2", "ACGTACGA", list(target = "ITS", window = 21, min_peak_ratio = 2.0, assay_id = "assay_ITS"))
state_qc <- state_assign
state_qc$results <- list(S1 = r1, S2 = r2)
state_qc$summary <- rbind(r1$summary, r2$summary)
cs <- stage3_empty_consensus_set()
# Build lightweight independent analysis records without full pair algorithm.
cs$records <- list(
  ISO1_ITS = list(
    consensus_id = "ISO1_ITS", final_name = "ISO1_ITS", isolate = "ISO1", locus = "ITS",
    status = "READY", sequence = "ACGTACGT", source_read_ids = c("S1", "S2"),
    source_sequences = list(S1 = "ACGTACGT", S2 = "ACGTACGA"),
    curation = list(revision = 3L, audit_log = data.frame()),
    evidence = data.frame(Alignment_Column = 1L, Needs_Review = FALSE, Decision = "agree", stringsAsFactors = FALSE),
    alignment = list(forward_aligned = "ACGTACGT", reverse_aligned = "ACGTACGA")
  )
)
cs$summary <- data.frame(
  Consensus_ID = "ISO1_ITS", Final_Name = "ISO1_ITS", Isolate = "ISO1", Locus = "ITS",
  Status = "READY", Length = 8L, Revision = 3L, stringsAsFactors = FALSE
)
cs$schema <- "pitax-consensus-set-v1"
cs$algorithm <- "pitax-overlap-consensus-v1"
state_qc$consensus_set <- cs
# Force currentness helper path by matching stored source sequences.
snap_qc <- pitax_export_capture_snapshot(state_qc, source_stage = "qc", project_mode = "paired_consensus", app_version = "3.4.0")
# stage3_consensus_set_is_current compares result$seq to source_sequences: keep them equal.
trimmed <- pitax_export_extract_qc_trimmed_sequences(snap_qc, c("S1"))
analysis <- pitax_export_extract_consensus_analysis_sequences(snap_qc, c("S1"))
check(identical(trimmed$S1$sequence_role, "trimmed_curated_read"), "Trimmed role must be labeled")
check(identical(analysis$ISO1_ITS$sequence_role, "consensus"), "Consensus role must be labeled")

# --- Full FASTA names (per-read default), legacy rename-map precedence -----
sel_read <- pitax_export_default_selection("qc", "stage")
check(identical(sel_read$sequence_layout, "per_record"), "Per-read FASTA should be the default export layout")
sel_read$component_ids <- "qc_trimmed_sequences"
sel_read$include_plots <- FALSE
plan_read <- pitax_export_build_plan(snap_qc, sel_read)
read_paths <- vapply(Filter(function(x) identical(x$kind, "fasta"), plan_read$files), function(x) x$path, character(1))
check(setequal(read_paths, c("sequences/qc_trimmed_sequences/ISO1_ITS_F.fasta",
                             "sequences/qc_trimmed_sequences/ISO1_ITS_R.fasta")),
      "Per-record export paths must retain the complete final read names")
zip_read <- tempfile(fileext = ".zip")
pitax_export_execute(snap_qc, plan_read, zip_read)
read_dir <- tempfile("read_fasta_zip_")
dir.create(read_dir)
utils::unzip(zip_read, exdir = read_dir)
fasta_read <- file.path(read_dir, "sequences", "qc_trimmed_sequences", "ISO1_ITS_R.fasta")
check(file.exists(fasta_read), "ZIP must include the full-name per-read FASTA file")
check(identical(readLines(fasta_read, n = 1L), ">ISO1_ITS_R"),
      "The per-read FASTA header must use the complete final read name")

sel_grouped <- sel_read
sel_grouped$sequence_layout <- "per_isolate"
plan_grouped <- pitax_export_build_plan(snap_qc, sel_grouped)
grouped_paths <- vapply(Filter(function(x) identical(x$kind, "fasta"), plan_grouped$files), function(x) x$path, character(1))
check(identical(grouped_paths, "sequences/qc_trimmed_sequences/ISO1.fasta"),
      "Explicit per-isolate grouped exports should retain their isolate-only package filename")

snap_legacy <- snap_qc
snap_legacy$read_assignments$Final_Name[2] <- "320_TEF_R"
snap_legacy$rename$New_name[2] <- "320"
check(identical(pitax_export_final_name(snap_legacy, "S2"), "320_TEF_R"),
      "Explicit Final_Name must override an outdated isolate-only rename map")
check(identical(pitax_export_extract_qc_trimmed_sequences(snap_legacy, "S2")$S2$final_name, "320_TEF_R"),
      "Read sequence export must inherit the full explicit final name"

# --- Reusable, per-record custom FASTA headers -------------------------------
conf <- pitax_header_default_config()
check(identical(pitax_header_coerce_config(NULL)$fields, c("Organism", "Experiment", "Seq_Date")),
      "Legacy schema-6 projects must safely default custom header metadata")
check(identical(pitax_header_record_key("blast_query_sequences", "ISO1_ITS"),
                pitax_header_record_key("consensus_analysis_sequences", "ISO1_ITS")),
      "BLAST and analysis exports of one record must share its annotation key")
for (key in c("read::S1", "read::S2", "analysis::ISO1_ITS")) {
  conf <- pitax_header_set_value(conf, key, "Organism", "E_necator")
  conf <- pitax_header_set_value(conf, key, "Experiment", "Fungicide resistance X")
  conf <- pitax_header_set_value(conf, key, "Seq_Date", "2026-10-28")
}
rows <- pitax_header_edit_table(snap_qc, conf)
check(all(c("read::S1", "read::S2", "analysis::ISO1_ITS") %in% rows$Record_Key),
      "Header editor must enumerate stable read and consensus keys")
check(identical(rows$Experiment[rows$Record_Key == "read::S2"], "Fungicide resistance X"),
      "Metadata must be associated with the exact read key")
custom_records <- pitax_header_apply_to_records(
  pitax_export_extract_qc_trimmed_sequences(snap_qc), "qc_trimmed_sequences", conf)
check(identical(custom_records$S2$custom_header,
                'ISO1_ITS_R organism=E_necator experiment="Fungicide resistance X" date=2026-10-28'),
      "Compiled header must expand every placeholder without changing Final_Name")
check(identical(pitax_header_validate_template(">{Isolate} | {Experiment}", conf$fields),
                "{Isolate} | {Experiment}"), "Leading '>' must not be duplicated")
err_field <- tryCatch({pitax_header_validate_field("Isolate", conf$fields); NULL},
                      error = function(e) conditionMessage(e))
check(!is.null(err_field), "Reserved biological identity columns cannot become custom fields")
err_unknown <- tryCatch({pitax_header_validate_template("{Oops}", conf$fields); NULL},
                        error = function(e) conditionMessage(e))
check(!is.null(err_unknown), "Unknown custom-header tokens must fail visibly")
err_line <- tryCatch({pitax_header_set_value(conf, "read::S1", "Experiment", "bad\n>record"); NULL},
                     error = function(e) conditionMessage(e))
check(!is.null(err_line), "Metadata must not allow newline/FASTA injection")
err_missing <- tryCatch({
  pitax_header_render(custom_records$S2, "qc_trimmed_sequences", "S42", conf); NULL
}, error = function(e) conditionMessage(e))
check(!is.null(err_missing), "A missing metadata row must not silently export empty required values")

custom_sel <- pitax_export_default_selection("qc", "stage")
custom_sel$component_ids <- "qc_trimmed_sequences"
custom_sel$sequence_header <- "custom"
custom_sel$header_config <- conf
custom_sel$sequence_layout <- "per_record"
custom_sel$include_plots <- FALSE
custom_plan <- pitax_export_build_plan(snap_qc, custom_sel)
custom_zip <- tempfile(fileext = ".zip")
pitax_export_execute(snap_qc, custom_plan, custom_zip)
custom_dir <- tempfile("custom_fasta_")
dir.create(custom_dir)
utils::unzip(custom_zip, exdir = custom_dir)
custom_fasta <- file.path(custom_dir, "sequences", "qc_trimmed_sequences", "ISO1_ITS_R.fasta")
check(file.exists(custom_fasta), "Custom export must preserve normal full-name FASTA filenames")
check(identical(readLines(custom_fasta, n = 1L),
                '>ISO1_ITS_R organism=E_necator experiment="Fungicide resistance X" date=2026-10-28'),
      "Custom FASTA contents must have the rendered header exactly once")
manifest <- utils::read.csv(file.path(custom_dir, "sequences", "FASTA_header_manifest.csv"),
                            stringsAsFactors = FALSE)
check(nrow(manifest) == 2L && "Header" %in% names(manifest),
      "Custom export must include an auditable mapping of full generated headers")
run_custom <- paste(readLines(file.path(custom_dir, "RUN_INFO.txt")), collapse = "\n")
check(grepl("fasta_header_template:", run_custom, fixed = TRUE),
      "Custom header template must be preserved in export provenance")

paired_recs <- pitax_export_extract_consensus_analysis_sequences(snap_qc)
check(identical(pitax_header_apply_to_records(
  paired_recs, "consensus_analysis_sequences", conf)$ISO1_ITS$custom_header,
  'ISO1_ITS organism=E_necator experiment="Fungicide resistance X" date=2026-10-28'),
  "The same template must work for consensus FASTA exports")

dup <- conf
dup$draft_template <- "{Isolate} organism={Organism}"
dup_err <- tryCatch({
  pitax_header_apply_to_records(
    pitax_export_extract_qc_trimmed_sequences(snap_qc), "qc_trimmed_sequences", dup)
  NULL
}, error = function(e) conditionMessage(e))
check(!is.null(dup_err) && grepl("Duplicate", dup_err),
      "Non-unique primary FASTA identifiers must fail instead of silently colliding")

# Save/Load serialization must retain the user template and metadata table.
serialized <- tempfile(fileext = ".rds")
saveRDS(conf, serialized)
reloaded <- pitax_header_coerce_config(readRDS(serialized))
check(identical(reloaded$values, conf$values) &&
      identical(reloaded$draft_template, conf$draft_template),
      "FASTA header configuration must survive a project RDS round-trip")

# --- BLAST partial / waiting / failed / stale --------------------------------
state_blast <- state_qc
state_blast$blast_jobs <- data.frame(
  final_name = c("ISO1_ITS", "ISO1_ITS", "ISO1_ITS"),
  original_name = c("ISO1_ITS", "ISO1_ITS", "ISO1_ITS"),
  rid = c("RID_WAIT", "RID_FAIL", "RID_STALE"),
  rtoe = c("10", "10", "10"),
  database = c("nt", "nt", "nt"),
  hitlist_size = c(50L, 50L, 50L),
  consensus_revision = c(3L, 3L, 1L),
  status = c("WAITING", "FAILED", "STALE"),
  submitted_at = c("", "", ""), last_checked_at = c("", "", ""),
  auto_poll_enabled = c(FALSE, FALSE, FALSE), auto_poll_attempts = c(0L, 0L, 0L),
  next_poll_at = c("", "", ""), earliest_retrieve_at = c("", "", ""),
  manual_retrieval_required = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
state_blast$blast_hits <- data.frame()
snap_blast <- pitax_export_capture_snapshot(state_blast, source_stage = "blast", app_version = "3.4.0")
av_jobs <- pitax_export_assess_component("blast_jobs", snap_blast, list(sample_mode = "all"))
av_hits <- pitax_export_assess_component("blast_hits", snap_blast, list(sample_mode = "all"))
check(av_jobs$status %in% c("available", "partial"), "BLAST jobs with mixed statuses remain exportable as evidence")
check(av_hits$status %in% c("partial", "not_created"), "Waiting/failed BLAST without hits is partial/not created, not taxonomy")

# --- Ready beside not-ready samples -----------------------------------------
state_partial <- state_qc
state_partial$results <- list(S1 = r1)  # S2 not processed
state_partial$summary <- r1$summary
snap_partial <- pitax_export_capture_snapshot(state_partial, source_stage = "qc", app_version = "3.4.0")
sel_partial <- pitax_export_default_selection("qc", "stage")
sel_partial$sample_mode <- "selected"
sel_partial$selected_samples <- "S1"
sel_partial$component_ids <- c("qc_summary", "qc_trimmed_sequences")
plan_partial <- pitax_export_build_plan(snap_partial, sel_partial)
check("qc_summary" %in% plan_partial$included_components, "Ready sample export must not be blocked by other missing samples outside selection")

# --- Active vs historical curation ------------------------------------------
active <- pitax_export_extract_qc_active_curation(snap_qc, "S1")
hist <- pitax_export_extract_qc_historical_audit(snap_qc, "S1")
check(nrow(active) >= 1L && identical(active$Export_class[1], "active_curation"), "Active curation class must be marked")
check(nrow(hist) >= 1L && identical(hist$Export_class[1], "historical_audit"), "Historical audit class must be marked")

# --- Per-read settings and optional plot failure documentation --------------
check(!identical(r1$processing_settings$window, r2$processing_settings$window), "Fixture uses distinct per-read settings")
sel_plots <- pitax_export_default_selection("qc", "stage")
sel_plots$component_ids <- c("qc_metrics_plots", "qc_chromatograms", "qc_summary")
sel_plots$include_plots = TRUE
sel_plots$selected_samples <- c("S1", "S2")
sel_plots$sample_mode <- "selected"
snap_plots <- pitax_export_capture_snapshot(state_qc, source_stage = "qc", app_version = "3.4.0")
# Break chromatogram inputs for S2 to force optional failure.
snap_plots$results$S2$peak_pos <- NULL
snap_plots$results$S2$trace <- NULL
plan_plots <- pitax_export_build_plan(snap_plots, sel_plots)
zip_plots <- tempfile(fileext = ".zip")
# Ensure currentness doesn't block: neutralize consensus components.
plan_plots$included_components <- setdiff(plan_plots$included_components, grep("^consensus_", plan_plots$included_components, value = TRUE))
plan_plots$selection$component_ids <- plan_plots$included_components
plan_plots <- pitax_export_build_plan(snap_plots, plan_plots$selection)
# Avoid consensus currentness by not selecting consensus comps.
res_plots <- tryCatch(pitax_export_execute(snap_plots, plan_plots, zip_plots), error = function(e) e)
if (inherits(res_plots, "error")) {
  # If package failed for required table reasons, surface clearly.
  fail("Plot package execution failed unexpectedly: ", conditionMessage(res_plots))
}
check(file.exists(zip_plots), "QC plot package should still be produced when optional plots fail")
unzip_dir <- tempfile("plots_unzip_")
dir.create(unzip_dir)
utils::unzip(zip_plots, exdir = unzip_dir)
info2 <- paste(readLines(file.path(unzip_dir, "RUN_INFO.txt"), warn = FALSE), collapse = "\n")
check(grepl("processing_settings_by_source:", info2, fixed = TRUE), "RUN_INFO must record per-read processing settings")
check(grepl("S1:.*window=9", info2) || grepl("window=9", info2, fixed = TRUE), "RUN_INFO must retain read-specific settings")

# --- Stale active BLAST blocked from active hit export ----------------------
hits_active <- pitax_export_extract_blast_hits(snap_blast, NULL)
check(nrow(hits_active) == 0L, "Stale/waiting BLAST must not invent active hit rows")

# --- Excel split / combined layouts ----------------------------------------
sel_excel <- pitax_export_default_selection("rename", "stage")
sel_excel$component_ids <- "assign_identity_table"
sel_excel$excel_layout <- "sheet_per_component"
plan_excel <- pitax_export_build_plan(snap_assign, sel_excel)
zip_excel <- tempfile(fileext = ".zip")
pitax_export_execute(snap_assign, plan_excel, zip_excel)
check(file.exists(zip_excel), "Excel layout package must write")

# --- Catalog completeness ---------------------------------------------------
ids <- pitax_export_component_ids()
check(length(ids) >= 30L, "Catalog should cover all stage families")
check(all(c(
  "upload_file_inventory", "assign_identity_table", "qc_chromatograms",
  "consensus_analysis_sequences", "blast_jobs", "taxonomy_summary", "multilocus_evidence"
) %in% ids), "Catalog missing required component families")

# --- Shiny startup + dialog wiring ------------------------------------------
if (requireNamespace("shiny", quietly = TRUE)) {
  app_env <- new.env(parent = globalenv())
  sys.source("app.R", envir = app_env)
  check(is.function(app_env$server), "app.R must define server()")
  check(is.function(app_env$pitax_export_build_plan), "Export planner must be loaded through bootstrap")
  check(is.function(app_env$pitax_export_catalog), "Export catalog must be loaded through bootstrap")
  shiny::testServer(app_env$server, {
    session$setInputs(pipeline_step = "upload")
    session$setInputs(open_export_from_stage = 1L)
    session$flushReact()
    stopifnot(!is.null(export_dialog_state$snap))
    stopifnot(!is.null(export_dialog_state$plan))
    stopifnot(identical(export_dialog_state$stage_id, "upload"))
    av_qc <- export_dialog_state$plan$availability$qc_summary
    stopifnot(!is.null(av_qc))
    stopifnot(identical(av_qc$status, "not_created"))
  })
}

message("Shared export engine smoke tests passed.")
