ui <- fluidPage(
  tags$head(
    tags$title("PITAX — Taxonomic Identification Tool"),
    tags$link(rel = "stylesheet", type = "text/css", href = "pitax.css"),
    tags$script(src = "pitax.js")
  ),

  div(id = "shiny_activity_bar", `aria-hidden` = "true"),
  div(id = "app_loading_overlay", class = "app-loading-overlay", `aria-hidden` = "true",
      div(class = "app-loading-card",
          div(class = "app-loading-spinner"),
          div(id = "app_loading_text", "Loading workspace…")
      )
  ),

  div(class = "pipeline-container",
    div(class = "app-header app-shell-header",
      if (PITAX_LOGO_AVAILABLE)
        tags$img(src = "logo.png", class = "app-brand-logo", alt = "PITAX — Taxonomic Identification Tool")
      else
        tagList(
          div(class = "app-brand-mark", icon("flask")),
          div(class = "app-brand-copy", h2("PITAX"), p("Taxonomic Identification Tool"))
        ),
      div(class = "app-brand-copy",
        p("Sanger sequence analysis · quality review · BLAST · taxonomic interpretation")
      ),
      div(class = "app-version-badge", paste0("v", APP_VERSION))
    ),

    div(class = "project-bar project-session-card",
      div(class = "project-bar-title", icon("folder-open"), span("Project session")),
      downloadButton("save_project", "Save project", class = "btn-project"),
      fileInput("load_project", NULL, multiple = FALSE, accept = c(".sangerproject", ".rds"),
                buttonLabel = "Load project", placeholder = "No project selected"),
      div(class = "project-status", uiOutput("project_status"))
    ),

    tabsetPanel(id = "pipeline_step", type = "tabs",

      # --------------------------------------------------------
      # 1. Upload
      # --------------------------------------------------------
      tabPanel("1 · Upload", value = "upload",
        stage_heading("upload", "Upload chromatograms", "Keep the sequencer barcode and source filename unchanged. Biological identity is assigned later in Rename.", "Step 1 of 9"),
        stage_topbar(
          div(class = "stage-topbar-spacer"),
          actionButton("to_settings", "Continue to Assay", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "stage-grid stage-grid-upload",
          div(class = "panel-box upload-drop-card",
            card_title("Raw AB1 files", "Select one or more .ab1 chromatogram files. The original files remain the immutable source for all later QC and curation.", "upload"),
            fileInput("ab1_files", NULL, multiple = TRUE,
                      accept = c(".ab1", ".AB1"), buttonLabel = "Select AB1 files", placeholder = "or drag files here"),
            div(class = "compact-hint", icon("info-circle"), span("Multiple files can be processed in one run."))
          ),
          div(class = "panel-box stage-table-card",
            card_title("Files in this run", "Review the uploaded filenames before continuing.", "list"),
            DTOutput("uploaded_files_table")
          )
        ),
        div(class = "panel-box",
            card_title("Project read model", "Choose whether every chromatogram is an independent sequence or whether explicit Forward/Reverse reads should be paired into one isolate-level consensus.", "sitemap"),
          div(class = "project-mode-options",
            radioButtons(
              "project_mode", NULL,
              choices = c(
                "Simple reads — no Forward/Reverse matching" = "simple",
                "Paired reads — build Forward/Reverse consensus" = "paired_consensus"
              ),
              selected = "simple", inline = FALSE
            )
          ),
          div(class = "compact-hint", icon("info-circle"), span("Direction is still required in both modes so Reverse reads can be oriented correctly. The mode can be changed before Stage 3 is built."))
        ),
        pipeline_stage_footer(1)
      ),
      # --------------------------------------------------------
      # 2. Assay and trimming settings
      # --------------------------------------------------------
      tabPanel("2 · Assay", value = "settings",
        stage_heading("sliders", "Assay setup", "Define run-level locus/primer defaults and the automatic trimming rules. Trimming starts only after Rename.", "Step 2 of 9"),
        stage_topbar(
          actionButton("back_upload", "Back", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_rename", "Continue to Rename", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "stage-grid stage-grid-2",
          div(class = "panel-box settings-card",
            card_title("Assay information", "Primer metadata is retained for provenance. Experimental primer-site mapping is optional and does not control trimming.", "flask"),
            selectInput("target", "Target / Gene",
                        c("ITS", "LSU", "TEF1 / EF1-alpha", "RPB2", "Beta-tubulin", "CYP51", "SDHB", "IGS", "Other"),
                        selected = "ITS"),
            div(class = "form-grid-2",
              textInput("forward_primer", "Forward primer name", ""),
              textInput("reverse_primer", "Reverse primer name", "")
            ),
            div(class = "form-grid-2",
              textInput("forward_primer_seq", "Forward primer sequence (5'→3')", ""),
              textInput("reverse_primer_seq", "Reverse primer sequence (5'→3')", "")
            ),
            radioButtons(
              "sequencing_primer", "Primer used for this Sanger read",
              choices = c("Forward" = "Forward", "Reverse" = "Reverse", "Unknown / infer" = "Unknown"),
              selected = "Forward", inline = TRUE
            ),
            div(class = "inline-control-row",
              checkboxInput("enable_primer_mapping", "Experimental primer-site mapping", value = FALSE),
              info_tip("Primer-site inference can be unreliable because the sequencing primer itself may not appear in the base calls and the beginning of a Sanger read is often noisy.")
            ),
            div(class = "form-grid-2",
              numericInput("expected_amplicon_len", "Expected amplicon length (bp)", 650, min = 1),
              numericInput("absolute_max_base_index", "Maximum sequence position (bp)", 680, min = 50)
            )
          ),
          div(class = "panel-box settings-card",
            card_title("Automatic trimming", "These parameters define good-start detection and sustained signal-collapse detection. Full definitions are available in Help / About.", "cut"),
            div(class = "form-grid-2",
              numericInput("window", "Window size", 25, min = 5),
              numericInput("min_peak_ratio", "Minimum peak ratio", 3, min = 0, step = 0.1),
              numericInput("min_relative_signal", "Minimum relative signal", 0.20, min = 0, max = 1, step = 0.01),
              numericInput("min_len_before_collapse", "Minimum length before collapse search", 350, min = 1),
              numericInput("bad_run_windows", "Consecutive bad windows", 12, min = 1),
              numericInput("min_usable_len", "Minimum usable trimmed length", 400, min = 1)
            )
          )
        ),
        pipeline_stage_footer(2)
      ),
      # --------------------------------------------------------
      # 4. QC, chromatogram & sequence preview
      # --------------------------------------------------------
      tabPanel("4 · Trim & QC", value = "qc",
        stage_heading("bar-chart", "Trimming results, QC & curation", "Review the completed trim, inspect renamed chromatograms, and document manual sequence curation.", "Step 4 of 9"),
        stage_topbar(
          actionButton("back_rename_from_qc", "Back to Rename", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_consensus", "Continue to Consensus", icon = icon("arrow-right"), class = "btn-primary")
        ),
        uiOutput("qc_summary_cards"),
        div(class = "panel-box stage-table-card",
          card_title("Trimming results", "Run-level summary of automatic trimming for all uploaded chromatograms.", "table"),
          DTOutput("summary_table")
        ),
        div(class = "qc-inspection-grid",
          div(class = "panel-box qc-sidebar-card",
            card_title("Sequence inspection", "Select a sample to inspect its current curated sequence, QC metrics and optional primer mapping.", "search"),
            selectInput("inspect_sample", "Sample", choices = NULL),
            DTOutput("sequence_metrics"),
            conditionalPanel(
              condition = "input.enable_primer_mapping === true",
              div(class = "subsection-divider"),
              div(class = "subsection-title", "Primer-site mapping", info_tip("Experimental mapping only; it is not used as an automatic trim boundary.")),
              DTOutput("primer_match_table")
            )
          ),
          div(class = "qc-main-stack",
            div(class = "panel-box compact-plot-card",
              card_title("Read overview", "Automatic trim boundaries and optional primer context.", "arrows-h"),
              plotOutput("amplicon_overview", height = "190px"),
              uiOutput("expected_amplicon_note")
            ),
            div(class = "panel-box chromatogram-card",
              card_title("Chromatogram", "X axis = raw called-base position. Zoom, pan, or click a flagged position to inspect it.", "line-chart"),
              plotly::plotlyOutput("chromatogram_plot", height = "540px")
            )
          )
        ),
        div(class = "panel-box stage1-evidence-card",
          card_title("Stage 1 · AB1 evidence audit", "Observational audit of the established PITAX v2 read model against raw ABIF primary-call coordinates, basecaller quality and canonical A/C/G/T trace evidence. This panel does not alter trimming, curation, BLAST or taxonomy.", "microscope"),
          tags$details(
            class = "stage1-audit-details",
            tags$summary(class = "stage1-audit-toggle", "Open evidence audit"),
            div(class = "status-note",
                tags$strong("Validation mode only. "),
                "The active v2.14.2 trimming/QC path is still the decision path. The current build keeps the same-length PCON-only comparison observational; it does not apply that proposal to the processed sequence, FASTA or BLAST output."),
            div(class = "subsection-title", "Run-level audit"),
            DTOutput("ab1_evidence_run_table"),
            div(class = "blast-action-row",
              downloadButton("download_ab1_evidence_run", "Download run audit CSV", class = "btn-default"),
              downloadButton("download_ab1_evidence_detail", "Download selected-base audit CSV", class = "btn-default")
            ),
            div(class = "subsection-divider"),
            div(class = "subsection-title", "Selected sample · per-base evidence"),
            uiOutput("ab1_evidence_selected_note"),
            div(class = "subsection-title", "Legacy auto trim vs PCON-only comparison"),
            DTOutput("ab1_trim_comparison_table"),
            DTOutput("ab1_evidence_detail_table")
          )
        ),
        div(class = "panel-box peak-review-card",
          card_title("Ambiguous peak review", "Flags identify channel competition or locally non-dominant calls. Left-click centers the trace; right-click opens curation actions.", "flag"),
          div(class = "curation-toolbar",
            radioButtons(
              "peak_flag_scope", "Positions",
              choices = c("Trimmed sequence" = "trimmed", "Entire raw read" = "raw"),
              selected = "trimmed", inline = TRUE
            ),
            checkboxInput("show_peak_flags", "Show markers", value = TRUE),
            checkboxInput("show_reviewed_flags", "Show reviewed", value = FALSE),
            div(class = "curation-actions",
              actionButton("auto_correct_preview", "Auto-correct high-confidence", icon = icon("magic"), title = "Preview positions that meet the current auto-correction criteria. No base is changed before confirmation."),
              actionButton("auto_correct_settings", "Criteria", icon = icon("cog"), title = "Edit the thresholds used to propose automatic base corrections."),
              actionButton("curation_undo", "Undo", icon = icon("undo"), title = "Undo the most recent curation transaction for this sample."),
              actionButton("curation_redo", "Redo", icon = icon("repeat"), title = "Redo the most recently undone curation transaction."),
              actionButton("curation_history", "History", icon = icon("history"), title = "Open the full per-sample manual curation audit log."),
              actionButton("curation_reset", "Reset to auto", icon = icon("refresh"), title = "Restore automatic trim/base calls as the active curated sequence. The reset is logged and undoable.")
            )
          ),
          uiOutput("ambiguous_peak_summary"),
          DTOutput("ambiguous_peak_table")
        ),
        div(class = "panel-box",
          card_title("QC signal metrics", "Called-base signal and called-peak / second-peak ratio. These plots are retained in QC exports.", "area-chart"),
          plotOutput("qc_plot", height = "650px")
        ),
        conditionalPanel(
          condition = "input.enable_primer_mapping === true",
          div(class = "panel-box",
            card_title("Experimental primer alignment", "Detailed primer alignment output for manual inspection.", "code"),
            verbatimTextOutput("primer_alignment")
          )
        ),
        div(class = "panel-box sequence-preview-card",
          card_title("Curated sequence", "This exact sequence retains its resolved sample name and is carried forward to export and BLAST.", "file-text"),
          textAreaInput("trimmed_sequence_preview", NULL, value = "", rows = 8, width = "100%")
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy",
            card_title("Checkpoint B · Renamed and curated sequences", "Save the resolved naming and QC/curation state before export.", "save")
          ),
          downloadButton("download_trim_checkpoint", "Download checkpoint ZIP")
        ),
        pipeline_stage_footer(4)
      ),
      # --------------------------------------------------------
      # 3. Rename
      # --------------------------------------------------------
      tabPanel("3 · Rename & Assign", value = "rename",
        stage_heading("tags", "Rename and assign read identity", "Assign isolate, gene and Forward/Reverse direction here, after Upload and Assay and before trimming. The upload barcode remains unchanged.", "Step 3 of 9"),
        stage_topbar(
          actionButton("back_settings_from_rename", "Back to Assay", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("run_trimming", "Start trimming", icon = icon("play"), class = "btn-primary")
        ),
        div(class = "stage-grid stage-grid-rename",
          div(class = "panel-box",
            card_title("Assignment key", "Import XLSX/CSV with old_id, isolate, locus (or gene), and direction columns. old_id matches the uploaded barcode; prefix matching is supported.", "key"),
            fileInput("rename_key_file", NULL, multiple = FALSE, accept = c(".xlsx", ".csv"), buttonLabel = "Choose assignment key", placeholder = "XLSX or CSV"),
            actionButton("apply_rename_key", "Apply assignment key", icon = icon("check"), class = "btn-primary"),
            downloadButton("download_assignment_key_template", "Download key template"),
            uiOutput("rename_key_status")
          ),
          div(class = "panel-box",
            card_title("Batch assignment", "Apply isolate edits, gene and direction to checked rows. If no rows are checked, the action applies to all uploaded reads.", "edit"),
            div(class = "form-grid-2",
              textInput("batch_isolate_prefix", "Isolate prefix", ""),
              textInput("batch_isolate_suffix", "Isolate suffix", "")
            ),
            div(class = "form-grid-2",
              textInput("batch_isolate_find", "Find in isolate", ""),
              textInput("batch_isolate_replace", "Replace with", "")
            ),
            div(class = "form-grid-2",
              textInput("batch_locus", "Set gene / locus", ""),
              selectInput("batch_direction", "Set direction", c("No change" = "", "Forward" = "Forward", "Reverse" = "Reverse"), selected = "")
            ),
            actionButton("apply_assignment_batch", "Apply to checked / all", class = "btn-primary"),
            actionButton("reset_assignments", "Clear checked / all")
          )
        ),
        div(class = "panel-box stage-table-card",
          card_title("Final read / FASTA names and biological identity", "The generated name is shown beside the immutable barcode. Edit Isolate and Gene, then choose Forward or Reverse from the list. Changes are committed together, without refreshing the table after every cell.", "list-alt"),
          uiOutput("assignment_editor"),
          div(class = "assignment-editor-actions",
            actionButton("save_assignment_edits", "Apply table changes", icon = icon("check"), class = "btn-primary"),
            div(class = "compact-hint", icon("info-circle"), span("The table stays stable while you edit; generated names update after Apply."))
          ),
          uiOutput("rename_validation")
        ),
        div(class = "panel-box stage2-assignment-card",
          card_title("Stage 2 · Architecture preview", "Built directly from the explicit identity fields. Project mode determines whether matching Forward/Reverse reads will later be merged or remain independent.", "sitemap"),
          uiOutput("architecture_summary")
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint A · Renamed and assigned reads", "Save resolved read names and biological identity before trimming.", "save")),
          downloadButton("download_rename_checkpoint", "Download checkpoint ZIP")
        ),
        pipeline_stage_footer(3)
      ),
      # --------------------------------------------------------
      # 5. Forward / Reverse consensus
      # --------------------------------------------------------
      tabPanel("5 · Analysis Sequence", value = "consensus",
        stage_heading("random", "Stage 3 · Analysis sequence", "Create the sequence used for export and BLAST according to the project read model selected at Upload.", "Step 5 of 9"),
        stage_topbar(
          actionButton("back_qc_from_consensus", "Back to QC", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("build_consensus", "Build / rebuild", icon = icon("cogs"), class = "btn-success"),
          actionButton("to_export", "Continue to Export", icon = icon("arrow-right"), class = "btn-primary")
        ),
        uiOutput("consensus_mode_note"),
        div(class = "stage-grid stage-grid-2",
          div(class = "panel-box settings-card",
            card_title("Conservative overlap rules", "These are run-level review thresholds, not universal biological constants. A weak overlap is blocked rather than joined.", "sliders"),
            div(class = "form-grid-2",
              numericInput("consensus_min_overlap", "Minimum overlap (bp)", 40, min = 10, step = 5),
              numericInput("consensus_min_identity", "Minimum overlap identity (%)", 85, min = 0, max = 100, step = 1),
              numericInput("consensus_quality_delta", "Minimum quality advantage", 10, min = 0, step = 1),
              numericInput("consensus_strong_quality", "Strong base quality", 20, min = 0, step = 1)
            ),
            div(class = "compact-hint", "A mismatch is resolved automatically only when one base has sufficiently stronger basecaller quality. Otherwise PITAX retains an IUPAC ambiguity call and requires review.")
          ),
          div(class = "panel-box",
            card_title("Consensus gate", "Single reads remain valid representatives. Paired reads must have a reliable overlap and no unresolved review positions before downstream analysis.", "check-circle"),
            uiOutput("consensus_gate_status"),
            div(class = "button-row",
              downloadButton("download_consensus_fasta", "Analysis FASTA"),
              downloadButton("download_consensus_checkpoint", "Consensus checkpoint ZIP")
            )
          )
        ),
        div(class = "panel-box stage-table-card",
          card_title("Analysis sequence summary", "Simple mode keeps one oriented sequence per read. Paired mode creates one row per Isolate–Locus and labels unpaired reads as single-read representatives.", "table"),
          DTOutput("consensus_summary_table")
        ),
        div(class = "panel-box consensus-review-card",
          card_title("Consensus Review & Curation", "Resolve only evidence conflicts that need a biological decision. Forward, Reverse and the automatic IUPAC call remain linked to their raw chromatograms and every action creates an auditable revision.", "edit"),
          div(class = "consensus-review-grid",
            div(
              selectInput("consensus_review_position", "Review position", choices = NULL),
              uiOutput("consensus_review_call_ui"),
              textInput("consensus_review_note", "Review note", placeholder = "Optional rationale"),
              div(class = "consensus-review-actions",
                actionButton("apply_consensus_review", "Apply decision", icon = icon("check"), class = "btn-primary"),
                actionButton("consensus_review_undo", "Undo", icon = icon("undo")),
                actionButton("consensus_review_redo", "Redo", icon = icon("repeat")),
                actionButton("consensus_review_history", "History", icon = icon("history"))
              ),
              uiOutput("consensus_review_state")
            ),
            div(class = "consensus-trace-grid",
              plotly::plotlyOutput("consensus_forward_conflict_plot", height = "315px"),
              plotly::plotlyOutput("consensus_reverse_conflict_plot", height = "315px")
            )
          )
        ),
        div(class = "qc-inspection-grid",
          div(class = "panel-box qc-sidebar-card",
            card_title("Consensus inspection", "Select an isolate/locus to inspect the oriented overlap and every automatic or unresolved base decision.", "search"),
            selectInput("consensus_sample", "Isolate / locus", choices = NULL),
            uiOutput("consensus_selected_status"),
            DTOutput("consensus_selected_metrics")
          ),
          div(class = "qc-main-stack",
            div(class = "panel-box sequence-preview-card",
              card_title("Oriented overlap", "The Reverse read is shown after reverse-complementing. Gaps are alignment characters only and do not modify either curated source read.", "exchange"),
              verbatimTextOutput("consensus_alignment_text")
            ),
            div(class = "panel-box sequence-preview-card",
              card_title("Analysis sequence", "This is the oriented independent read or paired consensus that proceeds to FASTA, BLAST and later Stage 4 locus profiles after the gate is green.", "file-text"),
              verbatimTextOutput("consensus_sequence_text")
            )
          )
        ),
        div(class = "panel-box stage-table-card",
          card_title("Per-column provenance", "Every consensus position records the Forward and Reverse calls, available basecaller quality, raw chromatogram positions, decision rule and review state.", "list-alt"),
          DTOutput("consensus_evidence_table")
        ),
        pipeline_stage_footer(5)
      ),
      # --------------------------------------------------------
      # 6. Export
      # --------------------------------------------------------
      tabPanel("6 · Export", value = "export",
        stage_heading("download", "Export analysis sequences", "Create working FASTA files or an auditable package containing the analysis sequence and original read evidence.", "Step 6 of 9"),
        stage_topbar(
          actionButton("back_consensus_from_export", "Back to Consensus", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_blast", "Continue to NCBI BLAST", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "panel-box export-summary-card",
          card_title("Run summary", "Final pre-export overview of the processed sequence set.", "check-circle"),
          uiOutput("export_summary")
        ),
        div(class = "export-tile-grid",
          div(class = "export-tile", div(class = "export-tile-icon", icon("search")), h4("BLAST FASTA"), p("Processed sequences ready for sequence search."), downloadButton("download_blast_fasta", "Download FASTA")),
          div(class = "export-tile", div(class = "export-tile-icon", icon("file-text")), h4("FASTA + metadata"), p("Processed sequences with run metadata in FASTA headers."), downloadButton("download_full_fasta", "Download FASTA")),
          div(class = "export-tile", div(class = "export-tile-icon", icon("table")), h4("Summary CSV"), p("Compact processing and QC summary for downstream review."), downloadButton("download_summary_csv", "Download CSV")),
          div(class = "export-tile export-tile-primary", div(class = "export-tile-icon", icon("archive")), h4("Complete results package"), p("Sequences, QC evidence, settings and curation records in one ZIP."), downloadButton("download_all_zip", "Download results ZIP"))
        ),
        pipeline_stage_footer(6)
      ),
      # --------------------------------------------------------
      # 7. NCBI BLAST
      # --------------------------------------------------------
      tabPanel("7 · NCBI BLAST", value = "blast",
        stage_heading("search", "NCBI BLAST workspace", "Submit analysis sequences, retrieve accession-level hits, and keep each RID linked to the active sequence revision.", "Step 7 of 9"),
        stage_topbar(
          actionButton("back_export", "Back to Export", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_taxonomy", "Continue to Taxonomic Summary", icon = icon("arrow-right"), class = "btn-primary"),
          actionButton("reset_pipeline", "Start new run")
        ),
        div(class = "panel-box blast-query-card",
          card_title("Query workspace", "The sequence shown here is the current curated and renamed sequence. Changing a curated sequence after BLAST marks its previous result stale.", "file-code-o"),
          div(class = "blast-query-grid",
            div(class = "blast-query-controls",
              selectInput("blast_sample", "Sequence", choices = NULL),
              div(class = "form-grid-2",
                selectInput("blast_database", "NCBI database", choices = c("core_nt", "nt"), selected = "core_nt"),
                numericInput("blast_hitlist", tagList("Results to retrieve", info_tip("Controls BLAST retrieval depth, not a voting threshold. All returned accessions are audited; only molecularly close alternatives can reduce the supported taxonomic rank.")), 25, min = 1, max = 100)
              ),
              div(class = "compact-hint", icon("info-circle"), span("More weak tail hits do not lower confidence. A close competing taxon can lower it even if it appears late in the returned list.")),
              div(class = "button-row",
                actionButton("copy_blast_sequence", "Copy sequence", icon = icon("copy")),
                actionButton("open_ncbi_blast", "Open NCBI BLAST", icon = icon("external-link")),
                downloadButton("download_selected_blast", "Selected FASTA")
              )
            ),
            div(class = "sequence-code-panel",
              div(class = "subsection-title", "Processed query sequence"),
              uiOutput("blast_sequence_preview_ui")
            )
          )
        ),
        div(class = "panel-box blast-submit-card",
          card_title("Automated submission", "Each sample receives its own NCBI Request ID (RID). Automated contacts are paced to respect NCBI service limits.", "cloud-upload"),
          div(class = "blast-action-row blast-primary-actions",
            actionButton("submit_ncbi_blast", "Submit selected", icon = icon("paper-plane"), class = "btn-primary"),
            actionButton("submit_all_ncbi_blast", "Submit all", icon = icon("paper-plane"), class = "btn-success"),
            actionButton("retrieve_ncbi_blast", "Retrieve selected", icon = icon("refresh")),
            actionButton("retrieve_all_ncbi_blast", "Retrieve all", icon = icon("download"))
          ),
          div(class = "blast-status-strip", uiOutput("blast_batch_status"), uiOutput("blast_job_status")),
          DTOutput("blast_jobs_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Preliminary top-hit overview", "One row per retrieved sequence. This is a quick orientation only; multi-hit evidence is interpreted in the next stage.", "eye"),
          DTOutput("blast_identification_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Retrieved accession-level hits", "All parsed unique accessions for the selected RID. Multiple HSPs belonging to the same accession are aggregated.", "database"),
          DTOutput("blast_hits_table"),
          tags$details(
            class = "raw-response-details",
            tags$summary("Raw NCBI response"),
            verbatimTextOutput("blast_raw_preview")
          )
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint D · BLAST workspace", "Export BLAST job metadata and accession-level hits.", "save")),
          div(class = "button-row",
            downloadButton("download_blast_jobs", "Job/results CSV"),
            downloadButton("download_blast_hits", "All BLAST hits CSV")
          )
        ),
        pipeline_stage_footer(7)
      ),
      # --------------------------------------------------------
      # 8. Taxonomic interpretation
      # --------------------------------------------------------
      tabPanel("8 · Taxonomic summary", value = "taxonomy",
        stage_heading("sitemap", "Taxonomic interpretation", "Identify the best molecular match, inspect close alternatives and report the most conservative supported taxonomic level.", "Step 8 of 9"),
        stage_topbar(
          actionButton("back_blast", "Back to NCBI BLAST", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("run_taxonomy", "Analyze selected", icon = icon("play"), class = "btn-primary"),
          actionButton("run_taxonomy_all", "Analyze all retrieved", icon = icon("tasks"), class = "btn-success"),
          actionButton("to_multilocus", "Continue to Multi-locus", icon = icon("arrow-right")),
          actionButton("reset_pipeline_tax", "Start new run")
        ),
        div(class = "taxonomy-workspace",
          div(class = "taxonomy-workspace-left",
            selectInput("tax_sample", "Sequence", choices = NULL),
            info_tip("Only sequences with retrieved BLAST hits appear here. Each NCBI accession is counted once. The decision starts from Identity + query coverage, then checks close alternative taxa, sequence evidence and reference context.")
          ),
          div(class = "taxonomy-workspace-status", uiOutput("taxonomy_status"))
        ),
        uiOutput("taxonomy_hero"),
        uiOutput("taxonomy_summary_cards"),
        div(class = "taxonomy-main-grid",
          div(class = "taxonomy-card",
            div(class = "taxonomy-card-title", "Taxonomic interpretation", info_tip("The compact decision row shows the final recommendation. Detailed decision components remain available in exports and Help / About.")),
            DTOutput("taxonomy_summary_table"),
            uiOutput("taxonomy_locus_note")
          ),
          div(class = "taxonomy-card",
            div(class = "taxonomy-card-title", "Species evidence profile", info_tip("One row per resolved species. Best Identity and coverage describe that species' strongest comparable accession; accession count is database context only, not a vote.")),
            DTOutput("taxonomy_counts_table"),
            div(class = "taxonomy-agreement-note",
              div(class = "tax-callout-icon", icon("info-circle")),
              div("Best molecular match is chosen from the near-best query-coverage band and then by Identity. Close alternatives are species with nearly the same Identity and coverage; database abundance is shown separately and does not decide the identification.")
            )
          )
        ),
        div(class = "taxonomy-full-card",
          card_title("BLAST score landscape", "Each point is one unique NCBI accession. Points are colored by genus when several genera are present, or by species when the uncertainty is within one genus.", "line-chart"),
          tags$details(class = "taxonomy-explain-details",
            tags$summary(icon("info-circle"), "How to read this graph"),
            div(class = "taxonomy-explain-body",
              "Read left to right by BLAST score rank. The connecting line shows the score landscape, while point colors show which taxa occupy that landscape. When several genera occur, colors represent genus; when all resolved hits are within one genus, colors represent species. Hover to compare Identity, query coverage and accession. Similar colors at similar heights can indicate repeated support, but the number of accessions is not treated as a majority vote."
            )
          ),
          plotly::plotlyOutput("taxonomy_score_plot", height = "390px")
        ),
        div(class = "taxonomy-full-card",
          card_title("Taxonomy-enriched BLAST hits", "Inspect accession-level evidence and the taxonomy attached to each usable BLAST hit.", "database"),
          DTOutput("taxonomy_hits_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Team identification summary", "One row per processed sequence, combining QC, BLAST and final taxonomic interpretation.", "users"),
          DTOutput("team_summary_table"),
          div(class = "button-row export-inline-actions",
            downloadButton("download_team_summary_csv", "Team summary CSV"),
            downloadButton("download_team_summary_xlsx", "Team summary Excel")
          )
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint E · Taxonomic interpretation", "Save the final taxonomic evidence and interpretation package.", "save")),
          div(class = "button-row",
            downloadButton("download_taxonomy_summary", "Summary CSV"),
            downloadButton("download_taxonomy_hits", "Enriched hits CSV"),
            downloadButton("download_taxonomy_checkpoint", "Checkpoint ZIP")
          )
        ),
        pipeline_stage_footer(8)
      ),
      # --------------------------------------------------------
      # 9. Multi-locus isolate profile
      # --------------------------------------------------------
      tabPanel("9 · Multi-locus", value = "multilocus",
        stage_heading("th", "Stage 4 · Multi-locus isolate profile", "Combine separately completed single-locus PITAX projects by explicit Isolate code while retaining every locus-specific sequence and taxonomic limitation.", "Step 9 of 9"),
        stage_topbar(
          actionButton("back_taxonomy_from_multilocus", "Back to Taxonomic Summary", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("build_multilocus_profile", "Build / rebuild profile", icon = icon("cogs"), class = "btn-success")
        ),
        div(class = "multilocus-contract",
          strong("Scientific contract: "),
          "PITAX does not vote across loci. Concordance can strengthen a supported rank; a conflicting locus remains visible and prevents a combined call. Projects are joined only by the explicit Isolate field, never by parsing a filename."
        ),
        div(class = "panel-box multilocus-visual-card",
          card_title("Multi-locus interpretation", "Navigate between isolates and compare their locus-level molecular evidence without replacing the complete audit tables below.", "dashboard"),
          div(class = "multilocus-visual-toolbar",
            selectInput("multilocus_isolate", "Isolate", choices = NULL),
            uiOutput("multilocus_overview_cards")
          ),
          uiOutput("multilocus_profile_hero"),
          div(class = "multilocus-selected-grid",
            div(class = "multilocus-plot-card",
              div(class = "subsection-title", "Identity and query coverage by locus"),
              plotly::plotlyOutput("multilocus_evidence_plot", height = "330px")
            ),
            div(
              div(class = "subsection-title", "Locus evidence cards"),
              uiOutput("multilocus_locus_cards")
            )
          )
        ),
        div(class = "multilocus-grid",
          div(class = "panel-box multilocus-source-card",
            card_title("Completed single-locus projects", "Add saved .sangerproject files from other loci. The current session can be included without saving and re-uploading it.", "folder-open"),
            checkboxInput("multilocus_include_current", "Include the current PITAX project", value = TRUE),
            fileInput("multilocus_projects", NULL, multiple = TRUE, accept = c(".sangerproject", ".rds"),
                      buttonLabel = "Add locus projects", placeholder = "Select one or more completed projects"),
            div(class = "compact-hint", icon("info-circle"), span("Each source project must contain exactly one Gene/Locus and a current Stage 3 analysis sequence. Isolate codes must match across projects.")),
            uiOutput("multilocus_gate_status")
          ),
          div(class = "panel-box stage-table-card",
            card_title("Source provenance", "Every imported project remains traceable by filename, saved version, locus and file fingerprint.", "archive"),
            DTOutput("multilocus_sources_table")
          )
        ),
        div(class = "panel-box stage-table-card",
          card_title("Isolate profiles", "Profile status reports cross-locus concordance, missing evidence or conflict. It never counts loci as majority votes.", "sitemap"),
          DTOutput("multilocus_profiles_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Per-locus evidence", "One row per Isolate–Locus with sequence revision, BLAST/taxonomy result and reference provenance retained separately.", "list-alt"),
          DTOutput("multilocus_evidence_table")
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint F · Multi-locus profile", "Export the isolate profile, every locus-specific evidence row and all included sequences.", "save")),
          div(class = "button-row",
            downloadButton("download_multilocus_profiles", "Profiles CSV"),
            downloadButton("download_multilocus_evidence", "Per-locus evidence CSV"),
            downloadButton("download_multilocus_fasta", "Multi-locus FASTA"),
            downloadButton("download_multilocus_checkpoint", "Checkpoint ZIP")
          )
        ),
        pipeline_stage_footer(9)
      ),
      # --------------------------------------------------------
      # Help / About (documentation, not a pipeline stage)
      # --------------------------------------------------------
      tabPanel("Help / About", value = "help",
        panel_box(
          section_title(paste0("Sanger Sequence Pipeline v", APP_VERSION)),
          p(class="about-lead",
            "Documentation for the laboratory workflow, the BLAST/taxonomy interpretation logic, and the scientific sources used to guide the application. Published evidence and application-specific heuristics are labeled separately."),
          div(class="help-flow",
              "AB1 upload + project mode  →  Assay settings  →  Rename & read assignment  →  Trim & QC  →  Analysis sequence  →  Export  →  NCBI BLAST  →  Taxonomic interpretation  →  Multi-locus profile")
        ),

        div(class="about-section",
          tabsetPanel(
            id="about_tabs",
            type="tabs",

            tabPanel("Overview",
              div(class="help-card",
                h3("What the application does"),
                p("The application converts raw Sanger AB1 chromatograms into auditable processed sequences and then supports sequence-based identification using NCBI BLAST and multi-hit taxonomic interpretation."),
                div(class="about-callout",
                  strong("Core principle: "),
                  "the application separates what the database hits agree on from how strong the sequence evidence is. A weak species-level result can therefore fall back to a supported genus instead of becoming automatically Unresolved."
                )
              ),
              fluidRow(
                column(6,
                  div(class="help-card",
                    h3("Workflow"),
                    p(strong("1. Upload"), " — raw AB1 chromatograms and the project read model; source barcodes remain unchanged."),
                    p(strong("2. Assay"), " — run-level primer defaults and trimming parameters; no trimming starts yet."),
                    p(strong("3. Rename & Assign"), " — import a key, batch-edit or manually assign explicit identity fields and PITAX-generated <Isolate>_<Locus>_<F/R> names."),
                    p(strong("4. Trim & QC"), " — start trimming explicitly, then review Quality Control plots, chromatograms and processed sequences."),
                    p(strong("5. Analysis Sequence"), " — keep independent reads in Simple mode or create an auditable F/R consensus in Paired mode; Reverse reads are oriented without changing the source."),
                    p(strong("6. Export"), " — isolate-level FASTA plus source-read QC and consensus evidence checkpoints."),
                    p(strong("7. NCBI BLAST"), " — submit isolate-level sequences and retrieve accession-level hits."),
                    p(strong("8. Taxonomic summary"), " — compare competitive hits and report identification plus confidence."),
                    p(strong("9. Multi-locus profile"), " — import separately completed locus projects, join them by explicit Isolate code and retain concordance or conflict without flat voting.")
                  )
                ),
                column(6,
                  div(class="help-card",
                    h3("Auditability"),
                    p("Read assignments, isolate/locus/read links, QC plots, automatic trimming parameters, manual curation audit log, BLAST hits, taxonomy-enriched hits, decision components and application version are retained in exports."),
                    p(strong("Save / Load project"), " stores the workflow state in a .sangerproject file so completed processing and retrieved BLAST results do not need to be repeated."),
                    p("The final identification is decision-support. It does not replace locus-specific taxonomy, type/reference inspection, morphology or additional loci when those are required.")
                  )
                )
              ),
              div(class="help-card",
                h3("Stage 2 project architecture"),
                p("PITAX stores Project → Isolate → Locus → Read as explicit linked objects. One isolate may have several loci, and one locus may have a single read or separate Forward and Reverse reads."),
                p("The upload barcode is an immutable technical source ID and is never parsed as biological identity. Isolate, Gene/Locus and Direction are edited as separate fields in Rename before trimming; PITAX generates the final <Isolate>_<Locus>_<F/R> label from those fields. Duplicate AB1 basenames are blocked because they cannot be represented safely by the established source-read key."),
                div(class="about-callout",
                  strong("Stage boundary: "),
                  "Stage 2 stores identity and pairing only. Stage 3 builds a separate isolate-level sequence with overlap evidence and never rewrites the source reads."
                )
              ),
              div(class="help-card",
                h3("Stage 3 isolate-level sequence"),
                p("A sequencing run must contain one Gene/Locus. PITAX reverse-complements Reverse reads in a derived view, aligns each explicit pair and records every call with source positions and available basecaller quality."),
                p("Single reads remain valid representatives. Weak pair overlaps and unresolved contradictions block downstream analysis; they are not hidden by concatenation or majority voting."),
                div(class="about-callout", strong("Accepted validation boundary: "), "The linked conflict-review mechanics, revisions and controlled AB1 fixtures are accepted for continued development. An independently sequenced paired-AB1 set remains a deferred biological validation item.")
              ),
              div(class="help-card",
                h3("Stage 4 multi-locus profile"),
                p("Each Gene/Locus is processed in its own PITAX project. Stage 4 imports those completed projects and joins their evidence only when the explicit Isolate code matches."),
                p("Every locus retains its sequence, consensus revision, BLAST RID, taxonomic result, limitation and reference context. Two concordant loci may support the same rank; a conflicting locus is never removed by a numerical majority."),
                div(class="about-callout", strong("Current alpha boundary: "), "alpha.9.2 retains the evidence-preserving profile, visual multi-isolate review, provenance, duplicate Isolate/Locus blocking, concordance/conflict states and exports. Taxon-specific marker recommendations still require a separately reviewed literature layer.")
              )
            ),

            tabPanel("Trimming & QC",
              div(class="help-card",
                h3("Trimming terminology"),
                div(class="method-step",
                  div(class="method-num", "1"),
                  div(h4("Peak ratio"), p("Ratio between the signal of the called base and the second-highest channel at that base call. Larger values indicate a clearer dominant peak."))
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(h4("Good start"), p("The first sustained window that satisfies the configured signal/peak-ratio rule. Bases before this point are excluded from the processed sequence."))
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(h4("Collapse"), p("The first sustained region where relative signal or peak-ratio behavior indicates deterioration. The processed read is normally ended before this point."))
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(h4("Minimum usable length"), p("A workflow threshold used to flag a processed sequence as SHORT_AFTER_TRIMMING rather than OK."))
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(h4("Ambiguous peak review"),
                      p("The app flags individual called-base positions where another A/C/G/T dye channel competes strongly with the current call, where the current call is not locally dominant, or where the call is ambiguous."),
                      p("The comparison uses a narrow trace window around each called peak and suppresses weak-signal competition using a run-specific signal baseline. By default, a called/competitor peak ratio ≤1.25 is Strong and ≤1.75 is Moderate; the local signal must be at least 20% of the retained-region median called signal. These are application review heuristics, not biological diagnoses."),
                      p("Flags are review targets. Left-click centers the chromatogram; right-click opens manual curation actions. No sequence-changing action is applied without an explicit confirmation step."))
                ),
                div(class="method-step",
                  div(class="method-num", "6"),
                  div(h4("Manual curation and provenance"),
                      p("The raw AB1 trace and original automatic base calls are never overwritten. The app maintains a curated sequence layer on top of the automatic trim. A user can change a base, set a two-base IUPAC ambiguity code, trim from either side through a flagged position, or mark the call as reviewed and unchanged."),
                      p("Every confirmed action is written to the Manual Curation audit log with sample, raw base position, before/after values, method, chromatogram evidence, timestamp, transaction ID and revision. Undo and redo are logged as well. Automatic trim boundaries are retained separately so the curated result can always be compared with the automatic result."),
                      p("For a flagged position inside the retained read, the curation menu shows both left- and right-trim actions. A directional recommendation is displayed only when the flag is near a current trim edge: within 15% of the retained length, bounded to 8–40 bases. The recommendation never performs a trim by itself and still requires explicit confirmation."))
                ),
                div(class="method-step",
                  div(class="method-num", "7"),
                  div(h4("High-confidence bulk correction"),
                      p("The bulk correction tool always shows a preview before applying changes. By default, a position is proposed only when the current called base is not locally dominant, the alternative channel is the strongest local channel, alternative/current signal is at least 1.80, alternative/third-channel signal is at least 2.00, the alternative peak maximum lies within ±2 trace samples of the called-base peak position, and the alternative signal is at least 50% of the retained-region median called signal."),
                      p("These thresholds are application heuristics and can be edited from QC & Chromatogram using the Criteria button next to Auto-correct. The active values are saved in the project/run settings and therefore remain part of the analysis provenance. The complete batch is applied as one undoable transaction after user confirmation. Ambiguous double peaks are not bulk-corrected merely because two channels are similar."))
                )
              ),
              div(class="help-card",
                h3("QC evidence"),
                p(strong("QC"), " means Quality Control. Called-base signal and peak-ratio plots are retained in checkpoint/final ZIP exports."),
                p("The ambiguous-peak table links directly to the interactive chromatogram. Right-clicking a row opens the curation actions, while the History view exposes the complete per-sample audit trail."),
                p("If a curated sequence is changed after BLAST has already been run, the previous BLAST/taxonomic result for that sample is marked stale and removed from active interpretation until BLAST is rerun on the new curated sequence.")
              )
            ),

            tabPanel("NCBI BLAST",
              div(class="help-card",
                h3("How BLAST evidence is represented"),
                div(class="method-step",
                  div(class="method-num", "1"),
                  div(h4("RID — Request ID"), p("Identifier returned by NCBI for an individual BLAST job. The app keeps Sample ↔ RID ↔ retrieved hits linked."))
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(h4("Identity"), p("Percentage of identical positions in the representative alignment for the accession."))
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(h4("Query coverage"), p("Percentage of the processed query covered by the subject hit. Multiple local segments from the same accession are combined rather than counted as separate database hits."))
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(h4("E-value and Bit score"), p("E-value describes the expected chance occurrence of a match of similar quality; Bit score is a normalized alignment score used here to compare the leading taxon with competing hits."))
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(h4("HSP — High-scoring Segment Pair"), p("A local alignment segment within a BLAST hit. Multiple HSPs for the same accession are aggregated so that one NCBI record contributes one accession-level hit."))
                )
              ),
              div(class="about-callout",
                strong("Why accession-level aggregation matters: "),
                "without aggregation, one NCBI record containing several HSPs could receive several votes and distort both query-coverage statistics and taxonomic consensus."
              )
            ),

            tabPanel("Taxonomic algorithm",
              div(class="help-card",
                h3("Evidence-first decision logic"),
                p(class="about-lead",
                  "The v2.14 engine uses one decision tree. It starts from the closest molecular match, checks whether other named taxa are practically indistinguishable, and then reports the lowest taxonomic rank that remains defensible. Database abundance is kept as context rather than used as a vote."),

                div(class="method-step",
                  div(class="method-num", "1"),
                  div(
                    h4("One accession = one hit", span(class="evidence-badge evidence-published", "BLAST structure")),
                    p("Multiple HSPs from the same NCBI accession are aggregated before taxonomic interpretation. A record therefore contributes one accession-level piece of evidence regardless of how many local alignment segments it contains.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(
                    h4("Best molecular match", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("Partial 100% matches are not allowed to outrank near-full-length evidence merely because their aligned segment is short. The app first prefers the >=90% coverage tier (or >=80% if needed), then keeps hits within 2 percentage points of the best query coverage in that tier. Identity is ranked first inside that near-best-coverage band; coverage, Bit score and E-value are tie-breakers.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(
                    h4("Close alternatives", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("A different named taxon is considered a close alternative when its best comparable accession is within 0.5 identity percentage points of the best match and is no more than 2 query-coverage points below it. These narrow windows are review heuristics, not universal species-delimitation thresholds.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(
                    h4("Conservative taxonomic fallback", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("If no close species alternative remains and sequence/reference evidence is strong, a species recommendation can be made. If several close species remain but they all belong to the same genus, species is reported as unresolved while the genus can remain a high-confidence recommendation. If close alternatives extend across genera, genus-level confidence is withheld and an LCA may be used as a fallback.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(
                    h4("Species evidence profile", span(class="evidence-badge evidence-heuristic", "Audit view")),
                    p("For every resolved species the app reports its strongest comparable accession, Identity, query coverage, Bit score, reference quality and the number of unique accessions returned for that species. The accession count describes database representation only; 23 records do not automatically defeat a better molecular match represented by one or two records.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "6"),
                  div(
                    h4("Sequence evidence and locus discrimination", span(class="evidence-badge evidence-published", "Biological context")),
                    p("Sequence quality and taxonomic discrimination are separate questions. A read can have excellent Identity and coverage yet fail to distinguish several species. In that situation the sequence evidence can remain High while the locus is flagged as having poor species-level discrimination. For ITS, broad literature benchmarks remain context for sequence-evidence strength rather than hard species cutoffs.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "7"),
                  div(
                    h4("Reference context", span(class="evidence-badge evidence-published", "Reference quality")),
                    p("Type-material wording, RefSeq status, ordinary GenBank records and unresolved/environmental annotations are retained separately. Reference quality can change confidence, but it does not replace the direct Identity/coverage comparison.")
                  )
                )
              ),
              div(class="about-callout",
                strong("Interpretation principle: "),
                "Best molecular match first; close alternatives determine how far the identification can safely be taken; database frequency is supporting context only. A poor species-level result can therefore be a property of the locus rather than a failure of the sequence."
              )
            ),

            tabPanel("Scientific references",
              p(class="about-lead",
                "The links below are the main scientific and NCBI sources used to justify the biological context of the workflow. External links open in a new browser tab."),

              div(class="reference-card",
                h4("Schoch et al. 2012 — ITS as the primary fungal DNA barcode"),
                p("Compared candidate fungal barcode regions and proposed the nuclear ribosomal ITS region as the universal DNA barcode marker for Fungi."),
                tags$a(href="https://www.pnas.org/doi/10.1073/pnas.1117018109", target="_blank", rel="noopener noreferrer", "Open PNAS article ↗")
              ),

              div(class="reference-card",
                h4("Vu et al. 2018 — large-scale fungal barcode thresholds"),
                p("Large-scale analysis of filamentous fungal DNA barcodes. The broad ITS benchmarks used as context in the app (approximately 99.6% species and 94.3% genus) come from this work."),
                tags$a(href="https://pmc.ncbi.nlm.nih.gov/articles/PMC6020082/", target="_blank", rel="noopener noreferrer", "Open full text at NCBI/PMC ↗")
              ),

              div(class="reference-card",
                h4("Garnica et al. 2016 — lineage-dependent ITS thresholds"),
                p("Shows that optimal ITS similarity thresholds vary among lineages and that no single identity cutoff is universally reliable even within a diverse fungal genus."),
                tags$a(href="https://academic.oup.com/femsec/article/92/4/fiw045/2197947", target="_blank", rel="noopener noreferrer", "Open FEMS Microbiology Ecology article ↗")
              ),

              div(class="reference-card",
                h4("Stielow et al. 2015 — secondary fungal DNA barcodes"),
                p("Evaluated secondary fungal barcode loci and universal primers, supporting the use of additional loci when ITS alone does not provide sufficient resolution."),
                tags$a(href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4713107/", target="_blank", rel="noopener noreferrer", "Open full text at NCBI/PMC ↗")
              ),

              div(class="reference-card",
                h4("NCBI RefSeq Targeted Loci — fungal ITS"),
                p("Describes NCBI's curated fungal ITS reference collection. Sequences are mostly derived from type material and are maintained with specimen/taxonomic context."),
                tags$a(href="https://www.ncbi.nlm.nih.gov/refseq/targetedloci/", target="_blank", rel="noopener noreferrer", "Open NCBI RefSeq Targeted Loci ↗")
              ),

              div(class="reference-card",
                h4("NCBI BLAST documentation"),
                p("Technical definitions for BLAST output fields including Bit score, HSPs and query-coverage measures used by the application parser."),
                tags$a(href="https://www.ncbi.nlm.nih.gov/books/NBK279684/", target="_blank", rel="noopener noreferrer", "Open NCBI BLAST+ manual ↗")
              ),

              div(class="reference-card",
                h4("UNITE — fungal ITS reference and Species Hypotheses"),
                p("A fungal ITS-centered identification resource and Species Hypotheses system. It is included here as scientific context and a possible future curated-reference extension; the current app does not query UNITE directly."),
                tags$a(href="https://unite.ut.ee/", target="_blank", rel="noopener noreferrer", "Open UNITE ↗")
              ),

              div(class="about-callout",
                strong("Important distinction: "),
                "the near-best coverage band and the close-match windows (0.5 Identity percentage points and 2 coverage points) are application-specific review heuristics, not published species-delimitation thresholds. Published sources guide the biological context, while these rules should be calibrated against known isolates."
              )
            ),

            tabPanel("Version history",
              div(class="help-card",
                h3("CHANGELOG"),
                p("Every version records changes to processing, BLAST parsing, taxonomic interpretation and reporting so exported results can be traced back to the software logic used at the time."),
                verbatimTextOutput("changelog_text")
              )
            )
          )
        )
      )
    )
  )
)
