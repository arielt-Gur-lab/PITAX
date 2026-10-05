# ============================================================
# Shared export component catalog
# Explicit definitions only; table presence in state is not enough.
# ============================================================

pitax_export_col <- function(id, label, type = c("text", "number", "logical", "datetime"),
                             link_key = FALSE, display_default = TRUE) {
  type <- match.arg(type)
  list(
    id = as.character(id)[1],
    label = as.character(label)[1],
    type = type,
    link_key = isTRUE(link_key),
    display_default = isTRUE(display_default)
  )
}

pitax_export_component <- function(id, label, stages, kind = c("table", "sequence", "plot", "raw_file", "text"),
                                   approval = c("approved_result", "evidence_audit", "historical"),
                                   row_meaning = "",
                                   columns = list(),
                                   link_keys = character(),
                                   split_by = character(),
                                   sequence_role = NULL,
                                   plot_type = NULL,
                                   default_selected = FALSE) {
  kind <- match.arg(kind)
  approval <- match.arg(approval)
  list(
    id = as.character(id)[1],
    label = as.character(label)[1],
    stages = unique(as.character(stages)),
    kind = kind,
    approval = approval,
    row_meaning = as.character(row_meaning)[1],
    columns = columns,
    link_keys = unique(as.character(link_keys)),
    split_by = unique(as.character(split_by)),
    sequence_role = sequence_role,
    plot_type = plot_type,
    default_selected = isTRUE(default_selected)
  )
}

pitax_export_catalog <- function() {
  list(
    pitax_export_component(
      "upload_file_inventory", "Upload file inventory", "upload", "table", "evidence_audit",
      row_meaning = "One source AB1 file known to the session.",
      columns = list(
        pitax_export_col("Source_ID", "Source ID", "text", TRUE),
        pitax_export_col("File", "File name", "text", TRUE),
        pitax_export_col("Size_KB", "Size (KB)", "number"),
        pitax_export_col("Bytes_available", "Original bytes available", "logical")
      ),
      link_keys = c("Source_ID", "File"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "upload_raw_ab1", "Original AB1 files", "upload", "raw_file", "evidence_audit",
      row_meaning = "Original sequencer bytes when still present on disk.",
      link_keys = c("Source_ID", "File")
    ),
    pitax_export_component(
      "assay_profiles_draft", "Assay profiles (draft)", "settings", "table", "evidence_audit",
      row_meaning = "Current assay profile definitions in the editor, which may differ from settings used during trimming.",
      columns = list(
        pitax_export_col("Assay_ID", "Assay ID", "text", TRUE),
        pitax_export_col("Assay_Name", "Assay name", "text"),
        pitax_export_col("Locus_ID", "Locus ID", "text", TRUE),
        pitax_export_col("Locus_Display_Name", "Locus display name", "text"),
        pitax_export_col("Forward_Primer_Name", "Forward primer name", "text"),
        pitax_export_col("Forward_Primer_Sequence", "Forward primer sequence", "text"),
        pitax_export_col("Reverse_Primer_Name", "Reverse primer name", "text"),
        pitax_export_col("Reverse_Primer_Sequence", "Reverse primer sequence", "text"),
        pitax_export_col("Expected_Amplicon_Length", "Expected amplicon length", "number"),
        pitax_export_col("Maximum_Sequence_Position", "Maximum sequence position", "number")
      ),
      link_keys = c("Assay_ID", "Locus_ID"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "project_defaults_draft", "Project trimming defaults (draft)", "settings", "table", "evidence_audit",
      row_meaning = "Shared project trimming defaults currently shown in Assay setup.",
      columns = list(
        pitax_export_col("Setting", "Setting", "text", TRUE),
        pitax_export_col("Value", "Value", "text")
      ),
      link_keys = "Setting",
      default_selected = TRUE
    ),
    pitax_export_component(
      "processing_settings_applied", "Processing settings used during trimming", "settings", "table", "evidence_audit",
      row_meaning = "Per-read settings snapshot stored on each processed result.",
      columns = list(
        pitax_export_col("Source_ID", "Source ID", "text", TRUE),
        pitax_export_col("Setting", "Setting", "text", TRUE),
        pitax_export_col("Value", "Value", "text")
      ),
      link_keys = c("Source_ID", "Setting"),
      split_by = c("sample", "isolate", "locus")
    ),
    pitax_export_component(
      "assign_identity_table", "Final read / FASTA names and biological identity", "rename", "table", "approved_result",
      row_meaning = "One assigned read identity row.",
      columns = list(
        pitax_export_col("Source_ID", "Source ID", "text", TRUE),
        pitax_export_col("File", "Source file", "text"),
        pitax_export_col("Isolate", "Isolate", "text", TRUE),
        pitax_export_col("Assay_ID", "Assay ID", "text", TRUE),
        pitax_export_col("Locus", "Locus", "text", TRUE),
        pitax_export_col("Direction", "Direction", "text"),
        pitax_export_col("Primer", "Primer", "text"),
        pitax_export_col("Final_Name", "Final name", "text", TRUE),
        pitax_export_col("Inference", "Inference", "text", display_default = FALSE),
        pitax_export_col("Notes", "Notes", "text", display_default = FALSE)
      ),
      link_keys = c("Source_ID", "Isolate", "Assay_ID", "Locus", "Final_Name"),
      split_by = c("isolate", "locus", "assay"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "assign_rename_map", "Rename map", "rename", "table", "evidence_audit",
      row_meaning = "Source ID to final FASTA name mapping.",
      columns = list(
        pitax_export_col("Original_name", "Original name", "text", TRUE),
        pitax_export_col("New_name", "New name", "text", TRUE)
      ),
      link_keys = c("Original_name", "New_name")
    ),
    pitax_export_component(
      "assign_architecture", "Project architecture tables", "rename", "table", "evidence_audit",
      row_meaning = "Project, isolate, locus and read architecture rows.",
      columns = list(
        pitax_export_col("Table", "Architecture table", "text", TRUE),
        pitax_export_col("Key", "Row key", "text", TRUE),
        pitax_export_col("Payload", "Serialized row", "text")
      ),
      link_keys = c("Table", "Key")
    ),
    pitax_export_component(
      "qc_summary", "Trim and QC summary", "qc", "table", "approved_result",
      row_meaning = "One processed source-read QC summary.",
      columns = list(
        pitax_export_col("sample_id", "Source ID", "text", TRUE),
        pitax_export_col("final_name", "Final name", "text", TRUE),
        pitax_export_col("target", "Locus / target", "text"),
        pitax_export_col("raw_length", "Raw length", "number"),
        pitax_export_col("trimmed_length", "Trimmed length", "number"),
        pitax_export_col("trim_start", "Trim start", "number"),
        pitax_export_col("trim_end", "Trim end", "number"),
        pitax_export_col("status", "Status", "text"),
        pitax_export_col("reason", "Reason", "text"),
        pitax_export_col("manual_curation", "Manual curation", "text"),
        pitax_export_col("curation_revision", "Curation revision", "number")
      ),
      link_keys = c("sample_id", "final_name"),
      split_by = c("sample", "isolate", "locus"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "qc_ab1_run_evidence", "AB1 run evidence", "qc", "table", "evidence_audit",
      row_meaning = "Run-level AB1 evidence summary for one source read.",
      link_keys = c("Sample", "Final_Name"),
      split_by = c("sample")
    ),
    pitax_export_component(
      "qc_ab1_base_evidence", "AB1 per-base evidence", "qc", "table", "evidence_audit",
      row_meaning = "One base-call evidence row from the AB1 audit detail.",
      link_keys = c("Sample_ID", "Position"),
      split_by = c("sample")
    ),
    pitax_export_component(
      "qc_peak_flags", "Ambiguous-peak flags", "qc", "table", "evidence_audit",
      row_meaning = "One flagged base position in the retained sequence.",
      link_keys = c("Sample", "Position"),
      split_by = c("sample")
    ),
    pitax_export_component(
      "qc_active_curation", "Active QC curation changes", "qc", "table", "evidence_audit",
      row_meaning = "Current applied manual curation actions retained on the active sequence.",
      link_keys = c("Sample", "Revision", "Position"),
      split_by = c("sample")
    ),
    pitax_export_component(
      "qc_historical_audit", "Historical QC curation audit log", "qc", "table", "historical",
      row_meaning = "Historical audit events, including undo and superseded actions, marked as historical.",
      link_keys = c("Sample", "Transaction_ID", "Revision"),
      split_by = c("sample")
    ),
    pitax_export_component(
      "qc_trimmed_sequences", "Trimmed and curated source sequences", "qc", "sequence", "approved_result",
      row_meaning = "One trimmed/curated source-read sequence.",
      sequence_role = "trimmed_curated_read",
      link_keys = c("sample_id", "final_name"),
      split_by = c("sample", "isolate", "locus"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "qc_metrics_plots", "QC metrics plots", "qc", "plot", "evidence_audit",
      row_meaning = "Called-base signal and peak-ratio metrics plot for one source read.",
      plot_type = "qc_metrics",
      link_keys = "sample_id",
      split_by = "sample"
    ),
    pitax_export_component(
      "qc_chromatograms", "Chromatogram PNGs", "qc", "plot", "evidence_audit",
      row_meaning = "Full-window chromatogram rendering for one source read.",
      plot_type = "chromatogram",
      link_keys = "sample_id",
      split_by = "sample"
    ),
    pitax_export_component(
      "consensus_summary", "Analysis sequence summary", c("consensus", "export"), "table", "approved_result",
      row_meaning = "One isolate-level analysis sequence summary row.",
      columns = list(
        pitax_export_col("Consensus_ID", "Consensus ID", "text", TRUE),
        pitax_export_col("Final_Name", "Final name", "text", TRUE),
        pitax_export_col("Isolate", "Isolate", "text", TRUE),
        pitax_export_col("Locus", "Locus", "text", TRUE),
        pitax_export_col("Status", "Status", "text"),
        pitax_export_col("Length", "Length", "number"),
        pitax_export_col("Revision", "Revision", "number"),
        pitax_export_col("Overlap", "Overlap", "number"),
        pitax_export_col("Identity_percent", "Identity percent", "number"),
        pitax_export_col("Review_positions", "Review positions", "number")
      ),
      link_keys = c("Consensus_ID", "Isolate", "Locus", "Final_Name"),
      split_by = c("isolate", "locus"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "consensus_analysis_sequences", "Analysis-oriented sequences", c("consensus", "export"), "sequence", "approved_result",
      row_meaning = "One analysis sequence used for BLAST (independent read or consensus).",
      sequence_role = "analysis_oriented",
      link_keys = c("Consensus_ID", "Final_Name"),
      split_by = c("isolate", "locus", "sample"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "consensus_column_evidence", "Consensus per-column evidence", "consensus", "table", "evidence_audit",
      row_meaning = "One alignment column decision for a consensus record.",
      link_keys = c("Consensus_ID", "Alignment_Column"),
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "consensus_pairwise_alignments", "Pairwise alignments", "consensus", "text", "evidence_audit",
      row_meaning = "Forward/Reverse overlap alignment text for one consensus.",
      link_keys = "Consensus_ID",
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "consensus_conflicts", "Consensus conflicts", "consensus", "table", "evidence_audit",
      row_meaning = "One unresolved or reviewed conflict column.",
      link_keys = c("Consensus_ID", "Alignment_Column"),
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "consensus_active_curation", "Active consensus curation", "consensus", "table", "evidence_audit",
      row_meaning = "Active consensus review decisions retained on the current revision.",
      link_keys = c("Consensus_ID", "Revision"),
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "consensus_historical_audit", "Historical consensus audit log", "consensus", "table", "historical",
      row_meaning = "Historical consensus audit events marked as historical.",
      link_keys = c("Consensus_ID", "Revision", "Timestamp"),
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "blast_jobs", "BLAST jobs and status", "blast", "table", "evidence_audit",
      row_meaning = "One NCBI BLAST job bound to an analysis sequence revision.",
      columns = list(
        pitax_export_col("final_name", "Final name", "text", TRUE),
        pitax_export_col("original_name", "Analysis ID", "text", TRUE),
        pitax_export_col("rid", "RID", "text", TRUE),
        pitax_export_col("status", "Status", "text"),
        pitax_export_col("consensus_revision", "Consensus revision", "number"),
        pitax_export_col("database", "Database", "text"),
        pitax_export_col("submitted_at", "Submitted at", "datetime"),
        pitax_export_col("last_checked_at", "Last checked at", "datetime")
      ),
      link_keys = c("original_name", "rid", "final_name"),
      split_by = c("sample", "isolate"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "blast_hits", "BLAST hits by accession", "blast", "table", "evidence_audit",
      row_meaning = "One accession-level BLAST hit. Top hit is not an approved taxonomic identity.",
      columns = list(
        pitax_export_col("final_name", "Final name", "text", TRUE),
        pitax_export_col("original_name", "Analysis ID", "text", TRUE),
        pitax_export_col("rid", "RID", "text", TRUE),
        pitax_export_col("rank", "Rank", "number"),
        pitax_export_col("accession", "Accession", "text", TRUE),
        pitax_export_col("organism", "Organism label from BLAST", "text"),
        pitax_export_col("identity_percent", "Identity percent", "number"),
        pitax_export_col("query_coverage_percent", "Query coverage percent", "number"),
        pitax_export_col("evalue", "E-value", "number"),
        pitax_export_col("bit_score", "Bit score", "number"),
        pitax_export_col("match_support", "Match support", "text")
      ),
      link_keys = c("original_name", "rid", "accession"),
      split_by = c("sample", "isolate")
    ),
    pitax_export_component(
      "blast_raw", "Stored BLAST raw results", "blast", "text", "evidence_audit",
      row_meaning = "Raw NCBI result text retained for one RID.",
      link_keys = "rid",
      split_by = "sample"
    ),
    pitax_export_component(
      "blast_query_sequences", "BLAST query sequences", "blast", "sequence", "evidence_audit",
      row_meaning = "Analysis sequence used as the BLAST query.",
      sequence_role = "blast_query",
      link_keys = c("original_name", "final_name"),
      split_by = c("sample", "isolate", "locus")
    ),
    pitax_export_component(
      "taxonomy_team_summary", "Team taxonomic summary", "taxonomy", "table", "approved_result",
      row_meaning = "One compact team-facing identification summary row for an analysis sequence.",
      link_keys = c("Sample", "Original_sample", "RID"),
      split_by = c("sample", "isolate"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "taxonomy_summary", "Detailed taxonomic summaries", "taxonomy", "table", "approved_result",
      row_meaning = "One computed taxonomic interpretation for an analysis sequence.",
      link_keys = c("final_name", "original_name", "rid"),
      split_by = c("sample", "isolate")
    ),
    pitax_export_component(
      "taxonomy_hits", "Taxonomy-enriched hits", "taxonomy", "table", "evidence_audit",
      row_meaning = "One enriched accession hit used by taxonomic interpretation.",
      link_keys = c("original_name", "rid", "accession"),
      split_by = c("sample", "isolate")
    ),
    pitax_export_component(
      "taxonomy_species_evidence", "Species evidence counts", "taxonomy", "table", "evidence_audit",
      row_meaning = "One species-level evidence aggregate.",
      link_keys = c("original_name", "rid", "taxon"),
      split_by = c("sample", "isolate")
    ),
    pitax_export_component(
      "taxonomy_score_plot", "Taxonomy score landscape plots", "taxonomy", "plot", "evidence_audit",
      row_meaning = "Bit-score landscape for one analyzed sequence.",
      plot_type = "taxonomy_scores",
      link_keys = "original_name",
      split_by = "sample"
    ),
    pitax_export_component(
      "multilocus_profiles", "Multi-locus isolate profiles", "multilocus", "table", "approved_result",
      row_meaning = "One isolate-level multi-locus profile.",
      link_keys = "Isolate",
      split_by = "isolate",
      default_selected = TRUE
    ),
    pitax_export_component(
      "multilocus_evidence", "Multi-locus per-locus evidence", "multilocus", "table", "evidence_audit",
      row_meaning = "One Isolate + Locus evidence row. Loci remain separate.",
      link_keys = c("Isolate", "Locus", "Consensus_ID"),
      split_by = c("isolate", "locus"),
      default_selected = TRUE
    ),
    pitax_export_component(
      "multilocus_sources", "Multi-locus source projects", "multilocus", "table", "evidence_audit",
      row_meaning = "One source project included in the multi-locus profile.",
      link_keys = c("Source", "Source_MD5")
    ),
    pitax_export_component(
      "multilocus_sequences", "Multi-locus sequences", "multilocus", "sequence", "approved_result",
      row_meaning = "One locus sequence belonging to an isolate profile. Loci are never concatenated.",
      sequence_role = "multilocus_locus",
      link_keys = c("Isolate", "Locus"),
      split_by = c("isolate", "locus")
    ),
    pitax_export_component(
      "multilocus_plot", "Multi-locus evidence plots", "multilocus", "plot", "evidence_audit",
      row_meaning = "Identity and query-coverage plot for one isolate.",
      plot_type = "multilocus_evidence",
      link_keys = "Isolate",
      split_by = "isolate"
    )
  )
}

pitax_export_component_ids <- function() {
  vapply(pitax_export_catalog(), function(x) x$id, character(1))
}

pitax_export_get_component <- function(id) {
  id <- as.character(id)[1]
  for (comp in pitax_export_catalog()) {
    if (identical(comp$id, id)) return(comp)
  }
  NULL
}

pitax_export_components_for_stage <- function(stage_id) {
  stage_id <- as.character(stage_id)[1]
  Filter(function(comp) stage_id %in% comp$stages || identical(stage_id, "export"), pitax_export_catalog())
}

pitax_export_default_component_ids <- function(stage_id, scope = c("stage", "project")) {
  scope <- match.arg(scope)
  comps <- if (identical(scope, "project") || identical(stage_id, "export")) {
    pitax_export_catalog()
  } else {
    Filter(function(comp) stage_id %in% comp$stages, pitax_export_catalog())
  }
  ids <- character()
  for (comp in comps) {
    if (isTRUE(comp$default_selected)) ids <- c(ids, comp$id)
  }
  unique(ids)
}
