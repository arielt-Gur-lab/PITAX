  # ---------------- Export records ----------------
  read_export_records <- reactive({
    req(rv$results,rv$rename)
    records <- list()
    for(original_name in names(rv$results)) {
      r <- rv$results[[original_name]]
      idx <- match(original_name,rv$rename$Original_name)
      r$final_name <- if(!is.na(idx)) rv$rename$New_name[idx] else original_name
      records[[original_name]] <- r
    }
    records
  })
  read_export_summary_df <- reactive({
    req(rv$summary,rv$rename)
    df <- rv$summary
    df$final_name <- rv$rename$New_name[match(df$sample_id,rv$rename$Original_name)]
    df
  })
  export_records <- reactive({
    req(is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
    stage3_analysis_records(rv$consensus_set)
  })
  export_summary_df <- reactive({
    req(is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
    stage3_analysis_summary(rv$consensus_set)
  })

  observeEvent(input$to_export, {
    e <- stage3_consensus_gate_error(rv$consensus_set, rv$results)
    if(!is.null(e)){showNotification(e,type="error",duration=10);return()}
    updateTabsetPanel(session,"pipeline_step",selected="export")
  })

  project_export_stem <- function() {
    if (is.list(rv$architecture) && is.data.frame(rv$architecture$loci) && "Locus" %in% names(rv$architecture$loci)) {
      locus_names <- unique(trimws(as.character(rv$architecture$loci$Locus)))
      locus_names <- locus_names[nzchar(locus_names)]
      if (length(locus_names) > 1L) return("PITAX_multi_locus")
      if (length(locus_names) == 1L) return(clean_fasta_name(locus_names))
    }
    if (!is.null(rv$settings)) return(clean_fasta_name(stage2_scalar_text(rv$settings$target, "PITAX_project")))
    "PITAX_project"
  }

  output$export_summary <- renderUI({
    req(rv$summary,rv$settings, is.null(stage3_consensus_gate_error(rv$consensus_set, rv$results)))
    consensus_summary <- rv$consensus_set$summary
    loci <- unique(as.character(consensus_summary$Locus))
    architecture_summary <- stage2_architecture_summary(rv$architecture)
    tagList(
      p(strong("Gene / locus: "), paste(loci, collapse = ", ")),
      p(strong("Source architecture: "), paste0(architecture_summary$Isolates, " isolate(s), ", architecture_summary$Reads, " read(s)")),
      p(strong("Forward primer: "),ifelse(rv$settings$forward_primer=="","Not specified",rv$settings$forward_primer)),
      p(strong("Reverse primer: "),ifelse(rv$settings$reverse_primer=="","Not specified",rv$settings$reverse_primer)),
      p(strong("Source reads processed: "),nrow(rv$summary)),
      p(strong("Isolate-level sequences for export: "),sum(consensus_summary$Length > 0, na.rm = TRUE))
    )
  })

  # ---------------- Checkpoints ----------------
  output$download_trim_checkpoint <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_checkpoint_B_qc.zip"),
    content=function(file){
      write_checkpoint_zip(file,"qc",read_export_records(),read_export_summary_df(),rv$settings,rv$rename,results=rv$results,
                           read_assignments=rv$read_assignments,architecture=rv$architecture)
    })
  output$download_rename_checkpoint <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_checkpoint_A_rename.zip"),
    content=function(file) {
      req(is.null(stage2_identity_error(rv$read_assignments)), is.list(rv$architecture), is.list(rv$settings))
      write_assignment_checkpoint_zip(file, rv$rename, rv$read_assignments, rv$architecture, rv$settings)
    })

  output$download_blast_fasta <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_BLAST.fasta"),
    content=function(file) writeLines(make_fasta(export_records(),FALSE),file))
  output$download_full_fasta <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_isolate_level_sequences.fasta"),
    content=function(file) writeLines(make_fasta(export_records(),TRUE,export_summary_df()),file))
  output$download_summary_csv <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_isolate_level_summary.csv"),
    content=function(file) write.csv(export_summary_df(),file,row.names=FALSE,fileEncoding="UTF-8"))
  output$download_all_zip <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_Sanger_pipeline_results.zip"),
    content=function(file) write_checkpoint_zip(file,"final",export_records(),export_summary_df(),rv$settings,rv$rename,results=rv$results,
                                                read_assignments=rv$read_assignments,architecture=rv$architecture,consensus_set=rv$consensus_set))

