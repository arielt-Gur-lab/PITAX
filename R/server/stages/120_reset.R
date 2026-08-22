  # ---------------- Reset ----------------
  reset_pipeline_state <- function() {
    rv$results <- list(); rv$summary <- NULL; rv$rename <- NULL; rv$settings <- NULL
    rv$read_assignments <- stage2_empty_assignments(); rv$architecture <- NULL; rv$consensus_set <- stage3_empty_consensus_set(); rv$multilocus_profile <- stage4_empty_profile(); rv$project_migration_log <- ""
    rv$blast_jobs <- rv$blast_jobs[0,]; rv$blast_raw <- list(); rv$blast_ids <- data.frame(); rv$blast_hits <- data.frame()
    rv$ncbi_last_contact <- as.POSIXct(NA); rv$blast_batch_status_text <- "No batch operation has been run yet."
    rv$taxonomy_summary <- data.frame(); rv$taxonomy_hits <- data.frame(); rv$taxonomy_counts <- data.frame()
    rv$taxonomy_status_text <- "No taxonomic analysis has been run yet."
    rv$taxonomy_batch_status_text <- "No batch taxonomic analysis has been run yet."
    rv$project_status_text <- "Current session has not been saved as a project."
    rv$project_loaded_name <- ""
    updateTabsetPanel(session,"pipeline_step",selected="upload")
  }
  observeEvent(input$reset_pipeline, reset_pipeline_state())
  observeEvent(input$reset_pipeline_tax, reset_pipeline_state())
