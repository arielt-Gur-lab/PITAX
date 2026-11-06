# Test the ACTUAL PITAX FASTA header Shiny stage observers in a small app.
# This catches browser-side reactive synchronization bugs not visible in model tests.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
source("R/export/availability.R", local = FALSE)
source("R/export/header_templates.R", local = FALSE)
# shinyApp(ui, server) has no appDir, so ./www is NOT served automatically.
# Explicitly serve the same production JavaScript through a test-only prefix.
stopifnot(file.exists("www/pitax-header-profiles.js"))
shiny::addResourcePath("pitaxassets", normalizePath("www"))

e2e_snap <- list(
  read_assignments = data.frame(
    Source_ID = c("S1", "S2"),
    Final_Name = c("ISO1_TEF_R", "ISO2_TEF_R"),
    Isolate = c("ISO1", "ISO2"), Locus = c("TEF", "TEF"),
    Direction = c("Reverse", "Reverse"), stringsAsFactors = FALSE
  ),
  results = list(), consensus_set = list(records = list()),
  multilocus_profile = list()
)

ui <- fluidPage(
  tags$head(tags$script(src = "pitaxassets/pitax-header-profiles.js")),
  titlePanel("PITAX metadata editor integration"),
  textAreaInput("export_header_template", "Header template",
                value = pitax_header_default_template(), width = "100%"),
  textInput("export_header_template_name", "Template name", value = "Research metadata"),
  selectInput("export_header_saved_template", "Saved templates", choices = "Research metadata"),
  actionButton("export_header_save_template", "Save / update template"),
  actionButton("export_header_load_template", "Load selected template"),
  textInput("export_header_new_field", "New metadata field", value = ""),
  actionButton("export_header_add_field", "Create and use field"),
  textOutput("export_header_registered_fields"),
  selectInput("export_header_existing_field", "All registered metadata fields",
              choices = c("Organism", "Experiment", "Seq_Date")),
  actionButton("export_header_activate_field", "Use selected field in header"),
  actionButton("export_header_remove_from_header", "Remove from active header"),
  actionButton("export_header_delete_field", "Delete field permanently"),
  selectInput("export_header_bulk_field", "Field to set", choices = "Organism"),
  textInput("export_header_bulk_value", "Value", value = ""),
  selectInput("export_header_bulk_scope", "Apply to",
              choices = c("All records" = "all", "Selected" = "selected"), selected = "all"),
  actionButton("export_header_apply_bulk", "Apply metadata value"),
  textOutput("export_header_active_template"),
  textOutput("export_header_persistence_status"),
  textOutput("export_header_validation"),
  textOutput("e2e_active_columns"),
  textOutput("e2e_persisted_fields"),
  textOutput("e2e_saved_templates"),
  DTOutput("export_header_values_table")
)

server <- function(input, output, session) {
  rv <- reactiveValues(fasta_header_config = pitax_header_default_config())
  export_dialog_state <- reactiveValues(snap = e2e_snap, stage_id = "qc")
  pitax_collect_export_selection_from_input <- function(stage_id) list()
  sys.source("R/server/stages/85_fasta_header_builder.R", envir = environment())

  output$e2e_active_columns <- renderText({
    paste(names(pitax_header_active_edit_table(e2e_snap, export_header_state$config)),
          collapse = ",")
  })
  output$e2e_persisted_fields <- renderText({
    paste(rv$fasta_header_config$fields, collapse = ",")
  })
  output$e2e_saved_templates <- renderText({
    paste(names(rv$fasta_header_config$templates), collapse = ",")
  })
}

shiny::runApp(shinyApp(ui, server), host = "127.0.0.1", port = 47778)
