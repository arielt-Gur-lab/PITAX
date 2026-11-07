# Real Shiny app hosting the structured Header Builder stage + assets.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
source("R/export/availability.R", local = FALSE)
source("R/export/header_templates.R", local = FALSE)

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

# Serve the same www/ assets as production so Playwright can verify HTTP 200.
shiny::addResourcePath("pitax-www", normalizePath("www", winslash = "/", mustWork = TRUE))

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "pitax-www/pitax-header-builder.css"),
    tags$script(src = "pitax-www/pitax-header-builder.js")
  ),
  titlePanel("PITAX header builder integration"),
  selectInput("export_header_saved_template", "Saved profiles",
              choices = c("Sequence Library default" = "prof_default")),
  textInput("export_header_template_name", "Profile name", value = "Sequence Library default"),
  selectInput("export_header_tag_style", "Tag style",
              choices = c("bare" = "bare", "bracket" = "bracket"), selected = "bare"),
  actionButton("export_header_save_template", "Save profile"),
  actionButton("export_header_load_template", "Load profile"),
  actionButton("export_header_save_browser", "Save to this browser"),
  actionButton("export_header_duplicate_profile", "Duplicate"),
  actionButton("export_header_delete_profile", "Delete profile"),
  textInput("export_header_new_field", "New custom field", value = ""),
  actionButton("export_header_add_field", "Create custom field"),
  selectInput("export_header_existing_field", "Custom fields", choices = character(0)),
  actionButton("export_header_delete_field", "Delete custom field"),
  selectInput("export_header_bulk_field", "Field to set", choices = character(0)),
  textInput("export_header_bulk_value", "Value", value = ""),
  selectInput("export_header_bulk_scope", "Apply to",
              choices = c("All records" = "all", "Selected" = "selected"), selected = "all"),
  selectizeInput("export_header_bulk_isolate", "Isolate filter", choices = character(0)),
  selectizeInput("export_header_bulk_locus", "Locus filter", choices = character(0)),
  actionButton("export_header_apply_bulk", "Apply metadata value"),
  div(id = "fhb_persistence_status", class = "fhb-status"),
  div(id = "fhb_library"),
  div(id = "fhb_construction", class = "fhb-construction"),
  div(id = "fhb_preview_template", class = "fhb-preview-box"),
  div(id = "fhb_structure_validation", class = "fhb-validation"),
  textOutput("export_header_registered_fields"),
  textOutput("export_header_active_template"),
  textOutput("export_header_validation"),
  textOutput("e2e_active_columns"),
  textOutput("e2e_persisted_fields"),
  textOutput("e2e_saved_profiles"),
  textOutput("e2e_browser_status"),
  textAreaInput("export_header_template", NULL, value = "", rows = 1),
  DTOutput("export_header_values_table"),
  DTOutput("export_header_preview_table")
)

server <- function(input, output, session) {
  rv <- reactiveValues(fasta_header_config = pitax_header_default_config())
  export_dialog_state <- reactiveValues(snap = e2e_snap, stage_id = "qc")
  pitax_collect_export_selection_from_input <- function(stage_id) list()
  sys.source("R/server/stages/85_fasta_header_builder.R", envir = environment())

  output$e2e_active_columns <- renderText({
    paste(pitax_header_active_metadata_fields(export_header_state$config), collapse = ",")
  })
  output$e2e_persisted_fields <- renderText({
    paste(export_header_state$config$draft$ordered_fields, collapse = ",")
  })
  output$e2e_saved_profiles <- renderText({
    paste(vapply(export_header_state$config$profiles, function(p) p$name, character(1)),
          collapse = ",")
  })
  output$e2e_browser_status <- renderText({
    p <- export_header_state$config$persistence
    paste0("browser=", isTRUE(p$browser_saved), "; status=", p$last_status %||% "")
  })
}

shiny::runApp(shinyApp(ui, server), host = "127.0.0.1", port = 47778)
