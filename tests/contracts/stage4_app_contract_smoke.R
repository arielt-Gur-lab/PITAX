# PITAX v3.0.0-alpha.9.2 - Stage 4 Shiny integration contract.

source(file.path("tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
stage4_text <- pitax_read_text("R", "domain", "multilocus", "stage4_multilocus.R")

must_contain <- function(text, marker) {
  if (!grepl(marker, text, fixed = TRUE)) stop(paste("Missing Stage 4 integration marker:", marker), call. = FALSE)
}

must_contain(app_text, 'pitax_source(file.path("R", "domain", "multilocus", "stage4_multilocus.R"), local = TRUE)')
must_contain(app_text, "PROJECT_SCHEMA_VERSION <- 6L")
must_contain(app_text, 'tabPanel("Multi-locus", value = "multilocus"')
must_contain(app_text, 'actionButton("build_multilocus_profile", "Build / rebuild profile"')
must_contain(app_text, 'rv$multilocus_imports')
must_contain(app_text, 'output$multilocus_import_queue <- renderUI({')
must_contain(app_text, 'for (i in seq_along(imports)) {')
must_contain(app_text, 'stage4_build_profile(projects, source_names = source_names, source_md5 = source_md5)')
# Build must use the accumulated queue, not re-read only the latest fileInput value.
build_pos <- regexpr('observeEvent(input$build_multilocus_profile, {', app_text, fixed = TRUE)[1]
if (build_pos < 1L) stop("Missing build_multilocus_profile observer.", call. = FALSE)
build_tail <- substr(app_text, build_pos, nchar(app_text))
next_observe <- regexpr('\n  observe({', build_tail, fixed = TRUE)[1]
if (next_observe < 1L) stop("Could not locate the isolate-selector observe after Build.", call. = FALSE)
build_block <- substr(build_tail, 1, next_observe)
if (grepl('uploads <- input$multilocus_projects', build_block, fixed = TRUE)) {
  stop("Build must not re-read input$multilocus_projects; use rv$multilocus_imports.", call. = FALSE)
}
must_contain(build_block, 'imports <- rv$multilocus_imports')
must_contain(app_text, 'multilocus_profile = rv$multilocus_profile')
must_contain(app_text, 'if (source_schema == 5L) {')
must_contain(app_text, 'assay_migrate_schema5_state(st)')
must_contain(app_text, 'schema5_migration_error')
must_contain(app_text, 'rv$multilocus_profile <- stage4_ensure_profile')
must_contain(app_text, 'stage4_current_project_matches(profile, make_project_bundle(), "Current session")')
must_contain(app_text, 'write_stage4_checkpoint_zip(file, rv$multilocus_profile)')
must_contain(app_text, 'selectInput("multilocus_isolate", "Isolate"')
must_contain(app_text, 'output$multilocus_overview_cards <- renderUI({')
must_contain(app_text, 'output$multilocus_profile_hero <- renderUI({')
must_contain(app_text, 'output$multilocus_locus_cards <- renderUI({')
must_contain(app_text, 'stage4_display_call(row$Taxonomy_Status, row$Recommended_Identification)')
must_contain(app_text, 'text = ~Locus, hovertext = ~Hover')
must_contain(app_text, 'stage4_isolate_evidence(rv$multilocus_profile, input$multilocus_isolate)')
must_contain(stage4_text, 'algorithm = "pitax-evidence-preserving-profile-v1"')
must_contain(stage4_text, 'stage4_profile_overview <- function(profile)')
must_contain(stage4_text, 'stage4_isolate_evidence <- function(profile, isolate)')
must_contain(stage4_text, 'Duplicate Isolate/Locus evidence is not allowed')
must_contain(stage4_text, '"GENUS_CONFLICT"')
must_contain(stage4_text, '"CONCORDANT_SPECIES"')
must_contain(stage4_text, 'The conflict is retained rather than voted away')

cat("v3.0.0-alpha.9.2 Stage 4 app integration contract passed.\n")
