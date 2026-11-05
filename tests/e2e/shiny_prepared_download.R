# Minimal running Shiny harness using PITAX production ZIP-button and file-transfer
# helpers. This tests actual HTTP browser downloading, not just generation.
suppressPackageStartupMessages(library(shiny))
source("R/export/dialog_ui.R")
source("R/export/execute.R")

ui <- fluidPage(
  titlePanel("PITAX ZIP download integration"),
  actionButton("prepare", "1. Prepare ZIP"),
  uiOutput("export_prepared_download_ui", inline = TRUE),
  textOutput("status")
)

server <- function(input, output, session) {
  prepared <- reactiveVal(NULL)
  output$status <- renderText({
    if (is.null(prepared())) "Not ready" else "ZIP ready"
  })
  output$export_prepared_download_ui <- renderUI({
    pitax_export_prepared_download_button(prepared())
  })
  observeEvent(input$prepare, {
    folder <- tempfile("test_zip_contents_")
    dir.create(folder)
    writeLines(c(">ISO1_TEF_R organism=En", "ACGTACGT"), file.path(folder, "sample.fasta"))
    archive <- tempfile("test_prepared_", fileext = ".zip")
    zip::zipr(archive, "sample.fasta", root = folder, mode = "mirror")
    stopifnot(file.exists(archive), file.info(archive)$size > 0)
    prepared(archive)
  })
  output$download_export_package <- downloadHandler(
    filename = function() "PITAX_browser_test.zip",
    contentType = "application/zip",
    content = function(file) pitax_export_transfer_prepared_zip(prepared(), file)
  )
  outputOptions(output, "download_export_package", suspendWhenHidden = FALSE)
  session$onSessionEnded(function() {
    p <- isolate(prepared())
    if (!is.null(p) && file.exists(p)) unlink(p)
  })
}

shiny::runApp(shinyApp(ui, server), host = "127.0.0.1", port = 47777)
