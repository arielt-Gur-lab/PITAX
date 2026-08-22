# ============================================================
# PITAX v3.0.0-alpha.10.2
# Application entry point
# ============================================================

pitax_source <- function(path, local) {
  target_environment <- if (isTRUE(local)) parent.frame() else local
  arguments <- list(file = path, local = target_environment)
  if (identical(.Platform$OS.type, "windows")) arguments$encoding <- "UTF-8"
  do.call(base::source, arguments)
}

pitax_source(file.path("R", "app", "bootstrap.R"), local = TRUE)
pitax_source(file.path("R", "ui", "components.R"), local = TRUE)
pitax_source(file.path("R", "ui", "app_ui.R"), local = TRUE)
pitax_source(file.path("R", "server", "app_server.R"), local = TRUE)

shinyApp(ui = ui, server = server)
