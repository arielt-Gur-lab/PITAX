# ============================================================
# PITAX v3.0.0-alpha.9.2
# Application entry point
# ============================================================

source(file.path("R", "app", "bootstrap.R"), local = TRUE)
source(file.path("R", "ui", "components.R"), local = TRUE)
source(file.path("R", "ui", "app_ui.R"), local = TRUE)
source(file.path("R", "server", "app_server.R"), local = TRUE)

shinyApp(ui = ui, server = server)
