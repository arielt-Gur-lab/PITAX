# ============================================================
# PITAX v3.0.0-alpha.10.1
# Application entry point
# ============================================================

source(file.path("R", "app", "bootstrap.R"), local = TRUE, encoding = "UTF-8")
source(file.path("R", "ui", "components.R"), local = TRUE, encoding = "UTF-8")
source(file.path("R", "ui", "app_ui.R"), local = TRUE, encoding = "UTF-8")
source(file.path("R", "server", "app_server.R"), local = TRUE, encoding = "UTF-8")

shinyApp(ui = ui, server = server)
