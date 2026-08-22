required_cran <- c("shiny", "DT", "zip", "readxl", "httr2", "plotly", "xml2", "rentrez", "taxize", "openxlsx")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("sangerseqR", quietly = TRUE)) BiocManager::install("sangerseqR")

library(shiny)
library(DT)
library(sangerseqR)
options(shiny.maxRequestSize = 500 * 1024^2)

# Static browser assets are served automatically by Shiny from www/.
PITAX_LOGO_AVAILABLE <- file.exists(file.path("www", "logo.png"))

source(file.path("R", "domain", "sanger", "ab1_evidence.R"), local = TRUE)
source(file.path("R", "domain", "assignment", "stage2_architecture.R"), local = TRUE)
source(file.path("R", "domain", "consensus", "stage3_consensus.R"), local = TRUE)
source(file.path("R", "domain", "multilocus", "stage4_multilocus.R"), local = TRUE)
source(file.path("R", "domain", "sanger", "core_sanger.R"), local = TRUE)
source(file.path("R", "domain", "sanger", "sequence_tools.R"), local = TRUE)
source(file.path("R", "export", "export_tools.R"), local = TRUE)
source(file.path("R", "services", "taxonomy_tools.R"), local = TRUE)

APP_VERSION <- tryCatch(trimws(readLines("VERSION.txt", warn = FALSE)[1]), error = function(e) "3.0.0-alpha.9.2")
APP_VERSION <- ifelse(is.na(APP_VERSION) || !nzchar(APP_VERSION), "3.0.0-alpha.9.2", APP_VERSION)
PROJECT_SCHEMA_VERSION <- 5L
