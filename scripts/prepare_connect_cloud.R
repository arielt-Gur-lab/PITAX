# PITAX - prepare a Posit Connect Cloud deployment manifest.
# Run with source("scripts/prepare_connect_cloud.R") from the project root
# after the app runs successfully locally.

if (!file.exists("app.R")) {
  stop("Run this script from the PITAX project root (the folder containing app.R).")
}

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

rsconnect::writeManifest(appDir = ".")
cat("Created manifest.json for Posit Connect Cloud.\n")
cat("Commit manifest.json to GitHub, then publish the repository as Shiny for R.\n")
