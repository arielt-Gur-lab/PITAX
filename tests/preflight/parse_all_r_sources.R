# PITAX preflight: parse every R source before running the numbered test suite.
# This catches syntax errors such as invalid string escapes before any smoke test executes.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "VERSION.txt"))) {
  stop("Preflight must be run from the PITAX Application root.", call. = FALSE)
}

r_files <- unique(c(
  file.path(project_root, "app.R"),
  list.files(file.path(project_root, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(project_root, "tests"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(project_root, "scripts"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
))
r_files <- r_files[file.exists(r_files)]
if (!length(r_files)) stop("No R sources were found for preflight parsing.", call. = FALSE)

failures <- character()
for (path in r_files) {
  result <- tryCatch(
    {
      parse(file = path, keep.source = TRUE)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(result)) {
    rel <- substring(normalizePath(path, winslash = "/", mustWork = TRUE), nchar(project_root) + 2L)
    failures <- c(failures, paste0(rel, ": ", result))
  }
}

if (length(failures)) {
  stop(
    "R-source preflight failed:\n - ",
    paste(failures, collapse = "\n - "),
    call. = FALSE
  )
}

cat("Preflight parse passed for ", length(r_files), " R source files.\n", sep = "")
