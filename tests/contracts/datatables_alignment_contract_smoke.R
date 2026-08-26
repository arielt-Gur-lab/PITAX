# PITAX DataTables alignment regression contract.

source(file.path("tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
app_lines <- strsplit(app_text, "\n", fixed = TRUE)[[1]]
css_text <- pitax_read_text("www", "pitax.css")
blast_text <- pitax_read_text("R", "server", "stages", "90_blast.R")

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
must_contain <- function(marker) {
  assert_true(grepl(marker, app_text, fixed = TRUE), paste("Missing DataTables alignment marker:", marker))
}

assert_true(!grepl("table.dataTable { width:100% !important; }", app_text, fixed = TRUE),
            "A global !important width overrides DataTables' computed header/body widths.")
must_contain("$.fn.dataTable.tables({ visible: true, api: true })")
must_contain("api.columns.adjust()")
must_contain("shown.bs.tab")
must_contain("init.dt draw.dt")
must_contain("scroll.pitaxAlignment")
must_contain(".dataTables_scrollBody table.dataTable { box-sizing:border-box; }")

# BLAST jobs table: autoWidth=FALSE is intentional; CSS keeps scrollX header/body aligned.
assert_true(grepl("#blast_jobs_table .dataTables_scrollHeadInner", css_text, fixed = TRUE),
            "BLAST jobs table needs dedicated scrollHeadInner CSS sync.")
assert_true(grepl("autoWidth = FALSE", blast_text, fixed = TRUE),
            "BLAST jobs table must keep autoWidth=FALSE with CSS sync.")

scroll_lines <- app_lines[grepl("scrollX[[:space:]]*=[[:space:]]*TRUE", app_lines)]
assert_true(length(scroll_lines) > 0L, "No horizontally scrolling DataTables were found.")

# Same-line scrollX+autoWidth=TRUE remains the default. Bare "scrollX = TRUE," lines
# are allowed only for the BLAST multiline options block (autoWidth=FALSE nearby).
for (line in scroll_lines) {
  if (grepl("autoWidth[[:space:]]*=[[:space:]]*TRUE", line)) next
  if (grepl("autoWidth[[:space:]]*=[[:space:]]*FALSE", line)) next
  if (grepl("^\\s*scrollX[[:space:]]*=[[:space:]]*TRUE\\s*,?\\s*$", line)) next
  stop("Unexpected scrollX line without autoWidth pairing: ", line, call. = FALSE)
}

cat("DataTables alignment contracts passed.\n")
