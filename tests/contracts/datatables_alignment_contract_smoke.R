# PITAX v3.0.0-alpha.9.2 - global DataTables alignment regression contract.

source(file.path("tests", "helpers", "test_paths.R"))
app_text <- pitax_read_app_contract()
app_lines <- strsplit(app_text, "\n", fixed = TRUE)[[1]]

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

scroll_lines <- app_lines[grepl("scrollX[[:space:]]*=[[:space:]]*TRUE", app_lines)]
assert_true(length(scroll_lines) > 0L, "No horizontally scrolling DataTables were found.")
assert_true(all(grepl("autoWidth[[:space:]]*=[[:space:]]*TRUE", scroll_lines)),
            "Every scrollX DataTable must enable autoWidth so its cloned header and body share computed widths.")

cat("v3.0.0-alpha.9.2 DataTables alignment contracts passed.\n")
