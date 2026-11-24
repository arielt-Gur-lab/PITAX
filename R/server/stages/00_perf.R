  # Early server stage: opt-in development diagnostics. Never persisted in
  # saved projects and never writes a file unless the user explicitly exports.
  perf_store <- pitax_perf_new(max_events = 600L)
  perf_revision <- reactiveVal(0L)
  perf_publish <- function() {
    perf_revision(isolate(perf_revision()) + 1L)
    invisible(NULL)
  }
  perf_measure <- function(operation, expr, details = "") {
    expr <- substitute(expr)
    frame <- parent.frame()
    if (!isTRUE(perf_store$enabled)) return(eval(expr, frame))
    start <- unname(proc.time()[["elapsed"]])
    on.exit({
      recorded <- pitax_perf_append(
        perf_store, "Server", operation,
        (unname(proc.time()[["elapsed"]]) - start) * 1000, details
      )
      if (isTRUE(recorded)) perf_publish()
    }, add = TRUE)
    eval(expr, frame)
  }
  observeEvent(input$pitax_perf_enabled, {
    perf_store$enabled <- isTRUE(input$pitax_perf_enabled)
    if (perf_store$enabled) perf_publish()
  }, ignoreInit = FALSE)

  observeEvent(input$pitax_perf_browser_event, {
    if (!isTRUE(perf_store$enabled)) return()
    event <- input$pitax_perf_browser_event
    action <- as.character(event$action)[1]
    duration <- suppressWarnings(as.numeric(event$duration_ms)[1])
    if (length(action) != 1L || is.na(action) ||
        !action %in% c("guided.row", "guided.start",
                       "guided.continue", "guided.flag",
                       "guided.confirm", "guided.cancel")) return()
    if (isTRUE(pitax_perf_append(
      perf_store, "Browser", paste0(action, ".to_shiny_idle"), duration,
      "Browser click to Shiny idle; includes round trip, not final paint"
    ))) perf_publish()
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_perf_clear, {
    perf_store$events <- pitax_perf_empty()
    perf_publish()
  }, ignoreInit = TRUE)

  output$pitax_perf_summary <- renderTable({
    req(isTRUE(input$pitax_perf_enabled))
    perf_revision()
    pitax_perf_summary(perf_store$events)
  }, striped = TRUE, spacing = "xs", digits = 2)

  output$pitax_perf_recent <- renderTable({
    req(isTRUE(input$pitax_perf_enabled))
    perf_revision()
    rows <- tail(perf_store$events, 25L)
    if (!nrow(rows)) return(rows)
    rows[nrow(rows):1L, , drop = FALSE]
  }, striped = TRUE, spacing = "xs", digits = 2)

  output$pitax_perf_download <- downloadHandler(
    filename = function() paste0("PITAX_performance_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      utils::write.csv(perf_store$events, file, row.names = FALSE, fileEncoding = "UTF-8")
    },
    contentType = "text/csv"
  )
