  # Early server stage: opt-in development diagnostics. Never persisted in
  # saved projects and never writes a file unless the user explicitly exports.
  perf_store <- pitax_perf_new(max_events = 600L)
  perf_revision <- reactiveVal(0L)
  # Coalesce table refreshes so Continue's reactive flush is not inflated by
  # re-rendering the diagnostics summary/recent tables in the same cycle.
  perf_publish_scheduled <- FALSE
  perf_publish <- function() {
    if (isTRUE(perf_publish_scheduled)) return(invisible(NULL))
    perf_publish_scheduled <<- TRUE
    session$onFlushed(function() {
      perf_publish_scheduled <<- FALSE
      perf_revision(isolate(perf_revision()) + 1L)
    }, once = TRUE)
    invisible(NULL)
  }
  perf_measure <- function(operation, expr, details = "", publish = TRUE) {
    expr <- substitute(expr)
    frame <- parent.frame()
    if (!isTRUE(perf_store$enabled)) return(eval(expr, frame))
    start <- unname(proc.time()[["elapsed"]])
    on.exit({
      recorded <- pitax_perf_append(
        perf_store, "Server", operation,
        (unname(proc.time()[["elapsed"]]) - start) * 1000, details
      )
      # Diagnostic table renders must not re-publish: that would invalidate
      # themselves and loop. Their rows appear on the next natural publish.
      if (isTRUE(recorded) && isTRUE(publish)) perf_publish()
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
    details <- as.character(event$details)[1]
    if (length(details) != 1L || is.na(details)) details <- ""
    if (!nzchar(details)) {
      details <- "Browser click to Shiny idle; includes round trip, not final paint"
    }
    if (isTRUE(pitax_perf_append(
      perf_store, "Browser", paste0(action, ".to_shiny_idle"), duration, details
    ))) perf_publish()

    # Optional phase splits from the browser probe (milliseconds).
    phases <- list(
      to_first_message = suppressWarnings(as.numeric(event$to_first_message_ms)[1]),
      to_last_tracked_value = suppressWarnings(as.numeric(event$to_last_value_ms)[1])
    )
    for (nm in names(phases)) {
      ms <- phases[[nm]]
      if (!length(ms) || !is.finite(ms) || ms < 0) next
      pitax_perf_append(
        perf_store, "Browser", paste0(action, ".", nm), ms,
        "Phase within click-to-idle; not final paint"
      )
    }
    if (any(vapply(phases, function(x) length(x) && is.finite(x), logical(1)))) {
      perf_publish()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$pitax_perf_clear, {
    perf_store$events <- pitax_perf_empty()
    perf_publish()
  }, ignoreInit = TRUE)

  output$pitax_perf_summary <- renderTable({
    req(isTRUE(input$pitax_perf_enabled))
    perf_revision()
    perf_measure("qc.perf_summary_render", {
      pitax_perf_summary(perf_store$events)
    }, publish = FALSE)
  }, striped = TRUE, spacing = "xs", digits = 2)

  output$pitax_perf_recent <- renderTable({
    req(isTRUE(input$pitax_perf_enabled))
    perf_revision()
    perf_measure("qc.perf_recent_render", {
      rows <- tail(perf_store$events, 25L)
      if (!nrow(rows)) return(rows)
      rows[nrow(rows):1L, , drop = FALSE]
    }, publish = FALSE)
  }, striped = TRUE, spacing = "xs", digits = 2)

  output$pitax_perf_download <- downloadHandler(
    filename = function() paste0("PITAX_performance_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      utils::write.csv(perf_store$events, file, row.names = FALSE, fileEncoding = "UTF-8")
    },
    contentType = "text/csv"
  )
