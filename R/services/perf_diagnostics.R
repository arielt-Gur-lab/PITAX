# Session-only, opt-in developer performance probe.
# Do not attach this to the saved .sangerproject schema. The log contains
# timings and operation names only: no read names, DNA, paths or BLAST data.
pitax_perf_new <- function(max_events = 600L) {
  e <- new.env(parent = emptyenv())
  e$enabled <- FALSE
  e$max_events <- max(10L, as.integer(max_events)[1])
  e$events <- pitax_perf_empty()
  e
}

pitax_perf_empty <- function() {
  data.frame(
    Time = character(), Source = character(), Operation = character(),
    Duration_ms = numeric(), Details = character(), stringsAsFactors = FALSE
  )
}

pitax_perf_append <- function(store, source, operation, milliseconds, details = "") {
  if (!is.environment(store) || !isTRUE(store$enabled)) return(invisible(FALSE))
  duration <- suppressWarnings(as.numeric(milliseconds)[1])
  if (!length(duration) || !is.finite(duration) || duration < 0) return(invisible(FALSE))
  source <- as.character(source)[1]
  operation <- as.character(operation)[1]
  if (!source %in% c("Server", "Browser") || is.na(operation) ||
      !nzchar(operation)) return(invisible(FALSE))
  # Bound every field; this must never become a sequence/metadata logger.
  operation <- substr(operation, 1L, 72L)
  details <- substr(as.character(details)[1], 1L, 100L)
  if (is.na(details)) details <- ""
  row <- data.frame(
    Time = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3 %z"),
    Source = source, Operation = operation,
    Duration_ms = round(duration, 2), Details = details,
    stringsAsFactors = FALSE
  )
  store$events <- rbind(store$events, row)
  over <- nrow(store$events) - store$max_events
  if (over > 0L) {
    store$events <- store$events[(over + 1L):nrow(store$events), , drop = FALSE]
    rownames(store$events) <- NULL
  }
  invisible(TRUE)
}

pitax_perf_measure <- function(store, operation, expr, on_record = NULL, details = "") {
  expr <- substitute(expr)
  frame <- parent.frame()
  if (!is.environment(store) || !isTRUE(store$enabled)) {
    return(eval(expr, frame))
  }
  start <- unname(proc.time()[["elapsed"]])
  on.exit({
    elapsed <- (unname(proc.time()[["elapsed"]]) - start) * 1000
    recorded <- pitax_perf_append(store, "Server", operation, elapsed, details)
    if (isTRUE(recorded) && is.function(on_record)) on_record()
  }, add = TRUE)
  eval(expr, frame)
}

pitax_perf_summary <- function(events) {
  if (!is.data.frame(events) || !nrow(events)) {
    return(data.frame(Operation = character(), N = integer(),
                      Median_ms = numeric(), P95_ms = numeric(),
                      Max_ms = numeric(), stringsAsFactors = FALSE))
  }
  groups <- split(events$Duration_ms, paste(events$Source, events$Operation, sep = " / "))
  out <- do.call(rbind, lapply(names(groups), function(label) {
    t <- groups[[label]]
    data.frame(
      Operation = label, N = length(t),
      Median_ms = round(stats::median(t), 1),
      P95_ms = round(as.numeric(stats::quantile(t, 0.95, names = FALSE)), 1),
      Max_ms = round(max(t), 1), stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$P95_ms, decreasing = TRUE), , drop = FALSE]
}
