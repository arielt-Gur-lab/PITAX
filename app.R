# ============================================================
# PITAX v3.0.0-alpha.6
# Stage 2: isolate / locus / read architecture
# ============================================================

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

# Serve the PITAX logo from the project root without exposing the rest of the source tree.
PITAX_LOGO_AVAILABLE <- file.exists("logo.png")
if (PITAX_LOGO_AVAILABLE) {
  PITAX_ASSET_DIR <- file.path(tempdir(), "pitax-assets")
  dir.create(PITAX_ASSET_DIR, recursive = TRUE, showWarnings = FALSE)
  file.copy("logo.png", file.path(PITAX_ASSET_DIR, "logo.png"), overwrite = TRUE)
  addResourcePath("pitax-assets", PITAX_ASSET_DIR)
}

source(file.path("R", "ab1_evidence.R"), local = TRUE)
source(file.path("R", "stage2_architecture.R"), local = TRUE)
source(file.path("R", "core_sanger.R"), local = TRUE)
source(file.path("R", "sequence_tools.R"), local = TRUE)
source(file.path("R", "export_tools.R"), local = TRUE)
source(file.path("R", "taxonomy_tools.R"), local = TRUE)

APP_VERSION <- tryCatch(trimws(readLines("VERSION.txt", warn = FALSE)[1]), error = function(e) "3.0.0-alpha.6")
APP_VERSION <- ifelse(is.na(APP_VERSION) || !nzchar(APP_VERSION), "3.0.0-alpha.6", APP_VERSION)
PROJECT_SCHEMA_VERSION <- 3L

# ============================================================
# UI helpers
# ============================================================

panel_box <- function(...) div(class = "panel-box", ...)
section_title <- function(x) div(class = "section-title", x)
info_tip <- function(text) tags$span(class = "info-tip", title = text, "ⓘ")

stage_heading <- function(icon_name, title, subtitle, badge = NULL) {
  div(
    class = "stage-heading",
    div(class = "stage-heading-icon", icon(icon_name)),
    div(
      class = "stage-heading-copy",
      h3(title),
      div(class = "stage-heading-subtitle", subtitle)
    ),
    if (!is.null(badge)) div(class = "stage-heading-badge", badge) else NULL
  )
}

card_title <- function(title, tip = NULL, icon_name = NULL) {
  div(
    class = "card-title-row",
    if (!is.null(icon_name)) div(class = "card-title-icon", icon(icon_name)) else NULL,
    div(class = "card-title-text", title),
    if (!is.null(tip)) info_tip(tip) else NULL
  )
}

rbind_fill <- function(a, b) {
  if (is.null(a) || !is.data.frame(a) || !nrow(a)) return(b)
  if (is.null(b) || !is.data.frame(b) || !nrow(b)) return(a)
  cols <- union(names(a), names(b))
  for (nm in setdiff(cols, names(a))) a[[nm]] <- NA
  for (nm in setdiff(cols, names(b))) b[[nm]] <- NA
  a <- a[, cols, drop = FALSE]
  b <- b[, cols, drop = FALSE]
  rbind(a, b)
}

stage_topbar <- function(...) {
  div(class = "stage-topbar", ...)
}

pipeline_stage_footer <- function(current_step) {
  labels <- c(
    "Upload",
    "Assay settings",
    "Rename & assign",
    "Trim & QC",
    "Export",
    "NCBI BLAST",
    "Taxonomic summary"
  )

  pieces <- list()
  for (i in seq_along(labels)) {
    state <- if (i < current_step) "done" else if (i == current_step) "current" else "future"
    pieces[[length(pieces) + 1]] <- div(
      class = paste("schema-stage", state),
      div(class = "schema-dot"),
      div(class = "schema-label", paste0(i, ". ", labels[[i]]))
    )
    if (i < length(labels)) {
      pieces[[length(pieces) + 1]] <- div(
        class = paste("schema-line", if (i < current_step) "done" else "")
      )
    }
  }

  div(
    class = "pipeline-schema-wrap",
    div(class = "pipeline-schema-title", "Workflow position"),
    div(class = "pipeline-schema", pieces)
  )
}

ui <- fluidPage(
  tags$head(
    tags$title("PITAX — Taxonomic Identification Tool"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('openUrl', function(message) {
        window.open(message.url, '_blank');
      });
      Shiny.addCustomMessageHandler('copyText', function(message) {
        navigator.clipboard.writeText(message.text);
      });

      (function() {
        var loaderShownAt = 0;
        var loaderFallback = null;
        function showStepLoader(text) {
          loaderShownAt = Date.now();
          $('#app_loading_text').text(text || 'Loading workspace…');
          $('#app_loading_overlay').addClass('visible').attr('aria-hidden', 'false');
          clearTimeout(loaderFallback);
          loaderFallback = setTimeout(hideStepLoader, 20000);
        }
        function hideStepLoader() {
          var elapsed = Date.now() - loaderShownAt;
          var wait = Math.max(0, 280 - elapsed);
          setTimeout(function() {
            $('#app_loading_overlay').removeClass('visible').attr('aria-hidden', 'true');
          }, wait);
          clearTimeout(loaderFallback);
        }
        Shiny.addCustomMessageHandler('showLoader', function(message) {
          showStepLoader(message && message.text ? message.text : 'Loading workspace…');
        });
        Shiny.addCustomMessageHandler('hideLoader', function() { hideStepLoader(); });
        $(document).on('shiny:busy', function() { $('#shiny_activity_bar').addClass('visible'); });
        $(document).on('shiny:idle', function() {
          $('#shiny_activity_bar').removeClass('visible');
          if ($('#app_loading_overlay').hasClass('visible')) hideStepLoader();
        });
        $(document).on('click', '#pipeline_step > li > a', function() {
          showStepLoader('Loading step…');
        });
      })();

    ")),
    tags$style(HTML("
      body { background:#f5f7fa; }
      .pipeline-container { max-width:1450px; margin:auto; }
      .app-header { background:white; padding:22px 28px; margin:20px 0 18px; border-radius:12px; border:1px solid #e5e7eb; }
      .app-header h2 { margin:0 0 5px; }
      .app-header p { margin:0; color:#6b7280; }
      .panel-box { background:white; padding:22px; border-radius:12px; border:1px solid #e5e7eb; margin-bottom:18px; }
      .section-title { font-size:20px; font-weight:600; margin-bottom:16px; }
      .summary-card { background:#f8fafc; border:1px solid #e5e7eb; border-radius:10px; padding:15px; margin-bottom:10px; }
      .summary-number { font-size:26px; font-weight:700; }
      .summary-label { color:#6b7280; }
      .settings-note { color:#6b7280; font-size:13px; }
      .status-ok { color:#15803d; font-weight:600; }
      .status-warning { color:#b45309; font-weight:600; }
      .status-error { color:#b91c1c; font-weight:600; }
      .blast-action-row { display:flex; flex-wrap:wrap; gap:8px; margin-bottom:10px; }
      .checkpoint { background:#fffaf0; border:1px solid #f3d7a0; }
      .mono { font-family:Consolas,Monaco,monospace; white-space:pre-wrap; word-break:break-all; background:#f8fafc; padding:12px; border-radius:8px; border:1px solid #e5e7eb; }
      table.dataTable td input { color:#111827 !important; background:#fff !important; }
      table.dataTable td input:focus { color:#111827 !important; background:#fff !important; }
      textarea.form-control { font-family:Consolas,Monaco,monospace; }
      .btn { border-radius:7px; margin-right:6px; margin-bottom:6px; }
      .experimental-box { background:#fffdf5; border:1px solid #f1dfae; border-radius:8px; padding:10px 12px; margin-top:8px; }
      .stage1-evidence-card { border-style:dashed; }
      .stage1-audit-details { margin-top:8px; }
      .stage1-audit-toggle { cursor:pointer; font-weight:600; padding:10px 0; color:#374151; }
      .stage1-audit-details[open] .stage1-audit-toggle { margin-bottom:10px; }
      .stage-topbar { background:#ffffff; border:1px solid #dbe3ee; border-radius:10px; padding:12px 14px; margin:0 0 16px; display:flex; align-items:center; gap:8px; flex-wrap:wrap; }
      .stage-topbar .btn { margin-bottom:0; }
      .pipeline-schema-wrap { margin:30px 4px 12px; padding:8px 4px 4px; }
      .pipeline-schema-title { font-size:11px; color:#64748b; font-weight:600; text-transform:uppercase; letter-spacing:.05em; margin-bottom:14px; }
      .pipeline-schema { display:flex; align-items:flex-start; width:100%; overflow-x:auto; padding:2px 2px 8px; }
      .schema-stage { min-width:92px; max-width:128px; display:flex; flex-direction:column; align-items:center; text-align:center; flex:0 0 auto; }
      .schema-dot { width:13px; height:13px; border-radius:50%; background:#cbd5e1; border:2px solid #f5f7fa; box-shadow:0 0 0 1px #cbd5e1; margin-top:0; }
      .schema-label { margin-top:8px; font-size:11px; line-height:1.25; color:#94a3b8; font-weight:500; }
      .schema-line { flex:1 1 42px; min-width:28px; height:2px; background:#cbd5e1; margin-top:6px; }
      .schema-stage.done .schema-dot { background:#65a30d; box-shadow:0 0 0 1px #65a30d; }
      .schema-stage.done .schema-label { color:#64748b; }
      .schema-line.done { background:#84cc16; }
      .schema-stage.current .schema-dot { width:17px; height:17px; margin-top:-2px; background:#2563eb; box-shadow:0 0 0 3px rgba(37,99,235,.14); border-color:#ffffff; }
      .schema-stage.current .schema-label { color:#1e3a8a; font-weight:700; }
      .taxonomy-hero { background:#f8fafc; border-left:4px solid #2563eb; padding:16px 18px; margin-bottom:16px; }
      .taxonomy-id { font-size:28px; font-weight:700; line-height:1.15; }
      .taxonomy-sub { color:#64748b; margin-top:4px; }
      .confidence-high { color:#15803d; font-weight:700; }
      .confidence-moderate { color:#a16207; font-weight:700; }
      .confidence-low { color:#b91c1c; font-weight:700; }
      .tax-note { background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:12px 14px; color:#475569; }
      table.dataTable td.dt-taxon { min-width:150px; max-width:230px; white-space:normal; line-height:1.25; vertical-align:top; font-style:italic; }
      .raw-response-details { margin-top:14px; }
      .raw-response-details summary { display:inline-block; cursor:pointer; user-select:none; padding:8px 12px; background:#f8fafc; border:1px solid #cbd5e1; border-radius:7px; color:#334155; font-weight:600; }
      .raw-response-details summary:hover { background:#eef2f7; }
      .raw-response-details[open] summary { margin-bottom:10px; background:#eaf2f8; border-color:#9dbbd3; }
      table.dataTable td.dt-organism { min-width:170px; max-width:260px; white-space:nowrap; line-height:1.25; vertical-align:top; font-style:italic; }
      table.dataTable td.dt-hit-title { min-width:300px; max-width:520px; white-space:normal; line-height:1.25; vertical-align:top; }
      table.dataTable td.dt-nowrap { white-space:nowrap; vertical-align:top; }
      table.dataTable td.dt-wrap { white-space:normal !important; line-height:1.3; vertical-align:top; }
      table.dataTable td.dt-wrap-compact { white-space:normal !important; line-height:1.18; vertical-align:top; overflow-wrap:anywhere; }
      table.dataTable td.dt-taxon-compact { white-space:normal !important; line-height:1.16; vertical-align:top; font-style:italic; overflow-wrap:anywhere; }
      table.dataTable td.dt-organism-compact { white-space:normal !important; line-height:1.16; vertical-align:top; font-style:italic; overflow-wrap:anywhere; }
      table.dataTable td.dt-compact-number { white-space:nowrap; text-align:center; vertical-align:middle; }
      table.dataTable td.dt-comment { min-width:340px; max-width:480px; white-space:normal !important; line-height:1.32; vertical-align:top; }
      .taxonomy-hits-compact.dataTable { table-layout:fixed !important; width:100% !important; font-size:11px; }
      .taxonomy-hits-compact.dataTable thead th { white-space:normal !important; line-height:1.15; padding-left:5px !important; padding-right:5px !important; }
      .taxonomy-hits-compact.dataTable tbody td { padding-left:5px !important; padding-right:5px !important; }
      .species-profile-table.dataTable td { white-space:normal; }
      .help-card { background:#ffffff; border:1px solid #e2e8f0; border-radius:10px; padding:20px 22px; margin-bottom:18px; line-height:1.55; }
      .help-card h3 { margin-top:0; margin-bottom:14px; font-size:19px; }
      .help-card p { margin-bottom:10px; }
      .help-flow { font-family:Consolas,Monaco,monospace; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:16px; white-space:pre-wrap; line-height:1.7; }
      #about_tabs { margin-bottom:20px; border-bottom:1px solid #cbd5e1; }
      #about_tabs > li > a { padding:10px 16px; font-weight:600; color:#475569; }
      #about_tabs > li.active > a { color:#0f172a; background:#f8fafc; border-bottom-color:#f8fafc; }
      .about-section { max-width:1080px; margin:0 auto; padding:4px 2px 12px 2px; }
      .about-lead { font-size:15px; color:#475569; line-height:1.65; margin-bottom:16px; }
      .method-step { display:grid; grid-template-columns:46px 1fr; gap:14px; padding:16px 0; border-bottom:1px solid #e2e8f0; }
      .method-step:last-child { border-bottom:none; }
      .method-num { width:36px; height:36px; border-radius:50%; display:flex; align-items:center; justify-content:center; background:#e2e8f0; color:#0f172a; font-weight:700; }
      .method-step h4 { margin:0 0 6px 0; font-size:16px; }
      .method-step p { margin:0; color:#475569; line-height:1.55; }
      .evidence-badge { display:inline-block; padding:3px 8px; border-radius:999px; font-size:11px; font-weight:700; margin-left:6px; vertical-align:1px; }
      .evidence-published { background:#dcfce7; color:#166534; }
      .evidence-heuristic { background:#fef3c7; color:#92400e; }
      .reference-card { background:#ffffff; border:1px solid #e2e8f0; border-left:4px solid #94a3b8; border-radius:8px; padding:15px 17px; margin-bottom:12px; }
      .landscape-explainer { background:#f8fafc; border-left:4px solid #64748b; border-radius:6px; padding:10px 14px; margin:4px 0 12px 0; color:#475569; line-height:1.5; }
      .landscape-explainer p { margin:0; }
      .reference-card h4 { margin:0 0 7px 0; font-size:15px; }
      .reference-card p { margin:0 0 8px 0; color:#475569; line-height:1.5; }
      .reference-card a { font-weight:600; }
      .about-callout { background:#f8fafc; border:1px solid #cbd5e1; border-radius:9px; padding:14px 16px; margin:14px 0; line-height:1.55; }
      .team-summary-note { background:#eef7ff; border:1px solid #bfdbfe; border-radius:8px; padding:12px 14px; color:#334155; }
      .peak-flag-summary { display:flex; gap:10px; flex-wrap:wrap; margin:8px 0 12px 0; }
      .peak-flag-pill { border:1px solid #dbe3ec; border-radius:999px; padding:5px 10px; background:#f8fafc; font-size:12px; color:#334155; }
      .peak-flag-pill strong { color:#111827; }
      .peak-flag-help { background:#fffaf0; border:1px solid #fed7aa; border-radius:8px; padding:10px 12px; color:#7c2d12; font-size:12px; }
      .info-tip { display:inline-flex; align-items:center; justify-content:center; width:18px; height:18px; margin-left:6px; border:1px solid #94a3b8; border-radius:50%; color:#64748b; font-size:11px; font-weight:700; cursor:help; vertical-align:1px; }
      .curation-toolbar { display:flex; gap:10px; align-items:flex-end; flex-wrap:wrap; margin:2px 0 8px; }
      .curation-toolbar .form-group { margin-bottom:0; }
      .curation-toolbar .shiny-input-radiogroup { min-width:250px; }
      .curation-actions { display:flex; gap:7px; flex-wrap:wrap; align-items:center; }
      .curation-history-note { color:#64748b; font-size:12px; margin-top:6px; }
      .manual-edit-badge { display:inline-block; padding:3px 8px; border-radius:999px; background:#eef2ff; color:#4338ca; font-size:11px; font-weight:700; }
      .status-note { background:#f8fafc; border-left:3px solid #94a3b8; padding:8px 10px; margin:8px 0; color:#475569; font-size:12px; }
      .project-bar { background:#ffffff; border:1px solid #dbe3ee; border-radius:10px; padding:12px 14px; margin:-2px 0 18px; display:flex; align-items:center; gap:12px; flex-wrap:wrap; }
      .project-bar-title { font-weight:700; color:#334155; margin-right:4px; }
      .project-bar .form-group { margin-bottom:0; }
      .project-bar .shiny-input-container { width:300px; max-width:100%; }
      .project-status { color:#64748b; font-size:12px; flex:1 1 260px; }
      #shiny_activity_bar { position:fixed; top:0; left:0; right:0; height:3px; z-index:10020; opacity:0; pointer-events:none; overflow:hidden; background:transparent; transition:opacity .12s ease; }
      #shiny_activity_bar.visible { opacity:1; }
      #shiny_activity_bar::after { content:''; display:block; width:34%; height:100%; background:#2563eb; animation:activitySlide 1.05s ease-in-out infinite; }
      @keyframes activitySlide { 0% { transform:translateX(-110%); } 100% { transform:translateX(320%); } }
      .app-loading-overlay { position:fixed; inset:0; z-index:10010; display:flex; align-items:center; justify-content:center; background:rgba(244,247,251,.82); backdrop-filter:blur(2px); opacity:0; visibility:hidden; pointer-events:none; transition:opacity .14s ease, visibility .14s ease; }
      .app-loading-overlay.visible { opacity:1; visibility:visible; pointer-events:auto; }
      .app-loading-card { display:flex; align-items:center; gap:14px; min-width:250px; padding:16px 20px; border-radius:14px; background:#fff; border:1px solid #dfe6ef; box-shadow:0 12px 34px rgba(15,23,42,.12); color:#334155; font-weight:650; }
      .app-loading-spinner { width:25px; height:25px; border-radius:50%; border:3px solid #dbeafe; border-top-color:#2563eb; animation:spinLoader .72s linear infinite; flex:0 0 auto; }
      @keyframes spinLoader { to { transform:rotate(360deg); } }

      /* v2.12 visual system ------------------------------------------------ */
      :root {
        --ui-bg:#f4f7fb;
        --ui-card:#ffffff;
        --ui-border:#dfe6ef;
        --ui-text:#132238;
        --ui-muted:#64748b;
        --ui-blue:#2563eb;
        --ui-blue-soft:#eff6ff;
        --ui-green:#159447;
        --ui-green-soft:#ecfdf3;
        --ui-amber:#d97706;
        --ui-amber-soft:#fff7e6;
        --ui-red:#c2413b;
        --ui-purple:#6d4bd2;
        --ui-cyan:#159bb7;
        --ui-shadow:0 5px 18px rgba(15,23,42,.055);
      }
      body { background:var(--ui-bg); color:var(--ui-text); font-family:Inter,-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Arial,sans-serif; }
      #changelog_text { white-space:pre-wrap !important; overflow-wrap:anywhere; word-break:normal; overflow-x:hidden; max-width:100%; font-family:'Aptos','Segoe UI',Arial,sans-serif; line-height:1.55; }
      #changelog_text pre { white-space:pre-wrap !important; overflow-wrap:anywhere; word-break:normal; overflow-x:hidden; max-width:100%; }
      .pipeline-container { max-width:1540px; padding:0 14px 28px; }
      .app-header { border:1px solid var(--ui-border); border-radius:16px; padding:24px 30px; box-shadow:var(--ui-shadow); background:linear-gradient(135deg,#ffffff 0%,#fbfdff 100%); }
      .app-header h2 { font-size:28px; font-weight:750; letter-spacing:-.02em; color:#0f1f35; }
      .panel-box { border:1px solid var(--ui-border); border-radius:14px; box-shadow:var(--ui-shadow); padding:22px 24px; }
      .section-title { color:#17263a; font-size:18px; font-weight:750; letter-spacing:-.01em; display:flex; align-items:center; gap:7px; }
      .stage-topbar { border-color:var(--ui-border); border-radius:12px; box-shadow:0 2px 8px rgba(15,23,42,.035); }
      .project-bar { border-color:var(--ui-border); border-radius:12px; box-shadow:0 2px 8px rgba(15,23,42,.035); }
      .form-control, .selectize-input { border-radius:8px !important; border-color:#cfd9e6 !important; box-shadow:none !important; }
      .form-control:focus, .selectize-input.focus { border-color:#7aa7f7 !important; box-shadow:0 0 0 3px rgba(37,99,235,.09) !important; }
      .btn { border-radius:8px; font-weight:600; box-shadow:none; }
      .btn-primary { background:#2563eb; border-color:#2563eb; }
      .btn-success { background:#168746; border-color:#168746; }
      #pipeline_step > .nav-tabs { border-bottom:1px solid #d8e1ed; margin-bottom:18px; display:flex; flex-wrap:wrap; gap:3px; }
      #pipeline_step > .nav-tabs > li > a { border:none; color:#64748b; font-weight:650; padding:10px 13px; border-radius:8px 8px 0 0; }
      #pipeline_step > .nav-tabs > li.active > a,
      #pipeline_step > .nav-tabs > li.active > a:hover { color:#1d4ed8; background:#ffffff; border:none; box-shadow:inset 0 -3px 0 #2563eb; }
      table.dataTable { border-collapse:separate !important; border-spacing:0 !important; }
      table.dataTable thead th { color:#334155; font-size:12px; font-weight:700; border-bottom:1px solid #dbe4ef !important; background:#fbfcfe; }
      table.dataTable tbody td { border-top:none !important; border-bottom:1px solid #e8edf4; vertical-align:middle; }
      table.dataTable tbody tr:hover { background:#f8fbff !important; }
      .dataTables_wrapper .dataTables_paginate .paginate_button.current { border:1px solid #2563eb !important; background:#fff !important; color:#2563eb !important; border-radius:6px; }


      /* v2.13 cross-stage visual system ----------------------------------- */
      body {
        background:var(--ui-bg);
        color:var(--ui-text);
        font-family:'Aptos','Segoe UI',Arial,sans-serif;
        font-size:14px;
      }
      button, input, select, .selectize-input, .btn, .nav, table, .dataTables_wrapper {
        font-family:'Aptos','Segoe UI',Arial,sans-serif !important;
      }
      textarea.form-control, .mono, pre, code { font-family:Consolas,'Courier New',monospace !important; }
      .pipeline-container { max-width:1580px; padding:0 18px 34px; }
      .app-shell-header { display:flex; align-items:center; gap:18px; min-height:104px; padding:16px 24px; }
      .app-brand-mark { width:48px; height:48px; border-radius:14px; background:#eff6ff; color:#2563eb; display:flex; align-items:center; justify-content:center; font-size:21px; flex:0 0 auto; }
      .app-brand-logo { display:block; width:min(420px,48vw); max-height:82px; object-fit:contain; object-position:left center; }
      .app-brand-copy { flex:1 1 auto; min-width:0; }
      .app-brand-copy h2 { margin:0 0 2px; font-size:25px; }
      .app-brand-copy p { margin:2px 0 0; font-size:12.5px; color:#6b7c92; letter-spacing:.01em; }
      .app-version-badge { padding:6px 10px; border-radius:999px; background:#f1f5f9; color:#52647a; font-size:12px; font-weight:750; }
      .project-session-card { padding:10px 13px; gap:9px; }
      .project-bar-title { display:flex; align-items:center; gap:7px; min-width:max-content; }
      .project-session-card .shiny-input-container { width:275px; }
      .project-session-card .input-group { margin-bottom:0; }
      .project-session-card .btn { margin-bottom:0; }
      .btn-project { background:#f8fafc; border-color:#d7e0ea; color:#334155; }
      #pipeline_step.nav-tabs { background:#fff; padding:5px; border:1px solid var(--ui-border); border-radius:12px; box-shadow:0 2px 8px rgba(15,23,42,.03); gap:3px; margin-bottom:18px; display:flex; flex-wrap:wrap; }
      #pipeline_step.nav-tabs > li > a { border:none; border-radius:8px; padding:9px 12px; font-size:12.5px; color:#64748b; font-weight:650; box-shadow:none !important; }
      #pipeline_step.nav-tabs > li.active > a,
      #pipeline_step.nav-tabs > li.active > a:hover { border:none; color:#1746a2; background:#eef5ff; box-shadow:none !important; }
      .stage-heading { display:grid; grid-template-columns:54px minmax(0,1fr) auto; align-items:center; gap:15px; padding:15px 18px; margin:0 0 12px; background:linear-gradient(135deg,#fff 0%,#fbfdff 100%); border:1px solid var(--ui-border); border-radius:14px; box-shadow:0 3px 12px rgba(15,23,42,.035); }
      .stage-heading-icon { width:48px; height:48px; border-radius:13px; background:#eff6ff; color:#2563eb; display:flex; align-items:center; justify-content:center; font-size:20px; }
      .stage-heading-copy h3 { margin:0 0 2px; color:#14243a; font-size:20px; font-weight:800; letter-spacing:-.015em; }
      .stage-heading-subtitle { color:#687b92; font-size:12.5px; line-height:1.35; }
      .stage-heading-badge { align-self:start; background:#f1f5f9; color:#62748b; border-radius:999px; padding:5px 9px; font-size:11px; font-weight:700; }
      .stage-topbar { min-height:49px; padding:7px 9px; background:#fff; position:sticky; top:0; z-index:20; }
      .stage-topbar-spacer { flex:1 1 auto; }
      .stage-topbar .btn { min-height:34px; padding:6px 11px; }
      .panel-box, .taxonomy-card, .taxonomy-full-card { background:#fff; }
      .panel-box { padding:19px 20px; margin-bottom:16px; }
      .card-title-row { display:flex; align-items:center; gap:8px; min-height:28px; margin:0 0 14px; color:#17263a; }
      .card-title-icon { width:31px; height:31px; border-radius:9px; display:flex; align-items:center; justify-content:center; background:#f1f5f9; color:#52647a; flex:0 0 auto; font-size:13px; }
      .card-title-text { font-size:16px; font-weight:800; letter-spacing:-.01em; }
      .card-title-row .info-tip { margin-left:auto; }
      .info-tip { cursor:help; color:#8aa0b8; font-size:14px; }
      .compact-hint { display:flex; align-items:center; gap:7px; color:#708197; font-size:11.5px; margin-top:8px; }
      .stage-grid { display:grid; gap:16px; align-items:start; }
      .stage-grid-2 { grid-template-columns:minmax(0,1fr) minmax(0,1fr); }
      .stage-grid-upload { grid-template-columns:minmax(300px,.72fr) minmax(0,1.55fr); }
      .stage-grid-rename { grid-template-columns:minmax(310px,.82fr) minmax(0,1.35fr); }
      .upload-drop-card { min-height:260px; }
      .upload-drop-card .shiny-input-container { width:100%; }
      .upload-drop-card .input-group { width:100%; }
      .upload-drop-card .btn-file { background:#eef5ff; border-color:#bdd2f5; color:#1d4ed8; font-weight:700; }
      .stage-table-card { min-width:0; overflow:hidden; }
      .form-grid-2, .form-grid-3 { display:grid; gap:10px 14px; align-items:end; }
      .form-grid-2 { grid-template-columns:repeat(2,minmax(0,1fr)); }
      .form-grid-3 { grid-template-columns:repeat(3,minmax(0,1fr)); }
      .form-grid-2 .form-group, .form-grid-3 .form-group { margin-bottom:6px; width:100%; }
      .settings-card .form-group > label, .blast-query-controls .form-group > label { color:#52647a; font-size:12px; font-weight:700; }
      .inline-control-row { display:flex; gap:8px; align-items:center; margin:5px 0 11px; }
      .inline-control-row .form-group { margin:0; }
      .subsection-divider { border-top:1px solid #e7edf4; margin:14px 0; }
      .subsection-title { color:#334155; font-size:12px; font-weight:750; margin:0 0 8px; display:flex; gap:6px; align-items:center; }
      .qc-inspection-grid { display:grid; grid-template-columns:minmax(290px,.42fr) minmax(0,1.58fr); gap:16px; align-items:start; }
      .qc-main-stack { min-width:0; }
      .qc-sidebar-card { position:sticky; top:62px; }
      .compact-plot-card { padding-bottom:10px; }
      .chromatogram-card { padding-bottom:14px; }
      .peak-review-card .curation-toolbar { background:#f8fafc; border:1px solid #e5ebf2; border-radius:10px; padding:9px 11px; }
      .peak-review-card .curation-toolbar .form-group { margin-bottom:0; }
      .sequence-preview-card textarea { background:#fbfcfe; border-color:#dfe6ef !important; }
      .checkpoint-modern { display:flex; align-items:center; gap:18px; justify-content:space-between; background:#fffdf7; }
      .checkpoint-modern .card-title-row { margin-bottom:0; }
      .checkpoint-copy { flex:1 1 auto; }
      .checkpoint-modern .btn { flex:0 0 auto; margin-bottom:0; }
      .export-summary-card { min-height:105px; }
      .export-tile-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:15px; margin-bottom:18px; }
      .export-tile { min-height:210px; background:#fff; border:1px solid var(--ui-border); border-radius:14px; padding:19px; box-shadow:var(--ui-shadow); display:flex; flex-direction:column; align-items:flex-start; }
      .export-tile-primary { border-color:#bdd2f5; background:linear-gradient(145deg,#fff 0%,#f7fbff 100%); }
      .export-tile-icon { width:42px; height:42px; border-radius:12px; background:#eef5ff; color:#2563eb; display:flex; align-items:center; justify-content:center; margin-bottom:14px; font-size:17px; }
      .export-tile h4 { margin:0 0 6px; font-size:15px; font-weight:800; color:#1e2f45; }
      .export-tile p { color:#718198; font-size:11.5px; line-height:1.45; flex:1 1 auto; margin-bottom:16px; }
      .export-tile .btn { width:100%; margin:0; }
      .blast-query-grid { display:grid; grid-template-columns:minmax(340px,.7fr) minmax(0,1.3fr); gap:18px; align-items:stretch; }
      .blast-query-controls { min-width:0; }
      .sequence-code-panel { background:#f8fafc; border:1px solid #e5ebf2; border-radius:10px; padding:13px 14px; min-width:0; }
      .sequence-code-viewer { background:#fff; border:1px solid #dfe6ef; border-radius:8px; min-height:190px; max-height:280px; overflow-y:auto; padding:11px 12px; margin:0; color:#24364b; font-family:Consolas,'Courier New',monospace; font-size:12px; line-height:1.5; white-space:pre-wrap; overflow-wrap:anywhere; word-break:break-all; }
      .sequence-code-empty { color:#94a3b8; font-family:Aptos,Arial,sans-serif; font-style:italic; }
      .button-row { display:flex; flex-wrap:wrap; gap:7px; align-items:center; }
      .button-row .btn, .blast-primary-actions .btn { margin:0; }
      .blast-primary-actions { padding-bottom:12px; border-bottom:1px solid #e8edf4; margin-bottom:10px; }
      .blast-status-strip { display:grid; grid-template-columns:1fr 1fr; gap:10px; margin-bottom:10px; }
      .blast-status-strip:empty { display:none; }
      .export-inline-actions { margin-top:14px; }
      .pipeline-schema-wrap { opacity:.9; margin-top:25px; }
      .pipeline-schema-title { font-size:10px; }
      .schema-label { font-size:10.5px; }
      .dataTables_wrapper { color:#334155; }
      table.dataTable { width:100% !important; }
      table.dataTable thead th { padding-top:10px !important; padding-bottom:10px !important; }
      table.dataTable tbody td { padding-top:9px !important; padding-bottom:9px !important; }
      .dataTables_wrapper .dataTables_filter input, .dataTables_wrapper .dataTables_length select { border:1px solid #cfd9e6; border-radius:7px; background:#fff; }
      .raw-response-details summary { font-size:11.5px; padding:6px 9px; }
      @media (max-width:1100px) {
        .stage-grid-2, .stage-grid-upload, .stage-grid-rename, .qc-inspection-grid, .blast-query-grid { grid-template-columns:1fr; }
        .qc-sidebar-card { position:static; }
        .export-tile-grid { grid-template-columns:repeat(2,minmax(0,1fr)); }
      }
      @media (max-width:700px) {
        .pipeline-container { padding-left:8px; padding-right:8px; }
        .stage-heading { grid-template-columns:44px 1fr; }
        .stage-heading-icon { width:40px; height:40px; }
        .stage-heading-badge { grid-column:2; justify-self:start; }
        .form-grid-2, .form-grid-3, .blast-status-strip { grid-template-columns:1fr; }
        .export-tile-grid { grid-template-columns:1fr; }
        .checkpoint-modern { align-items:flex-start; flex-direction:column; }
        .project-session-card .shiny-input-container { width:100%; }
      }

      /* Taxonomy dashboard ------------------------------------------------- */
      .taxonomy-workspace { display:flex; gap:18px; align-items:center; justify-content:space-between; padding:15px 18px; border:1px solid var(--ui-border); border-radius:13px; background:#fff; box-shadow:var(--ui-shadow); margin-bottom:18px; }
      .taxonomy-workspace-left { display:flex; align-items:flex-end; gap:8px; flex:1 1 560px; min-width:0; }
      .taxonomy-workspace-left .form-group { margin-bottom:0; width:min(620px,100%); }
      .taxonomy-workspace-status { flex:1 1 380px; min-width:260px; }
      .taxonomy-workspace-status .tax-note { margin:0; border:none; background:#f7f9fc; padding:10px 12px; }

      .taxonomy-result-hero { --accent:#2563eb; --soft:#eff6ff; position:relative; display:grid; grid-template-columns:94px minmax(0,1fr) minmax(300px,420px); gap:22px; background:#fff; border:1px solid var(--ui-border); border-left:5px solid var(--accent); border-radius:15px; box-shadow:var(--ui-shadow); overflow:hidden; margin:0 0 20px; padding:22px 24px 22px 0; }
      .taxonomy-result-hero.conf-high { --accent:#159447; --soft:#ecfdf3; }
      .taxonomy-result-hero.conf-moderate { --accent:#d97706; --soft:#fff7e6; }
      .taxonomy-result-hero.conf-low { --accent:#e17a00; --soft:#fff6e6; }
      .taxonomy-hero-icon-wrap { display:flex; align-items:center; justify-content:center; background:linear-gradient(90deg,var(--soft),rgba(255,255,255,0)); min-height:138px; }
      .taxonomy-hero-icon { width:54px; height:54px; border-radius:50%; display:flex; align-items:center; justify-content:center; background:var(--soft); color:var(--accent); font-size:26px; border:1px solid color-mix(in srgb, var(--accent) 18%, white); }
      .taxonomy-hero-main { min-width:0; padding:3px 0; }
      .taxonomy-hero-title { font-size:30px; line-height:1.15; font-weight:800; letter-spacing:-.025em; color:#14243a; margin-bottom:10px; }
      .taxonomy-hero-title .level { font-size:20px; font-weight:650; color:#46576e; }
      .taxonomy-hero-confidence { font-size:17px; color:var(--accent); font-weight:800; margin-left:8px; }
      .taxonomy-hero-lines { display:grid; gap:7px; color:#334155; font-size:14px; }
      .taxonomy-hero-line { display:flex; align-items:flex-start; gap:9px; line-height:1.35; }
      .taxonomy-hero-line .mini-icon { width:18px; color:#62748c; text-align:center; flex:0 0 18px; }
      .taxonomy-hero-line em { color:#1556c0; }
      .taxonomy-why { align-self:stretch; background:var(--soft); border-radius:10px; padding:15px 17px; color:#334155; line-height:1.45; }
      .taxonomy-why-title { color:color-mix(in srgb, var(--accent) 86%, #111827); font-weight:800; margin-bottom:8px; display:flex; align-items:center; justify-content:space-between; gap:10px; }
      .taxonomy-why-text { font-size:13px; }

      .taxonomy-metric-grid { display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:14px; margin-bottom:20px; }
      .tax-metric { --metric:#2563eb; --metric-soft:#eff6ff; min-height:154px; background:#fff; border:1px solid var(--ui-border); border-radius:13px; box-shadow:0 3px 12px rgba(15,23,42,.04); padding:17px 17px 15px; display:grid; grid-template-columns:48px 1fr; column-gap:12px; align-content:start; }
      .tax-metric.metric-purple { --metric:#6d4bd2; --metric-soft:#f3efff; }
      .tax-metric.metric-blue { --metric:#2563eb; --metric-soft:#eff6ff; }
      .tax-metric.metric-green { --metric:#159447; --metric-soft:#ecfdf3; }
      .tax-metric.metric-cyan { --metric:#159bb7; --metric-soft:#ecfbff; }
      .tax-metric.metric-amber { --metric:#d97706; --metric-soft:#fff7e6; }
      .tax-metric-icon { width:44px; height:44px; border-radius:50%; display:flex; align-items:center; justify-content:center; background:var(--metric-soft); color:var(--metric); font-size:20px; grid-row:1 / span 4; }
      .tax-metric-label { font-size:13px; line-height:1.32; color:#25364c; font-weight:650; min-height:34px; }
      .tax-metric-value { font-size:29px; line-height:1.1; font-weight:800; color:#16263b; margin-top:5px; letter-spacing:-.02em; }
      .tax-metric-value.metric-colored { color:var(--metric); }
      .tax-metric-sub { color:#64748b; font-size:12px; line-height:1.35; margin-top:6px; }
      .tax-meter { height:5px; background:#e6ebf1; border-radius:999px; overflow:hidden; margin-top:11px; }
      .tax-meter-fill { height:100%; background:var(--metric); border-radius:999px; }

      .taxonomy-main-grid { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1.12fr); gap:18px; align-items:start; margin-bottom:18px; }
      .taxonomy-card { background:#fff; border:1px solid var(--ui-border); border-radius:14px; box-shadow:var(--ui-shadow); padding:20px; min-width:0; }
      .taxonomy-card-title { font-size:17px; color:#17263a; font-weight:800; margin:0 0 14px; display:flex; align-items:center; gap:7px; }
      .taxonomy-card .dataTables_wrapper { font-size:12px; }
      .taxonomy-card .dataTables_wrapper table.dataTable { width:100% !important; }
      .tax-decision-table tbody tr { background:#f2fbf5 !important; }
      .tax-decision-table tbody td { font-weight:600; }
      .tax-callout { display:grid; grid-template-columns:24px 1fr; gap:9px; background:#f5f9ff; border:1px solid #e1ecfb; border-radius:9px; padding:11px 12px; margin-top:9px; color:#334155; font-size:12px; line-height:1.45; }
      .tax-callout-icon { color:#2563eb; font-size:15px; padding-top:1px; }
      .tax-callout strong { color:#1d4ed8; display:block; margin-bottom:2px; }
      .taxonomy-agreement-note { display:grid; grid-template-columns:22px 1fr; gap:8px; background:#f5f9ff; border-radius:8px; padding:10px 11px; margin-top:10px; color:#475569; font-size:11.5px; line-height:1.4; }
      .taxonomy-full-card { background:#fff; border:1px solid var(--ui-border); border-radius:14px; box-shadow:var(--ui-shadow); padding:20px 22px; margin-bottom:18px; }
      .taxonomy-full-card .section-title { margin-bottom:10px; }
      .taxonomy-explain-details { margin:0 0 10px; }
      .taxonomy-explain-details summary { cursor:pointer; color:#4b5f79; font-size:12px; font-weight:650; display:inline-flex; align-items:center; gap:6px; padding:6px 9px; border-radius:7px; background:#f7f9fc; border:1px solid #e2e8f0; }
      .taxonomy-explain-details[open] summary { margin-bottom:8px; }
      .taxonomy-explain-body { color:#52647a; font-size:12px; line-height:1.5; background:#fafcff; border-radius:8px; padding:10px 12px; }

      @media (max-width:1180px) {
        .taxonomy-metric-grid { grid-template-columns:repeat(3,minmax(0,1fr)); }
        .taxonomy-result-hero { grid-template-columns:82px minmax(0,1fr); padding-right:20px; }
        .taxonomy-why { grid-column:1 / -1; margin-left:20px; }
      }
      @media (max-width:900px) {
        .taxonomy-main-grid { grid-template-columns:1fr; }
        .taxonomy-metric-grid { grid-template-columns:repeat(2,minmax(0,1fr)); }
        .taxonomy-workspace { align-items:stretch; flex-direction:column; }
        .taxonomy-workspace-status { min-width:0; }
      }
      @media (max-width:620px) {
        .app-shell-header { min-height:88px; padding:14px 16px; gap:10px; }
        .app-brand-logo { width:min(300px,70vw); max-height:66px; }
        .app-version-badge { margin-left:auto; }
        .taxonomy-result-hero { grid-template-columns:1fr; padding:18px; }
        .taxonomy-hero-icon-wrap { min-height:0; justify-content:flex-start; background:none; }
        .taxonomy-why { grid-column:auto; margin-left:0; }
        .taxonomy-metric-grid { grid-template-columns:1fr; }
      }
    ")),
    tags$script(HTML("
      // The source keeps the large QC panel before Rename for maintainability,
      // while the user-facing workflow is Rename -> QC.
      function pitaxOrderWorkflowTabs() {
        var nav = $('#pipeline_step');
        var renameTab = nav.find('a[data-value=\"rename\"]').parent();
        var qcTab = nav.find('a[data-value=\"qc\"]').parent();
        if (renameTab.length && qcTab.length) renameTab.insertBefore(qcTab);
      }
      $(pitaxOrderWorkflowTabs);
      $(document).on('shiny:connected', pitaxOrderWorkflowTabs);
    "))
  ),

  div(id = "shiny_activity_bar", `aria-hidden` = "true"),
  div(id = "app_loading_overlay", class = "app-loading-overlay", `aria-hidden` = "true",
      div(class = "app-loading-card",
          div(class = "app-loading-spinner"),
          div(id = "app_loading_text", "Loading workspace…")
      )
  ),

  div(class = "pipeline-container",
    div(class = "app-header app-shell-header",
      if (PITAX_LOGO_AVAILABLE)
        tags$img(src = "pitax-assets/logo.png", class = "app-brand-logo", alt = "PITAX — Taxonomic Identification Tool")
      else
        tagList(
          div(class = "app-brand-mark", icon("flask")),
          div(class = "app-brand-copy", h2("PITAX"), p("Taxonomic Identification Tool"))
        ),
      div(class = "app-brand-copy",
        p("Sanger sequence analysis · quality review · BLAST · taxonomic interpretation")
      ),
      div(class = "app-version-badge", paste0("v", APP_VERSION))
    ),

    div(class = "project-bar project-session-card",
      div(class = "project-bar-title", icon("folder-open"), span("Project session")),
      downloadButton("save_project", "Save project", class = "btn-project"),
      fileInput("load_project", NULL, multiple = FALSE, accept = c(".sangerproject", ".rds"),
                buttonLabel = "Load project", placeholder = "No project selected"),
      div(class = "project-status", uiOutput("project_status"))
    ),

    tabsetPanel(id = "pipeline_step", type = "tabs",

      # --------------------------------------------------------
      # 1. Upload
      # --------------------------------------------------------
      tabPanel("1 · Upload", value = "upload",
        stage_heading("upload", "Upload and identify chromatograms", "Keep the sequencer barcode as the immutable source, then assign isolate, gene and direction before processing.", "Step 1 of 7"),
        stage_topbar(
          div(class = "stage-topbar-spacer"),
          actionButton("to_settings", "Continue to Assay", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "stage-grid stage-grid-upload",
          div(class = "panel-box upload-drop-card",
            card_title("Raw AB1 files", "Select one or more .ab1 chromatogram files. The original files remain the immutable source for all later QC and curation.", "upload"),
            fileInput("ab1_files", NULL, multiple = TRUE,
                      accept = c(".ab1", ".AB1"), buttonLabel = "Select AB1 files", placeholder = "or drag files here"),
            div(class = "compact-hint", icon("info-circle"), span("Multiple files can be processed in one run."))
          ),
          div(class = "panel-box stage-table-card",
            card_title("Files in this run", "Review the uploaded filenames before continuing.", "list"),
            DTOutput("uploaded_files_table")
          )
        ),
        div(class = "stage-grid stage-grid-rename",
          div(class = "panel-box",
            card_title("Assignment key", "Import XLSX/CSV with old_id, isolate, locus (or gene), and direction columns. old_id matches the uploaded barcode; prefix matching is supported.", "key"),
            fileInput("rename_key_file", NULL, multiple = FALSE, accept = c(".xlsx", ".csv"), buttonLabel = "Choose assignment key", placeholder = "XLSX or CSV"),
            actionButton("apply_rename_key", "Apply assignment key", icon = icon("check"), class = "btn-primary"),
            downloadButton("download_assignment_key_template", "Download key template"),
            uiOutput("rename_key_status")
          ),
          div(class = "panel-box",
            card_title("Batch assignment", "Apply isolate edits, gene and direction to selected rows. If no rows are selected, the action applies to all uploaded reads.", "edit"),
            div(class = "form-grid-2",
              textInput("batch_isolate_prefix", "Isolate prefix", ""),
              textInput("batch_isolate_suffix", "Isolate suffix", "")
            ),
            div(class = "form-grid-2",
              textInput("batch_isolate_find", "Find in isolate", ""),
              textInput("batch_isolate_replace", "Replace with", "")
            ),
            div(class = "form-grid-2",
              textInput("batch_locus", "Set gene / locus", ""),
              selectInput("batch_direction", "Set direction", c("No change" = "", "Forward" = "Forward", "Reverse" = "Reverse"), selected = "")
            ),
            actionButton("apply_assignment_batch", "Apply to selected / all", class = "btn-primary"),
            actionButton("reset_assignments", "Clear assignments")
          )
        ),
        div(class = "panel-box stage-table-card",
          card_title("Read identity and generated FASTA names", "Edit Isolate, Gene / locus and Direction directly. PITAX generates the final name; it never extracts biological identity from the upload barcode.", "list-alt"),
          DTOutput("assignment_upload_table"),
          uiOutput("upload_assignment_validation"),
          uiOutput("upload_architecture_summary")
        ),
        pipeline_stage_footer(1)
      ),
      # --------------------------------------------------------
      # 2. Assay and trimming settings
      # --------------------------------------------------------
      tabPanel("2 · Assay", value = "settings",
        stage_heading("sliders", "Assay setup", "Define run-level locus/primer defaults and the automatic trimming rules. Trimming starts only after Rename.", "Step 2 of 7"),
        stage_topbar(
          actionButton("back_upload", "Back", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_rename", "Continue to Rename", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "stage-grid stage-grid-2",
          div(class = "panel-box settings-card",
            card_title("Assay information", "Primer metadata is retained for provenance. Experimental primer-site mapping is optional and does not control trimming.", "flask"),
            selectInput("target", "Target / Gene",
                        c("ITS", "LSU", "TEF1 / EF1-alpha", "RPB2", "Beta-tubulin", "CYP51", "SDHB", "IGS", "Other"),
                        selected = "ITS"),
            div(class = "form-grid-2",
              textInput("forward_primer", "Forward primer name", ""),
              textInput("reverse_primer", "Reverse primer name", "")
            ),
            div(class = "form-grid-2",
              textInput("forward_primer_seq", "Forward primer sequence (5'→3')", ""),
              textInput("reverse_primer_seq", "Reverse primer sequence (5'→3')", "")
            ),
            radioButtons(
              "sequencing_primer", "Primer used for this Sanger read",
              choices = c("Forward" = "Forward", "Reverse" = "Reverse", "Unknown / infer" = "Unknown"),
              selected = "Forward", inline = TRUE
            ),
            div(class = "inline-control-row",
              checkboxInput("enable_primer_mapping", "Experimental primer-site mapping", value = FALSE),
              info_tip("Primer-site inference can be unreliable because the sequencing primer itself may not appear in the base calls and the beginning of a Sanger read is often noisy.")
            ),
            div(class = "form-grid-2",
              numericInput("expected_amplicon_len", "Expected amplicon length (bp)", 650, min = 1),
              numericInput("absolute_max_base_index", "Maximum sequence position (bp)", 680, min = 50)
            )
          ),
          div(class = "panel-box settings-card",
            card_title("Automatic trimming", "These parameters define good-start detection and sustained signal-collapse detection. Full definitions are available in Help / About.", "cut"),
            div(class = "form-grid-2",
              numericInput("window", "Window size", 25, min = 5),
              numericInput("min_peak_ratio", "Minimum peak ratio", 3, min = 0, step = 0.1),
              numericInput("min_relative_signal", "Minimum relative signal", 0.20, min = 0, max = 1, step = 0.01),
              numericInput("min_len_before_collapse", "Minimum length before collapse search", 350, min = 1),
              numericInput("bad_run_windows", "Consecutive bad windows", 12, min = 1),
              numericInput("min_usable_len", "Minimum usable trimmed length", 400, min = 1)
            )
          )
        ),
        pipeline_stage_footer(2)
      ),
      # --------------------------------------------------------
      # 4. QC, chromatogram & sequence preview
      # --------------------------------------------------------
      tabPanel("4 · Trim & QC", value = "qc",
        stage_heading("bar-chart", "Trimming results, QC & curation", "Review the completed trim, inspect renamed chromatograms, and document manual sequence curation.", "Step 4 of 7"),
        stage_topbar(
          actionButton("back_rename_from_qc", "Back to Rename", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_export", "Continue to Export", icon = icon("arrow-right"), class = "btn-primary")
        ),
        uiOutput("qc_summary_cards"),
        div(class = "panel-box stage-table-card",
          card_title("Trimming results", "Run-level summary of automatic trimming for all uploaded chromatograms.", "table"),
          DTOutput("summary_table")
        ),
        div(class = "qc-inspection-grid",
          div(class = "panel-box qc-sidebar-card",
            card_title("Sequence inspection", "Select a sample to inspect its current curated sequence, QC metrics and optional primer mapping.", "search"),
            selectInput("inspect_sample", "Sample", choices = NULL),
            DTOutput("sequence_metrics"),
            conditionalPanel(
              condition = "input.enable_primer_mapping === true",
              div(class = "subsection-divider"),
              div(class = "subsection-title", "Primer-site mapping", info_tip("Experimental mapping only; it is not used as an automatic trim boundary.")),
              DTOutput("primer_match_table")
            )
          ),
          div(class = "qc-main-stack",
            div(class = "panel-box compact-plot-card",
              card_title("Read overview", "Automatic trim boundaries and optional primer context.", "arrows-h"),
              plotOutput("amplicon_overview", height = "190px"),
              uiOutput("expected_amplicon_note")
            ),
            div(class = "panel-box chromatogram-card",
              card_title("Chromatogram", "X axis = raw called-base position. Zoom, pan, or click a flagged position to inspect it.", "line-chart"),
              plotly::plotlyOutput("chromatogram_plot", height = "540px")
            )
          )
        ),
        div(class = "panel-box stage1-evidence-card",
          card_title("Stage 1 · AB1 evidence audit", "Observational audit of the established PITAX v2 read model against raw ABIF primary-call coordinates, basecaller quality and canonical A/C/G/T trace evidence. This panel does not alter trimming, curation, BLAST or taxonomy.", "microscope"),
          tags$details(
            class = "stage1-audit-details",
            tags$summary(class = "stage1-audit-toggle", "Open evidence audit"),
            div(class = "status-note",
                tags$strong("Validation mode only. "),
                "The active v2.14.2 trimming/QC path is still the decision path. alpha.6 keeps the same-length PCON-only comparison observational; it does not apply that proposal to the processed sequence, FASTA or BLAST output."),
            div(class = "subsection-title", "Run-level audit"),
            DTOutput("ab1_evidence_run_table"),
            div(class = "blast-action-row",
              downloadButton("download_ab1_evidence_run", "Download run audit CSV", class = "btn-default"),
              downloadButton("download_ab1_evidence_detail", "Download selected-base audit CSV", class = "btn-default")
            ),
            div(class = "subsection-divider"),
            div(class = "subsection-title", "Selected sample · per-base evidence"),
            uiOutput("ab1_evidence_selected_note"),
            div(class = "subsection-title", "Legacy auto trim vs PCON-only comparison"),
            DTOutput("ab1_trim_comparison_table"),
            DTOutput("ab1_evidence_detail_table")
          )
        ),
        div(class = "panel-box peak-review-card",
          card_title("Ambiguous peak review", "Flags identify channel competition or locally non-dominant calls. Left-click centers the trace; right-click opens curation actions.", "flag"),
          div(class = "curation-toolbar",
            radioButtons(
              "peak_flag_scope", "Positions",
              choices = c("Trimmed sequence" = "trimmed", "Entire raw read" = "raw"),
              selected = "trimmed", inline = TRUE
            ),
            checkboxInput("show_peak_flags", "Show markers", value = TRUE),
            checkboxInput("show_reviewed_flags", "Show reviewed", value = FALSE),
            div(class = "curation-actions",
              actionButton("auto_correct_preview", "Auto-correct high-confidence", icon = icon("magic"), title = "Preview positions that meet the current auto-correction criteria. No base is changed before confirmation."),
              actionButton("auto_correct_settings", "Criteria", icon = icon("cog"), title = "Edit the thresholds used to propose automatic base corrections."),
              actionButton("curation_undo", "Undo", icon = icon("undo"), title = "Undo the most recent curation transaction for this sample."),
              actionButton("curation_redo", "Redo", icon = icon("repeat"), title = "Redo the most recently undone curation transaction."),
              actionButton("curation_history", "History", icon = icon("history"), title = "Open the full per-sample manual curation audit log."),
              actionButton("curation_reset", "Reset to auto", icon = icon("refresh"), title = "Restore automatic trim/base calls as the active curated sequence. The reset is logged and undoable.")
            )
          ),
          uiOutput("ambiguous_peak_summary"),
          DTOutput("ambiguous_peak_table")
        ),
        div(class = "panel-box",
          card_title("QC signal metrics", "Called-base signal and called-peak / second-peak ratio. These plots are retained in QC exports.", "area-chart"),
          plotOutput("qc_plot", height = "650px")
        ),
        conditionalPanel(
          condition = "input.enable_primer_mapping === true",
          div(class = "panel-box",
            card_title("Experimental primer alignment", "Detailed primer alignment output for manual inspection.", "code"),
            verbatimTextOutput("primer_alignment")
          )
        ),
        div(class = "panel-box sequence-preview-card",
          card_title("Curated sequence", "This exact sequence retains its resolved sample name and is carried forward to export and BLAST.", "file-text"),
          textAreaInput("trimmed_sequence_preview", NULL, value = "", rows = 8, width = "100%")
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy",
            card_title("Checkpoint B · Renamed and curated sequences", "Save the resolved naming and QC/curation state before export.", "save")
          ),
          downloadButton("download_trim_checkpoint", "Download checkpoint ZIP")
        ),
        pipeline_stage_footer(4)
      ),
      # --------------------------------------------------------
      # 3. Rename
      # --------------------------------------------------------
      tabPanel("3 · Rename & Assign", value = "rename",
        stage_heading("tags", "Review read identity before trimming", "Confirm the explicit isolate, gene and direction fields entered at Upload. The final FASTA name is generated from those fields.", "Step 3 of 7"),
        stage_topbar(
          actionButton("back_settings_from_rename", "Back to Assay", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("run_trimming", "Start trimming", icon = icon("play"), class = "btn-primary")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Final read / FASTA names", "The source barcode stays unchanged. Edit the three identity fields; PITAX generates <Isolate>_<Locus>_<F/R>.", "list-alt"),
          DTOutput("assignment_review_table"),
          uiOutput("rename_validation")
        ),
        div(class = "panel-box stage2-assignment-card",
          card_title("Stage 2 · Architecture preview", "Built directly from the explicit identity fields. A pair requires two distinct source reads assigned to the same isolate/locus, one Forward and one Reverse.", "sitemap"),
          uiOutput("architecture_summary")
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint A · Renamed and assigned reads", "Save resolved read names and biological identity before trimming.", "save")),
          downloadButton("download_rename_checkpoint", "Download checkpoint ZIP")
        ),
        pipeline_stage_footer(3)
      ),
      # --------------------------------------------------------
      # 5. Export
      # --------------------------------------------------------
      tabPanel("5 · Export", value = "export",
        stage_heading("download", "Export processed sequences", "Create working FASTA files or an auditable package of the complete processing run.", "Step 5 of 7"),
        stage_topbar(
          actionButton("back_qc_from_export", "Back to QC", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_blast", "Continue to NCBI BLAST", icon = icon("arrow-right"), class = "btn-primary")
        ),
        div(class = "panel-box export-summary-card",
          card_title("Run summary", "Final pre-export overview of the processed sequence set.", "check-circle"),
          uiOutput("export_summary")
        ),
        div(class = "export-tile-grid",
          div(class = "export-tile", div(class = "export-tile-icon", icon("search")), h4("BLAST FASTA"), p("Processed sequences ready for sequence search."), downloadButton("download_blast_fasta", "Download FASTA")),
          div(class = "export-tile", div(class = "export-tile-icon", icon("file-text")), h4("FASTA + metadata"), p("Processed sequences with run metadata in FASTA headers."), downloadButton("download_full_fasta", "Download FASTA")),
          div(class = "export-tile", div(class = "export-tile-icon", icon("table")), h4("Summary CSV"), p("Compact processing and QC summary for downstream review."), downloadButton("download_summary_csv", "Download CSV")),
          div(class = "export-tile export-tile-primary", div(class = "export-tile-icon", icon("archive")), h4("Complete results package"), p("Sequences, QC evidence, settings and curation records in one ZIP."), downloadButton("download_all_zip", "Download results ZIP"))
        ),
        pipeline_stage_footer(5)
      ),
      # --------------------------------------------------------
      # 6. NCBI BLAST
      # --------------------------------------------------------
      tabPanel("6 · NCBI BLAST", value = "blast",
        stage_heading("search", "NCBI BLAST workspace", "Submit curated sequences, retrieve accession-level hits, and keep each RID linked to its sequence revision.", "Step 6 of 7"),
        stage_topbar(
          actionButton("back_export", "Back to Export", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("to_taxonomy", "Continue to Taxonomic Summary", icon = icon("arrow-right"), class = "btn-primary"),
          actionButton("reset_pipeline", "Start new run")
        ),
        div(class = "panel-box blast-query-card",
          card_title("Query workspace", "The sequence shown here is the current curated and renamed sequence. Changing a curated sequence after BLAST marks its previous result stale.", "file-code-o"),
          div(class = "blast-query-grid",
            div(class = "blast-query-controls",
              selectInput("blast_sample", "Sequence", choices = NULL),
              div(class = "form-grid-2",
                selectInput("blast_database", "NCBI database", choices = c("core_nt", "nt"), selected = "core_nt"),
                numericInput("blast_hitlist", tagList("Maximum hits", info_tip("Controls how many BLAST hits NCBI returns. Taxonomic interpretation uses accession-level competitive evidence rather than flat Top-N voting.")), 25, min = 1, max = 100)
              ),
              div(class = "button-row",
                actionButton("copy_blast_sequence", "Copy sequence", icon = icon("copy")),
                actionButton("open_ncbi_blast", "Open NCBI BLAST", icon = icon("external-link")),
                downloadButton("download_selected_blast", "Selected FASTA")
              )
            ),
            div(class = "sequence-code-panel",
              div(class = "subsection-title", "Processed query sequence"),
              uiOutput("blast_sequence_preview_ui")
            )
          )
        ),
        div(class = "panel-box blast-submit-card",
          card_title("Automated submission", "Each sample receives its own NCBI Request ID (RID). Automated contacts are paced to respect NCBI service limits.", "cloud-upload"),
          div(class = "blast-action-row blast-primary-actions",
            actionButton("submit_ncbi_blast", "Submit selected", icon = icon("paper-plane"), class = "btn-primary"),
            actionButton("submit_all_ncbi_blast", "Submit all", icon = icon("paper-plane"), class = "btn-success"),
            actionButton("retrieve_ncbi_blast", "Retrieve selected", icon = icon("refresh")),
            actionButton("retrieve_all_ncbi_blast", "Retrieve all", icon = icon("download"))
          ),
          div(class = "blast-status-strip", uiOutput("blast_batch_status"), uiOutput("blast_job_status")),
          DTOutput("blast_jobs_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Preliminary top-hit overview", "One row per retrieved sequence. This is a quick orientation only; multi-hit evidence is interpreted in the next stage.", "eye"),
          DTOutput("blast_identification_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Retrieved accession-level hits", "All parsed unique accessions for the selected RID. Multiple HSPs belonging to the same accession are aggregated.", "database"),
          DTOutput("blast_hits_table"),
          tags$details(
            class = "raw-response-details",
            tags$summary("Raw NCBI response"),
            verbatimTextOutput("blast_raw_preview")
          )
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint C · BLAST workspace", "Export BLAST job metadata and accession-level hits.", "save")),
          div(class = "button-row",
            downloadButton("download_blast_jobs", "Job/results CSV"),
            downloadButton("download_blast_hits", "All BLAST hits CSV")
          )
        ),
        pipeline_stage_footer(6)
      ),
      # --------------------------------------------------------
      # 7. Taxonomic interpretation
      # --------------------------------------------------------
      tabPanel("7 · Taxonomic summary", value = "taxonomy",
        stage_heading("sitemap", "Taxonomic interpretation", "Identify the best molecular match, inspect close alternatives and report the most conservative supported taxonomic level.", "Step 7 of 7"),
        stage_topbar(
          actionButton("back_blast", "Back to NCBI BLAST", icon = icon("arrow-left")),
          div(class = "stage-topbar-spacer"),
          actionButton("run_taxonomy", "Analyze selected", icon = icon("play"), class = "btn-primary"),
          actionButton("run_taxonomy_all", "Analyze all retrieved", icon = icon("tasks"), class = "btn-success"),
          actionButton("reset_pipeline_tax", "Start new run")
        ),
        div(class = "taxonomy-workspace",
          div(class = "taxonomy-workspace-left",
            selectInput("tax_sample", "Sequence", choices = NULL),
            info_tip("Only sequences with retrieved BLAST hits appear here. Each NCBI accession is counted once. The decision starts from Identity + query coverage, then checks close alternative taxa, sequence evidence and reference context.")
          ),
          div(class = "taxonomy-workspace-status", uiOutput("taxonomy_status"))
        ),
        uiOutput("taxonomy_hero"),
        uiOutput("taxonomy_summary_cards"),
        div(class = "taxonomy-main-grid",
          div(class = "taxonomy-card",
            div(class = "taxonomy-card-title", "Taxonomic interpretation", info_tip("The compact decision row shows the final recommendation. Detailed decision components remain available in exports and Help / About.")),
            DTOutput("taxonomy_summary_table"),
            uiOutput("taxonomy_locus_note")
          ),
          div(class = "taxonomy-card",
            div(class = "taxonomy-card-title", "Species evidence profile", info_tip("One row per resolved species. Best Identity and coverage describe that species' strongest comparable accession; accession count is database context only, not a vote.")),
            DTOutput("taxonomy_counts_table"),
            div(class = "taxonomy-agreement-note",
              div(class = "tax-callout-icon", icon("info-circle")),
              div("Best molecular match is chosen from the near-best query-coverage band and then by Identity. Close alternatives are species with nearly the same Identity and coverage; database abundance is shown separately and does not decide the identification.")
            )
          )
        ),
        div(class = "taxonomy-full-card",
          card_title("BLAST score landscape", "Each point is one unique NCBI accession. Points are colored by genus when several genera are present, or by species when the uncertainty is within one genus.", "line-chart"),
          tags$details(class = "taxonomy-explain-details",
            tags$summary(icon("info-circle"), "How to read this graph"),
            div(class = "taxonomy-explain-body",
              "Read left to right by BLAST score rank. The connecting line shows the score landscape, while point colors show which taxa occupy that landscape. When several genera occur, colors represent genus; when all resolved hits are within one genus, colors represent species. Hover to compare Identity, query coverage and accession. Similar colors at similar heights can indicate repeated support, but the number of accessions is not treated as a majority vote."
            )
          ),
          plotly::plotlyOutput("taxonomy_score_plot", height = "390px")
        ),
        div(class = "taxonomy-full-card",
          card_title("Taxonomy-enriched BLAST hits", "Inspect accession-level evidence and the taxonomy attached to each usable BLAST hit.", "database"),
          DTOutput("taxonomy_hits_table")
        ),
        div(class = "panel-box stage-table-card",
          card_title("Team identification summary", "One row per processed sequence, combining QC, BLAST and final taxonomic interpretation.", "users"),
          DTOutput("team_summary_table"),
          div(class = "button-row export-inline-actions",
            downloadButton("download_team_summary_csv", "Team summary CSV"),
            downloadButton("download_team_summary_xlsx", "Team summary Excel")
          )
        ),
        div(class = "panel-box checkpoint checkpoint-modern",
          div(class = "checkpoint-copy", card_title("Checkpoint D · Taxonomic interpretation", "Save the final taxonomic evidence and interpretation package.", "save")),
          div(class = "button-row",
            downloadButton("download_taxonomy_summary", "Summary CSV"),
            downloadButton("download_taxonomy_hits", "Enriched hits CSV"),
            downloadButton("download_taxonomy_checkpoint", "Checkpoint ZIP")
          )
        ),
        pipeline_stage_footer(7)
      ),
      # --------------------------------------------------------
      # Help / About (documentation, not a pipeline stage)
      # --------------------------------------------------------
      tabPanel("Help / About", value = "help",
        panel_box(
          section_title(paste0("Sanger Sequence Pipeline v", APP_VERSION)),
          p(class="about-lead",
            "Documentation for the laboratory workflow, the BLAST/taxonomy interpretation logic, and the scientific sources used to guide the application. Published evidence and application-specific heuristics are labeled separately."),
          div(class="help-flow",
              "AB1 upload  →  Assay settings  →  Rename & read assignment  →  Start trimming  →  QC  →  Export  →  NCBI BLAST  →  Taxonomic interpretation")
        ),

        div(class="about-section",
          tabsetPanel(
            id="about_tabs",
            type="tabs",

            tabPanel("Overview",
              div(class="help-card",
                h3("What the application does"),
                p("The application converts raw Sanger AB1 chromatograms into auditable processed sequences and then supports sequence-based identification using NCBI BLAST and multi-hit taxonomic interpretation."),
                div(class="about-callout",
                  strong("Core principle: "),
                  "the application separates what the database hits agree on from how strong the sequence evidence is. A weak species-level result can therefore fall back to a supported genus instead of becoming automatically Unresolved."
                )
              ),
              fluidRow(
                column(6,
                  div(class="help-card",
                    h3("Workflow"),
                    p(strong("1. Upload"), " — raw AB1 chromatograms plus explicit isolate, gene and direction assignment; source barcodes remain unchanged."),
                    p(strong("2. Assay"), " — run-level primer defaults and trimming parameters; no trimming starts yet."),
                    p(strong("3. Rename & Assign"), " — review the explicit identity fields and PITAX-generated <Isolate>_<Locus>_<F/R> names."),
                    p(strong("4. Trim & QC"), " — start trimming explicitly, then review Quality Control plots, chromatograms and processed sequences."),
                    p(strong("5. Export"), " — processed FASTA, CSV/Excel summaries and checkpoints."),
                    p(strong("6. NCBI BLAST"), " — submit processed sequences and retrieve accession-level hits."),
                    p(strong("7. Taxonomic summary"), " — compare competitive hits and report identification plus confidence.")
                  )
                ),
                column(6,
                  div(class="help-card",
                    h3("Auditability"),
                    p("Read assignments, isolate/locus/read links, QC plots, automatic trimming parameters, manual curation audit log, BLAST hits, taxonomy-enriched hits, decision components and application version are retained in exports."),
                    p(strong("Save / Load project"), " stores the workflow state in a .sangerproject file so completed processing and retrieved BLAST results do not need to be repeated."),
                    p("The final identification is decision-support. It does not replace locus-specific taxonomy, type/reference inspection, morphology or additional loci when those are required.")
                  )
                )
              ),
              div(class="help-card",
                h3("Stage 2 project architecture"),
                p("PITAX stores Project → Isolate → Locus → Read as explicit linked objects. One isolate may have several loci, and one locus may have a single read or separate Forward and Reverse reads."),
                p("The upload barcode is an immutable technical source ID and is never parsed as biological identity. Isolate, Gene/Locus and Direction are edited as separate fields on Upload and reviewed before trimming; PITAX generates the final <Isolate>_<Locus>_<F/R> label from those fields. Duplicate AB1 basenames are blocked because they cannot be represented safely by the established source-read key."),
                div(class="about-callout",
                  strong("Stage boundary: "),
                  "alpha.6 records pairing metadata but does not align or merge Forward/Reverse reads. Consensus construction remains Stage 3."
                )
              )
            ),

            tabPanel("Trimming & QC",
              div(class="help-card",
                h3("Trimming terminology"),
                div(class="method-step",
                  div(class="method-num", "1"),
                  div(h4("Peak ratio"), p("Ratio between the signal of the called base and the second-highest channel at that base call. Larger values indicate a clearer dominant peak."))
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(h4("Good start"), p("The first sustained window that satisfies the configured signal/peak-ratio rule. Bases before this point are excluded from the processed sequence."))
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(h4("Collapse"), p("The first sustained region where relative signal or peak-ratio behavior indicates deterioration. The processed read is normally ended before this point."))
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(h4("Minimum usable length"), p("A workflow threshold used to flag a processed sequence as SHORT_AFTER_TRIMMING rather than OK."))
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(h4("Ambiguous peak review"),
                      p("The app flags individual called-base positions where another A/C/G/T dye channel competes strongly with the current call, where the current call is not locally dominant, or where the call is ambiguous."),
                      p("The comparison uses a narrow trace window around each called peak and suppresses weak-signal competition using a run-specific signal baseline. By default, a called/competitor peak ratio ≤1.25 is Strong and ≤1.75 is Moderate; the local signal must be at least 20% of the retained-region median called signal. These are application review heuristics, not biological diagnoses."),
                      p("Flags are review targets. Left-click centers the chromatogram; right-click opens manual curation actions. No sequence-changing action is applied without an explicit confirmation step."))
                ),
                div(class="method-step",
                  div(class="method-num", "6"),
                  div(h4("Manual curation and provenance"),
                      p("The raw AB1 trace and original automatic base calls are never overwritten. The app maintains a curated sequence layer on top of the automatic trim. A user can change a base, set a two-base IUPAC ambiguity code, trim from either side through a flagged position, or mark the call as reviewed and unchanged."),
                      p("Every confirmed action is written to the Manual Curation audit log with sample, raw base position, before/after values, method, chromatogram evidence, timestamp, transaction ID and revision. Undo and redo are logged as well. Automatic trim boundaries are retained separately so the curated result can always be compared with the automatic result."),
                      p("For a flagged position inside the retained read, the curation menu shows both left- and right-trim actions. A directional recommendation is displayed only when the flag is near a current trim edge: within 15% of the retained length, bounded to 8–40 bases. The recommendation never performs a trim by itself and still requires explicit confirmation."))
                ),
                div(class="method-step",
                  div(class="method-num", "7"),
                  div(h4("High-confidence bulk correction"),
                      p("The bulk correction tool always shows a preview before applying changes. By default, a position is proposed only when the current called base is not locally dominant, the alternative channel is the strongest local channel, alternative/current signal is at least 1.80, alternative/third-channel signal is at least 2.00, the alternative peak maximum lies within ±2 trace samples of the called-base peak position, and the alternative signal is at least 50% of the retained-region median called signal."),
                      p("These thresholds are application heuristics and can be edited from QC & Chromatogram using the Criteria button next to Auto-correct. The active values are saved in the project/run settings and therefore remain part of the analysis provenance. The complete batch is applied as one undoable transaction after user confirmation. Ambiguous double peaks are not bulk-corrected merely because two channels are similar."))
                )
              ),
              div(class="help-card",
                h3("QC evidence"),
                p(strong("QC"), " means Quality Control. Called-base signal and peak-ratio plots are retained in checkpoint/final ZIP exports."),
                p("The ambiguous-peak table links directly to the interactive chromatogram. Right-clicking a row opens the curation actions, while the History view exposes the complete per-sample audit trail."),
                p("If a curated sequence is changed after BLAST has already been run, the previous BLAST/taxonomic result for that sample is marked stale and removed from active interpretation until BLAST is rerun on the new curated sequence.")
              )
            ),

            tabPanel("NCBI BLAST",
              div(class="help-card",
                h3("How BLAST evidence is represented"),
                div(class="method-step",
                  div(class="method-num", "1"),
                  div(h4("RID — Request ID"), p("Identifier returned by NCBI for an individual BLAST job. The app keeps Sample ↔ RID ↔ retrieved hits linked."))
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(h4("Identity"), p("Percentage of identical positions in the representative alignment for the accession."))
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(h4("Query coverage"), p("Percentage of the processed query covered by the subject hit. Multiple local segments from the same accession are combined rather than counted as separate database hits."))
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(h4("E-value and Bit score"), p("E-value describes the expected chance occurrence of a match of similar quality; Bit score is a normalized alignment score used here to compare the leading taxon with competing hits."))
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(h4("HSP — High-scoring Segment Pair"), p("A local alignment segment within a BLAST hit. Multiple HSPs for the same accession are aggregated so that one NCBI record contributes one accession-level hit."))
                )
              ),
              div(class="about-callout",
                strong("Why accession-level aggregation matters: "),
                "without aggregation, one NCBI record containing several HSPs could receive several votes and distort both query-coverage statistics and taxonomic consensus."
              )
            ),

            tabPanel("Taxonomic algorithm",
              div(class="help-card",
                h3("Evidence-first decision logic"),
                p(class="about-lead",
                  "The v2.14 engine uses one decision tree. It starts from the closest molecular match, checks whether other named taxa are practically indistinguishable, and then reports the lowest taxonomic rank that remains defensible. Database abundance is kept as context rather than used as a vote."),

                div(class="method-step",
                  div(class="method-num", "1"),
                  div(
                    h4("One accession = one hit", span(class="evidence-badge evidence-published", "BLAST structure")),
                    p("Multiple HSPs from the same NCBI accession are aggregated before taxonomic interpretation. A record therefore contributes one accession-level piece of evidence regardless of how many local alignment segments it contains.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "2"),
                  div(
                    h4("Best molecular match", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("Partial 100% matches are not allowed to outrank near-full-length evidence merely because their aligned segment is short. The app first prefers the >=90% coverage tier (or >=80% if needed), then keeps hits within 2 percentage points of the best query coverage in that tier. Identity is ranked first inside that near-best-coverage band; coverage, Bit score and E-value are tie-breakers.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "3"),
                  div(
                    h4("Close alternatives", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("A different named taxon is considered a close alternative when its best comparable accession is within 0.5 identity percentage points of the best match and is no more than 2 query-coverage points below it. These narrow windows are review heuristics, not universal species-delimitation thresholds.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "4"),
                  div(
                    h4("Conservative taxonomic fallback", span(class="evidence-badge evidence-heuristic", "App rule")),
                    p("If no close species alternative remains and sequence/reference evidence is strong, a species recommendation can be made. If several close species remain but they all belong to the same genus, species is reported as unresolved while the genus can remain a high-confidence recommendation. If close alternatives extend across genera, genus-level confidence is withheld and an LCA may be used as a fallback.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "5"),
                  div(
                    h4("Species evidence profile", span(class="evidence-badge evidence-heuristic", "Audit view")),
                    p("For every resolved species the app reports its strongest comparable accession, Identity, query coverage, Bit score, reference quality and the number of unique accessions returned for that species. The accession count describes database representation only; 23 records do not automatically defeat a better molecular match represented by one or two records.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "6"),
                  div(
                    h4("Sequence evidence and locus discrimination", span(class="evidence-badge evidence-published", "Biological context")),
                    p("Sequence quality and taxonomic discrimination are separate questions. A read can have excellent Identity and coverage yet fail to distinguish several species. In that situation the sequence evidence can remain High while the locus is flagged as having poor species-level discrimination. For ITS, broad literature benchmarks remain context for sequence-evidence strength rather than hard species cutoffs.")
                  )
                ),
                div(class="method-step",
                  div(class="method-num", "7"),
                  div(
                    h4("Reference context", span(class="evidence-badge evidence-published", "Reference quality")),
                    p("Type-material wording, RefSeq status, ordinary GenBank records and unresolved/environmental annotations are retained separately. Reference quality can change confidence, but it does not replace the direct Identity/coverage comparison.")
                  )
                )
              ),
              div(class="about-callout",
                strong("Interpretation principle: "),
                "Best molecular match first; close alternatives determine how far the identification can safely be taken; database frequency is supporting context only. A poor species-level result can therefore be a property of the locus rather than a failure of the sequence."
              )
            ),

            tabPanel("Scientific references",
              p(class="about-lead",
                "The links below are the main scientific and NCBI sources used to justify the biological context of the workflow. External links open in a new browser tab."),

              div(class="reference-card",
                h4("Schoch et al. 2012 — ITS as the primary fungal DNA barcode"),
                p("Compared candidate fungal barcode regions and proposed the nuclear ribosomal ITS region as the universal DNA barcode marker for Fungi."),
                tags$a(href="https://www.pnas.org/doi/10.1073/pnas.1117018109", target="_blank", rel="noopener noreferrer", "Open PNAS article ↗")
              ),

              div(class="reference-card",
                h4("Vu et al. 2018 — large-scale fungal barcode thresholds"),
                p("Large-scale analysis of filamentous fungal DNA barcodes. The broad ITS benchmarks used as context in the app (approximately 99.6% species and 94.3% genus) come from this work."),
                tags$a(href="https://pmc.ncbi.nlm.nih.gov/articles/PMC6020082/", target="_blank", rel="noopener noreferrer", "Open full text at NCBI/PMC ↗")
              ),

              div(class="reference-card",
                h4("Garnica et al. 2016 — lineage-dependent ITS thresholds"),
                p("Shows that optimal ITS similarity thresholds vary among lineages and that no single identity cutoff is universally reliable even within a diverse fungal genus."),
                tags$a(href="https://academic.oup.com/femsec/article/92/4/fiw045/2197947", target="_blank", rel="noopener noreferrer", "Open FEMS Microbiology Ecology article ↗")
              ),

              div(class="reference-card",
                h4("Stielow et al. 2015 — secondary fungal DNA barcodes"),
                p("Evaluated secondary fungal barcode loci and universal primers, supporting the use of additional loci when ITS alone does not provide sufficient resolution."),
                tags$a(href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4713107/", target="_blank", rel="noopener noreferrer", "Open full text at NCBI/PMC ↗")
              ),

              div(class="reference-card",
                h4("NCBI RefSeq Targeted Loci — fungal ITS"),
                p("Describes NCBI's curated fungal ITS reference collection. Sequences are mostly derived from type material and are maintained with specimen/taxonomic context."),
                tags$a(href="https://www.ncbi.nlm.nih.gov/refseq/targetedloci/", target="_blank", rel="noopener noreferrer", "Open NCBI RefSeq Targeted Loci ↗")
              ),

              div(class="reference-card",
                h4("NCBI BLAST documentation"),
                p("Technical definitions for BLAST output fields including Bit score, HSPs and query-coverage measures used by the application parser."),
                tags$a(href="https://www.ncbi.nlm.nih.gov/books/NBK279684/", target="_blank", rel="noopener noreferrer", "Open NCBI BLAST+ manual ↗")
              ),

              div(class="reference-card",
                h4("UNITE — fungal ITS reference and Species Hypotheses"),
                p("A fungal ITS-centered identification resource and Species Hypotheses system. It is included here as scientific context and a possible future curated-reference extension; the current app does not query UNITE directly."),
                tags$a(href="https://unite.ut.ee/", target="_blank", rel="noopener noreferrer", "Open UNITE ↗")
              ),

              div(class="about-callout",
                strong("Important distinction: "),
                "the near-best coverage band and the close-match windows (0.5 Identity percentage points and 2 coverage points) are application-specific review heuristics, not published species-delimitation thresholds. Published sources guide the biological context, while these rules should be calibrated against known isolates."
              )
            ),

            tabPanel("Version history",
              div(class="help-card",
                h3("CHANGELOG"),
                p("Every version records changes to processing, BLAST parsing, taxonomic interpretation and reporting so exported results can be traced back to the software logic used at the time."),
                verbatimTextOutput("changelog_text")
              )
            )
          )
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  rv <- reactiveValues(
    results = list(), summary = NULL, rename = NULL, settings = NULL,
    read_assignments = stage2_empty_assignments(), architecture = NULL,
    project_migration_log = "",
    blast_jobs = data.frame(
      final_name=character(), original_name=character(), rid=character(), rtoe=character(),
      database=character(), hitlist_size=integer(), status=character(), submitted_at=character(), last_checked_at=character(), stringsAsFactors=FALSE
    ),
    ncbi_last_contact = as.POSIXct(NA),
    blast_batch_status_text = "No batch operation has been run yet.",
    blast_raw = list(), blast_ids = data.frame(), blast_hits = data.frame(),
    taxonomy_summary = data.frame(), taxonomy_hits = data.frame(), taxonomy_counts = data.frame(),
    taxonomy_status_text = "No taxonomic analysis has been run yet.",
    taxonomy_batch_status_text = "No batch taxonomic analysis has been run yet.",
    project_status_text = "Current session has not been saved as a project.",
    project_loaded_name = "",
    context_peak_flag = NULL,
    pending_curation = NULL,
    auto_correct_preview_df = data.frame()
  )

  ensure_blast_jobs_schema <- function(df) {
    template <- data.frame(
      final_name=character(), original_name=character(), rid=character(), rtoe=character(),
      database=character(), hitlist_size=integer(), status=character(), submitted_at=character(),
      last_checked_at=character(), stringsAsFactors=FALSE
    )
    if (!is.data.frame(df) || !nrow(df)) return(template)
    if (!"database" %in% names(df)) df$database <- ""
    if (!"hitlist_size" %in% names(df)) df$hitlist_size <- NA_integer_
    for (nm in setdiff(names(template), names(df))) df[[nm]] <- template[[nm]][NA_integer_]
    df <- df[, names(template), drop=FALSE]
    df$database <- as.character(df$database)
    df$hitlist_size <- suppressWarnings(as.integer(df$hitlist_size))
    df
  }

  settings_for_result <- function(result, fallback = rv$settings) {
    if (is.list(result) && is.list(result$processing_settings)) result$processing_settings else fallback
  }

  current_upload_source_ids <- function() {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return(NULL)
    vapply(input$ab1_files$name, stage2_read_stem, character(1))
  }

  initialize_current_read_assignments <- function() {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return(stage2_empty_assignments())
    stage2_make_read_assignments(
      input$ab1_files$name,
      default_locus = input$target,
      default_direction = input$sequencing_primer,
      forward_primer = input$forward_primer,
      reverse_primer = input$reverse_primer
    )
  }

  sync_summary_from_results <- function() {
    if (is.null(rv$summary) || !is.data.frame(rv$summary) || !nrow(rv$summary) || is.null(rv$results) || !length(rv$results)) return(invisible(NULL))
    for (nm in names(rv$results)) {
      r <- ensure_curation_state(rv$results[[nm]])
      r <- curation_rebuild(r, settings_for_result(r))
      rv$results[[nm]] <- r
      sm <- r$summary
      idx <- which(as.character(rv$summary$sample_id) == nm)
      if (!length(idx) || !nrow(sm)) next
      for (col in names(sm)) {
        if (!col %in% names(rv$summary)) rv$summary[[col]] <- NA
        rv$summary[idx[1], col] <- sm[1, col]
      }
    }
    invisible(NULL)
  }

  # ---------------- Project Save / Load ----------------
  output$project_status <- renderUI({
    span(rv$project_status_text)
  })

  make_project_bundle <- function() {
    list(
      format = "SangerSequencePipelineProject",
      schema_version = PROJECT_SCHEMA_VERSION,
      app_version = APP_VERSION,
      saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      active_tab = if (!is.null(input$pipeline_step)) input$pipeline_step else "upload",
      ui_state = list(
        inspect_sample = input$inspect_sample,
        blast_sample = input$blast_sample,
        tax_sample = input$tax_sample,
        blast_database = input$blast_database,
        blast_hitlist = input$blast_hitlist
      ),
      state = list(
        results = rv$results,
        summary = rv$summary,
        rename = rv$rename,
        settings = rv$settings,
        read_assignments = rv$read_assignments,
        architecture = rv$architecture,
        migration_log = rv$project_migration_log,
        blast_jobs = rv$blast_jobs,
        blast_raw = rv$blast_raw,
        blast_ids = rv$blast_ids,
        blast_hits = normalize_blast_hits_unique_accession(rv$blast_hits),
        blast_batch_status_text = rv$blast_batch_status_text,
        taxonomy_summary = rv$taxonomy_summary,
        taxonomy_hits = rv$taxonomy_hits,
        taxonomy_counts = rv$taxonomy_counts,
        taxonomy_status_text = rv$taxonomy_status_text,
        taxonomy_batch_status_text = rv$taxonomy_batch_status_text
      )
    )
  }

  output$save_project <- downloadHandler(
    filename = function() {
      paste0(project_export_stem(), "_", format(Sys.Date(), "%Y%m%d"), ".sangerproject")
    },
    content = function(file) {
      saveRDS(make_project_bundle(), file = file, compress = "xz")
      rv$project_status_text <- paste0("Project saved from v", APP_VERSION, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ".")
    }
  )

  rebuild_blast_ids <- function() {
    if (!is.data.frame(rv$blast_hits) || !nrow(rv$blast_hits)) {
      rv$blast_ids <- data.frame()
      return(invisible(NULL))
    }
    hits <- normalize_blast_hits_unique_accession(rv$blast_hits)
    rv$blast_hits <- hits
    if (!all(c("original_name", "rid") %in% names(hits))) {
      rv$blast_ids <- data.frame()
      return(invisible(NULL))
    }
    parts <- split(seq_len(nrow(hits)), paste(hits$original_name, hits$rid, sep = "\r"))
    tops <- lapply(parts, function(idx) {
      g <- hits[idx, , drop = FALSE]
      if ("rank" %in% names(g)) g <- g[order(suppressWarnings(as.numeric(g$rank)), na.last = TRUE), , drop = FALSE]
      g[1, , drop = FALSE]
    })
    top <- do.call(rbind, tops)
    keep <- intersect(c(
      "final_name","original_name","rid","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","match_support"
    ), names(top))
    rv$blast_ids <- top[, keep, drop = FALSE]
    rownames(rv$blast_ids) <- NULL
  }

  invalidate_downstream_for_sample <- function(original_name, reason = "Curated sequence changed") {
    if (is.null(original_name) || !nzchar(original_name)) return(invisible(NULL))
    stale_rids <- character()
    if (is.data.frame(rv$blast_jobs) && nrow(rv$blast_jobs) && "original_name" %in% names(rv$blast_jobs)) {
      idx <- which(rv$blast_jobs$original_name == original_name)
      if (length(idx)) {
        stale_rids <- unique(as.character(rv$blast_jobs$rid[idx]))
        rv$blast_jobs$status[idx] <- "STALE"
      }
    }
    if (is.data.frame(rv$blast_hits) && nrow(rv$blast_hits)) {
      if ("original_name" %in% names(rv$blast_hits)) {
        stale_rids <- unique(c(stale_rids, as.character(rv$blast_hits$rid[rv$blast_hits$original_name == original_name])))
        rv$blast_hits <- rv$blast_hits[rv$blast_hits$original_name != original_name, , drop = FALSE]
      } else if (length(stale_rids) && "rid" %in% names(rv$blast_hits)) {
        rv$blast_hits <- rv$blast_hits[!rv$blast_hits$rid %in% stale_rids, , drop = FALSE]
      }
    }
    rebuild_blast_ids()
    filter_stale_tax <- function(df) {
      if (!is.data.frame(df) || !nrow(df)) return(df)
      if ("original_name" %in% names(df)) return(df[df$original_name != original_name, , drop = FALSE])
      if (length(stale_rids) && "rid" %in% names(df)) return(df[!df$rid %in% stale_rids, , drop = FALSE])
      df
    }
    rv$taxonomy_summary <- filter_stale_tax(rv$taxonomy_summary)
    rv$taxonomy_hits <- filter_stale_tax(rv$taxonomy_hits)
    rv$taxonomy_counts <- filter_stale_tax(rv$taxonomy_counts)
    rv$blast_batch_status_text <- paste0("Sequence ", original_name, " changed after curation; previous BLAST results for this sample are stale and must be rerun.")
    rv$taxonomy_status_text <- paste0("Sequence ", original_name, " changed after curation; taxonomic interpretation was invalidated.")
    invisible(NULL)
  }

  commit_curated_result <- function(sample_name, new_result, label = "Manual curation") {
    old_result <- rv$results[[sample_name]]
    old_seq <- if (!is.null(old_result$seq)) as.character(old_result$seq) else ""
    new_seq <- if (!is.null(new_result$seq)) as.character(new_result$seq) else ""
    rv$results[[sample_name]] <- new_result
    sync_summary_from_results()
    if (identical(input$inspect_sample, sample_name)) {
      updateTextAreaInput(session, "trimmed_sequence_preview", value = new_seq)
    }
    if (!identical(old_seq, new_seq)) {
      invalidate_downstream_for_sample(sample_name, label)
    }
    rv$project_status_text <- paste0("Unsaved curation change: ", sample_name, " · ", label, ".")
    invisible(!identical(old_seq, new_seq))
  }

  restore_settings_inputs <- function(settings) {
    if (is.null(settings) || !is.list(settings)) return(invisible(NULL))
    if (!is.null(settings$target)) updateSelectInput(session, "target", selected = settings$target)
    if (!is.null(settings$forward_primer)) updateTextInput(session, "forward_primer", value = settings$forward_primer)
    if (!is.null(settings$forward_primer_seq)) updateTextInput(session, "forward_primer_seq", value = settings$forward_primer_seq)
    if (!is.null(settings$reverse_primer)) updateTextInput(session, "reverse_primer", value = settings$reverse_primer)
    if (!is.null(settings$reverse_primer_seq)) updateTextInput(session, "reverse_primer_seq", value = settings$reverse_primer_seq)
    if (!is.null(settings$sequencing_primer)) updateRadioButtons(session, "sequencing_primer", selected = settings$sequencing_primer)
    if (!is.null(settings$enable_primer_mapping)) updateCheckboxInput(session, "enable_primer_mapping", value = isTRUE(settings$enable_primer_mapping))
    if (!is.null(settings$expected_amplicon_len)) updateNumericInput(session, "expected_amplicon_len", value = settings$expected_amplicon_len)
    if (!is.null(settings$absolute_max_base_index)) updateNumericInput(session, "absolute_max_base_index", value = settings$absolute_max_base_index)
    if (!is.null(settings$window)) updateNumericInput(session, "window", value = settings$window)
    if (!is.null(settings$min_peak_ratio)) updateNumericInput(session, "min_peak_ratio", value = settings$min_peak_ratio)
    if (!is.null(settings$min_relative_signal)) updateNumericInput(session, "min_relative_signal", value = settings$min_relative_signal)
    if (!is.null(settings$min_len_before_collapse)) updateNumericInput(session, "min_len_before_collapse", value = settings$min_len_before_collapse)
    if (!is.null(settings$bad_run_windows)) updateNumericInput(session, "bad_run_windows", value = settings$bad_run_windows)
    if (!is.null(settings$min_usable_len)) updateNumericInput(session, "min_usable_len", value = settings$min_usable_len)
  }

  observeEvent(input$load_project, {
    req(input$load_project$datapath)
    obj <- tryCatch(readRDS(input$load_project$datapath), error = function(e) structure(list(error = conditionMessage(e)), class = "project_load_error"))
    if (inherits(obj, "project_load_error")) {
      showNotification(paste("Could not load project:", obj$error), type = "error", duration = 10)
      return()
    }
    if (!is.list(obj) || !identical(obj$format, "SangerSequencePipelineProject") || is.null(obj$state)) {
      showNotification("This file is not a valid Sanger Sequence Pipeline project.", type = "error", duration = 10)
      return()
    }
    source_schema <- if (is.null(obj$schema_version)) 1L else suppressWarnings(as.integer(obj$schema_version))
    if (length(source_schema) != 1L || is.na(source_schema) || source_schema < 1L) {
      showNotification("This project has an invalid schema version.", type = "error", duration = 10)
      return()
    }
    if (source_schema > PROJECT_SCHEMA_VERSION) {
      showNotification("This project was created by a newer project schema and cannot be loaded safely.", type = "error", duration = 10)
      return()
    }

    st <- obj$state
    if (source_schema < 2L) st <- stage2_migrate_v1_state(st)
    if (source_schema == 2L) st <- stage2_migrate_v2_state(st)
    loaded_assignments <- stage2_coerce_assignments(st$read_assignments)
    assignment_error <- stage2_validate_assignments(loaded_assignments)
    if (length(st$results) && !is.null(assignment_error)) {
      showNotification(paste("Project read architecture is invalid:", assignment_error), type = "error", duration = 10)
      return()
    }
    loaded_architecture <- if (nrow(loaded_assignments) && is.null(assignment_error)) tryCatch(
      stage2_build_architecture(
        loaded_assignments,
        project_id = if (is.list(st$architecture)) stage2_scalar_text(st$architecture$project_id, "project") else "project"
      ),
      error = function(e) NULL
    ) else NULL
    if (length(st$results) && nrow(loaded_assignments) && is.null(loaded_architecture)) {
      showNotification("Project read architecture could not be rebuilt safely.", type = "error", duration = 10)
      return()
    }
    rv$results <- if (!is.null(st$results)) st$results else list()
    if (length(rv$results)) {
      for (nm in names(rv$results)) {
        rv$results[[nm]] <- ensure_curation_state(rv$results[[nm]])
        loaded_settings <- if (is.list(rv$results[[nm]]$processing_settings)) rv$results[[nm]]$processing_settings else st$settings
        rv$results[[nm]] <- curation_rebuild(rv$results[[nm]], loaded_settings)
      }
    }
    rv$summary <- st$summary
    rv$rename <- st$rename
    rv$settings <- st$settings
    rv$read_assignments <- loaded_assignments
    rv$architecture <- loaded_architecture
    rv$project_migration_log <- stage2_scalar_text(st$migration_log)
    rv$blast_jobs <- ensure_blast_jobs_schema(if (is.data.frame(st$blast_jobs)) st$blast_jobs else NULL)
    rv$blast_raw <- if (is.list(st$blast_raw)) st$blast_raw else list()
    rv$blast_hits <- if (is.data.frame(st$blast_hits)) normalize_blast_hits_unique_accession(st$blast_hits) else data.frame()
    rv$blast_batch_status_text <- if (!is.null(st$blast_batch_status_text)) st$blast_batch_status_text else "Loaded project."
    rv$taxonomy_summary <- if (is.data.frame(st$taxonomy_summary)) st$taxonomy_summary else data.frame()
    rv$taxonomy_hits <- if (is.data.frame(st$taxonomy_hits)) st$taxonomy_hits else data.frame()
    rv$taxonomy_counts <- if (is.data.frame(st$taxonomy_counts)) st$taxonomy_counts else data.frame()
    rv$taxonomy_status_text <- if (!is.null(st$taxonomy_status_text)) st$taxonomy_status_text else "No taxonomic analysis has been run yet."
    rv$taxonomy_batch_status_text <- if (!is.null(st$taxonomy_batch_status_text)) st$taxonomy_batch_status_text else "No batch taxonomic analysis has been run yet."
    rv$ncbi_last_contact <- as.POSIXct(NA)
    sync_summary_from_results()
    rebuild_blast_ids()

    restore_settings_inputs(rv$settings)

    result_names <- names(rv$results)
    saved_ui <- obj$ui_state
    inspect_selected <- if (!is.null(saved_ui$inspect_sample) && saved_ui$inspect_sample %in% result_names) saved_ui$inspect_sample else if (length(result_names)) result_names[1] else character()
    sync_qc_sample_choices(inspect_selected)

    final_names <- character()
    if (length(result_names)) {
      result_labels <- result_names
      if (!is.null(rv$rename) && nrow(rv$rename)) {
        rename_idx <- match(result_names, rv$rename$Original_name)
        resolved <- !is.na(rename_idx) & nzchar(trimws(rv$rename$New_name[rename_idx]))
        result_labels[resolved] <- rv$rename$New_name[rename_idx[resolved]]
      }
      final_names <- setNames(result_names, result_labels)
    }
    blast_selected <- if (!is.null(saved_ui$blast_sample) && saved_ui$blast_sample %in% unname(final_names)) saved_ui$blast_sample else if (length(final_names)) unname(final_names)[1] else character()
    updateSelectInput(session, "blast_sample", choices = final_names, selected = blast_selected)
    if (!is.null(saved_ui$blast_database)) updateSelectInput(session, "blast_database", selected = saved_ui$blast_database)
    if (!is.null(saved_ui$blast_hitlist)) updateNumericInput(session, "blast_hitlist", value = saved_ui$blast_hitlist)

    if (nrow(rv$blast_hits)) {
      pairs <- unique(rv$blast_hits[, intersect(c("original_name","final_name"), names(rv$blast_hits)), drop = FALSE])
      if (all(c("original_name","final_name") %in% names(pairs))) {
        tax_choices <- setNames(as.character(pairs$original_name), as.character(pairs$final_name))
        tax_selected <- if (!is.null(saved_ui$tax_sample) && saved_ui$tax_sample %in% unname(tax_choices)) saved_ui$tax_sample else if (length(tax_choices)) unname(tax_choices)[1] else character()
        updateSelectInput(session, "tax_sample", choices = tax_choices, selected = tax_selected)
      }
    } else {
      updateSelectInput(session, "tax_sample", choices = character())
    }

    active <- if (!is.null(obj$active_tab) && obj$active_tab %in% c("upload","settings","qc","rename","export","blast","taxonomy","help")) obj$active_tab else if (nrow(rv$taxonomy_summary)) "taxonomy" else if (nrow(rv$blast_hits)) "blast" else if (length(rv$results)) "qc" else "upload"
    updateTabsetPanel(session, "pipeline_step", selected = active)

    rv$project_loaded_name <- input$load_project$name
    rv$project_status_text <- paste0(
      "Loaded ", input$load_project$name,
      " · saved with app v", ifelse(is.null(obj$app_version), "unknown", obj$app_version),
      if (!is.null(obj$saved_at)) paste0(" · saved ", obj$saved_at) else "",
      if (source_schema < PROJECT_SCHEMA_VERSION) paste0(" · migrated project schema ", source_schema, " → ", PROJECT_SCHEMA_VERSION, " in memory") else "",
      ". BLAST hits were normalized to one row per accession."
    )
    showNotification(
      if (source_schema < PROJECT_SCHEMA_VERSION) paste0("Older project loaded and migrated to schema ", PROJECT_SCHEMA_VERSION, ". Save it to persist the migration.") else "Project loaded successfully.",
      type = "message", duration = 8
    )
  })

  # ---------------- Upload ----------------
  output$uploaded_files_table <- renderDT({
    req(input$ab1_files)
    datatable(data.frame(File=input$ab1_files$name, Size_KB=round(input$ab1_files$size/1024,1)),
              rownames=FALSE, options=list(pageLength=15, dom="tip"))
  })

  observeEvent(input$ab1_files, {
    fresh <- initialize_current_read_assignments()
    previous <- stage2_coerce_assignments(rv$read_assignments)
    if (nrow(previous) && nrow(fresh)) {
      carry <- intersect(c("Isolate", "Locus", "Direction", "Primer", "Notes"), names(previous))
      for (i in seq_len(nrow(fresh))) {
        j <- match(fresh$Source_ID[i], previous$Source_ID)
        if (!is.na(j)) fresh[i, carry] <- previous[j, carry]
      }
    }
    rv$read_assignments <- fresh
    sync_assignment_state()
    rv$project_migration_log <- ""
  })

  sync_assignment_state <- function() {
    if (!nrow(rv$read_assignments)) {
      rv$rename <- stage2_default_rename_map(rv$read_assignments)
      rv$architecture <- NULL
      return(invisible(NULL))
    }
    rv$read_assignments <- stage2_sync_generated_names(
      rv$read_assignments,
      forward_primer = input$forward_primer,
      reverse_primer = input$reverse_primer
    )
    rv$rename <- stage2_default_rename_map(rv$read_assignments)
    assignment_error <- stage2_identity_error(rv$read_assignments)
    rv$architecture <- if (is.null(assignment_error)) stage2_build_architecture(rv$read_assignments) else NULL
    if (is.null(assignment_error) && length(rv$results)) {
      for (source_id in intersect(names(rv$results), rv$read_assignments$Source_ID)) {
        i <- match(source_id, rv$read_assignments$Source_ID)
        rv$results[[source_id]]$read_assignment <- as.list(rv$read_assignments[i, , drop = FALSE])
        if (is.list(rv$results[[source_id]]$processing_settings)) {
          rv$results[[source_id]]$processing_settings$target <- rv$read_assignments$Locus[i]
          rv$results[[source_id]]$processing_settings$sequencing_primer <- rv$read_assignments$Direction[i]
          if (rv$read_assignments$Direction[i] == "Forward" && nzchar(rv$read_assignments$Primer[i])) rv$results[[source_id]]$processing_settings$forward_primer <- rv$read_assignments$Primer[i]
          if (rv$read_assignments$Direction[i] == "Reverse" && nzchar(rv$read_assignments$Primer[i])) rv$results[[source_id]]$processing_settings$reverse_primer <- rv$read_assignments$Primer[i]
        }
      }
    }
    invisible(assignment_error)
  }

  observeEvent(list(input$target, input$sequencing_primer, input$forward_primer, input$reverse_primer), {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) return()
    sync_assignment_state()
  }, ignoreInit = TRUE)

  architecture_summary_ui <- function() {
    error <- stage2_identity_error(rv$read_assignments)
    if (!is.null(error)) return(div(class = "status-warning", "Architecture preview will appear after every read has explicit Isolate, Gene / locus and Forward / Reverse fields."))
    architecture <- tryCatch(stage2_build_architecture(rv$read_assignments), error = function(e) NULL)
    if (is.null(architecture)) return(NULL)
    sm <- stage2_architecture_summary(architecture)
    div(
      class = "compact-hint",
      paste0(
        "Architecture preview: ", sm$Isolates, " isolate(s) · ", sm$Loci, " locus/loci · ", sm$Reads, " read(s) · ",
        sm$Paired_loci, " Forward/Reverse pair(s) · ", sm$Single_read_loci, " single-read locus/loci."
      )
    )
  }
  output$architecture_summary <- renderUI(architecture_summary_ui())
  output$upload_architecture_summary <- renderUI(architecture_summary_ui())

  observeEvent(input$to_settings, {
    if (is.null(input$ab1_files) || nrow(input$ab1_files)==0) {
      showNotification("Please upload at least one AB1 file.", type="error"); return()
    }
    updateTabsetPanel(session, "pipeline_step", selected="settings")
  })
  observeEvent(input$back_upload, updateTabsetPanel(session,"pipeline_step",selected="upload"))
  observeEvent(input$back_settings_from_rename, updateTabsetPanel(session,"pipeline_step",selected="settings"))
  observeEvent(input$back_rename_from_qc, updateTabsetPanel(session,"pipeline_step",selected="rename"))
  observeEvent(input$back_qc_from_export, updateTabsetPanel(session,"pipeline_step",selected="qc"))
  observeEvent(input$back_export, updateTabsetPanel(session,"pipeline_step",selected="export"))
  observeEvent(input$back_blast, updateTabsetPanel(session,"pipeline_step",selected="blast"))

  current_settings_from_inputs <- function() {
    list(
      target=input$target,
      forward_primer=input$forward_primer,
      forward_primer_seq=sanitize_dna(input$forward_primer_seq),
      reverse_primer=input$reverse_primer,
      reverse_primer_seq=sanitize_dna(input$reverse_primer_seq),
      sequencing_primer=input$sequencing_primer,
      enable_primer_mapping=isTRUE(input$enable_primer_mapping),
      expected_amplicon_len=as.integer(input$expected_amplicon_len),
      absolute_max_base_index=as.integer(input$absolute_max_base_index),
      window=as.integer(input$window),
      min_peak_ratio=as.numeric(input$min_peak_ratio),
      min_relative_signal=as.numeric(input$min_relative_signal),
      min_len_before_collapse=as.integer(input$min_len_before_collapse),
      bad_run_windows=as.integer(input$bad_run_windows),
      min_usable_len=as.integer(input$min_usable_len),
      ambiguous_peak_strong_ratio=1.25,
      ambiguous_peak_moderate_ratio=1.75,
      ambiguous_peak_min_relative_signal=0.20,
      auto_correct_min_alt_to_called=1.80,
      auto_correct_min_alt_to_third=2.00,
      auto_correct_max_peak_offset=2L,
      auto_correct_min_relative_signal=0.50
    )
  }

  observeEvent(input$to_rename, {
    if (is.null(input$ab1_files) || !nrow(input$ab1_files)) {
      showNotification("Please upload at least one AB1 file.", type = "error")
      return()
    }
    source_ids <- current_upload_source_ids()
    if (!nrow(rv$read_assignments) || !identical(sort(rv$read_assignments$Source_ID), sort(source_ids))) {
      rv$read_assignments <- initialize_current_read_assignments()
    }
    rv$settings <- current_settings_from_inputs()
    sync_assignment_state()
    updateTabsetPanel(session, "pipeline_step", selected = "rename")
  })

  # ---------------- Trimming ----------------
  observeEvent(input$run_trimming, {
    req(input$ab1_files)
    name_error <- stage2_identity_error(rv$read_assignments)
    if (!is.null(name_error)) {
      showNotification(paste("Rename error:", name_error), type = "error", duration = 10)
      return()
    }
    sync_assignment_state()
    assignment_error <- stage2_validate_assignments(rv$read_assignments, current_upload_source_ids())
    if (!is.null(assignment_error)) {
      showNotification(paste("Read assignment error:", assignment_error), type = "error", duration = 10)
      return()
    }
    rv$read_assignments <- stage2_coerce_assignments(rv$read_assignments)
    rv$architecture <- tryCatch(
      stage2_build_architecture(rv$read_assignments),
      error = function(e) {
        showNotification(paste("Could not build project architecture:", conditionMessage(e)), type = "error", duration = 10)
        NULL
      }
    )
    if (is.null(rv$architecture)) return()
    settings <- current_settings_from_inputs()
    rv$settings <- settings
    all_results <- list(); summaries <- list(); files <- input$ab1_files

    withProgress(message="Processing AB1 files", value=0, {
      for (i in seq_len(nrow(files))) {
        sample_id <- sub("\\.ab1$", "", files$name[i], ignore.case=TRUE)
        assignment_idx <- match(sample_id, rv$read_assignments$Source_ID)
        assignment <- rv$read_assignments[assignment_idx, , drop = FALSE]
        read_settings <- settings
        read_settings$target <- assignment$Locus[1]
        read_settings$sequencing_primer <- assignment$Direction[1]
        if (assignment$Direction[1] == "Forward" && nzchar(trimws(assignment$Primer[1]))) read_settings$forward_primer <- assignment$Primer[1]
        if (assignment$Direction[1] == "Reverse" && nzchar(trimws(assignment$Primer[1]))) read_settings$reverse_primer <- assignment$Primer[1]
        incProgress(1/nrow(files), detail=paste("Processing", files$name[i]))
        result <- tryCatch(
          trim_one_ab1(files$datapath[i], sample_id, read_settings),
          error=function(e) structure(list(error=conditionMessage(e)), class="ab1_error")
        )
        if (inherits(result,"ab1_error")) {
          summaries[[sample_id]] <- make_failure_summary(sample_id, read_settings, result$error)
        } else {
          result$read_assignment <- as.list(assignment[1, , drop = FALSE])
          result$processing_settings <- read_settings
          result <- ensure_curation_state(result)
          result <- curation_rebuild(result, read_settings)
          all_results[[sample_id]] <- result
          summaries[[sample_id]] <- result$summary
        }
      }
    })

    rv$results <- all_results
    rv$summary <- Reduce(rbind_fill, summaries); rownames(rv$summary) <- NULL
    sync_qc_sample_choices()
    session$sendCustomMessage("showLoader", list(text = "Opening Trim & QC workspace…"))
    updateTabsetPanel(session,"pipeline_step",selected="qc")
  })

  qc_display_name <- function(original_name) {
    original_name <- as.character(original_name)[1]
    if (is.null(rv$rename) || !nrow(rv$rename)) return(original_name)
    idx <- match(original_name, rv$rename$Original_name)
    if (is.na(idx) || !nzchar(trimws(rv$rename$New_name[idx]))) original_name else as.character(rv$rename$New_name[idx])
  }

  sync_qc_sample_choices <- function(preferred = NULL) {
    keys <- names(rv$results)
    if (!length(keys)) {
      updateSelectInput(session, "inspect_sample", choices = character(), selected = character())
      return(invisible(NULL))
    }
    labels <- vapply(keys, qc_display_name, character(1))
    choices <- stats::setNames(keys, labels)
    current <- if (!is.null(preferred) && preferred %in% keys) preferred else isolate(input$inspect_sample)
    selected <- if (!is.null(current) && length(current) == 1L && current %in% keys) current else if (length(keys)) keys[1] else character()
    updateSelectInput(session, "inspect_sample", choices = choices, selected = selected)
  }

  # ---------------- QC summary ----------------
  output$qc_summary_cards <- renderUI({
    req(rv$summary)
    vals <- c(
      Total=nrow(rv$summary),
      OK=sum(rv$summary$status=="OK",na.rm=TRUE),
      Warnings=sum(rv$summary$status=="SHORT_AFTER_TRIMMING",na.rm=TRUE),
      Failed=sum(rv$summary$status %in% c("FAILED_TRIMMING","ERROR"),na.rm=TRUE)
    )
    fluidRow(lapply(names(vals), function(nm) column(3, div(class="summary-card", div(class="summary-number",vals[[nm]]), div(class="summary-label",nm)))))
  })

  output$summary_table <- renderDT({
    req(rv$summary)
    df <- rv$summary[,c("sample_id","target","raw_length","trimmed_length","trim_start","trim_end","collapse_index","reason","median_peak_ratio_trimmed","status")]
    df$sample_id <- vapply(df$sample_id, qc_display_name, character(1))
    names(df) <- c("Sample","Target","Raw length","Trimmed length","Start","End","Collapse","Reason","Median peak ratio","Status")
    datatable(df, rownames=FALSE, filter="top", options=list(pageLength=15,scrollX=TRUE))
  })

  selected_sample_key <- reactive({
    sid <- input$inspect_sample
    req(!is.null(sid), length(sid) == 1L, nzchar(sid), sid %in% names(rv$results))
    sid
  })

  selected_result <- reactive({
    sid <- selected_sample_key()
    r <- rv$results[[sid]]
    req(!is.null(r))
    r
  })

  selected_processing_settings <- reactive({
    settings <- settings_for_result(selected_result())
    req(is.list(settings))
    settings
  })

  observeEvent(input$inspect_sample, {
    r <- selected_result()
    updateTextAreaInput(session, "trimmed_sequence_preview", value = r$seq)
  })

  output$sequence_metrics <- renderDT({
    datatable(make_sequence_preview(selected_result()), rownames = FALSE, options = list(dom = "t"))
  })

  # ---------------- PITAX 3.0 Stage 1: AB1 evidence audit ----------------
  output$ab1_evidence_run_table <- renderDT({
    req(rv$results)
    df <- ab1_evidence_run_summary(rv$results)
    if (!nrow(df)) {
      return(datatable(data.frame(Message = "No AB1 evidence has been captured in this run."),
                       rownames = FALSE, selection = "none", options = list(dom = "t")))
    }
    display <- df
    display$Sample <- vapply(display$Sample, qc_display_name, character(1))
    names(display) <- c(
      "Sample", "Evidence", "Quality tag", "Quality coverage (%)",
      "Auto trim start", "Auto trim end", "Auto trim length", "Median quality · auto trim",
      "Q≥20 · auto trim (%)", "Q≥30 · auto trim (%)",
      "PCON window start", "PCON window end", "PCON window length", "PCON window shift",
      "Median quality · PCON window", "Q≥20 · PCON window (%)", "Q≥30 · PCON window (%)",
      "Primary position source",
      "Primary position coverage (%)", "Primary positions different (%)", "Median |position Δ|",
      "Legacy map", "Canonical A/C/G/T map", "Maps match",
      "Legacy call is max (%)", "Canonical call is max (%)",
      "Canonical call is max · auto trim (%)", "Median called/alternative ratio · auto trim"
    )
    datatable(
      display,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      options = list(pageLength = 12, scrollX = TRUE, dom = "tip")
    )
  })

  # The run-audit table and the Sample dropdown now share one source of truth.
  # Clicking an audit row selects the same sample throughout the QC workspace.
  observeEvent(input$ab1_evidence_run_table_rows_selected, {
    idx <- input$ab1_evidence_run_table_rows_selected
    if (is.null(idx) || length(idx) != 1L || idx < 1L) return()
    df <- ab1_evidence_run_summary(rv$results)
    if (idx > nrow(df)) return()
    key <- pitax_result_key_for_sample(rv$results, df$Sample[idx])
    if (nzchar(key)) updateSelectInput(session, "inspect_sample", selected = key)
  }, ignoreInit = TRUE)

  output$ab1_evidence_selected_note <- renderUI({
    key <- selected_sample_key()
    r <- selected_result()
    sid <- pitax_result_sample_id(r, key)
    display_name <- qc_display_name(key)
    sample_label <- if (!identical(display_name, sid)) paste0(display_name, " (original: ", sid, ")") else sid
    ev <- r$ab1_evidence
    if (is.null(ev)) {
      return(div(class = "status-note",
                 tags$strong(paste0("Selected sample: ", sample_label, ". ")),
                 "This sample has no Stage 1 evidence object. Reprocess the original AB1 with alpha.6 to create the audit."))
    }
    if (!is.null(ev$error) && nzchar(as.character(ev$error)[1])) {
      return(div(class = "status-error", tags$strong(paste0("Selected sample: ", sample_label, ". ")),
                 paste0("Evidence audit error: ", as.character(ev$error)[1],
                        ". The established trimming result was preserved.")))
    }
    sm <- ab1_evidence_result_summary(r)
    tagList(
      div(class = "compact-hint", tags$strong(paste0("Selected sample: ", sample_label))),
      div(class = "peak-flag-summary",
          div(class = "peak-flag-pill", tags$strong(sm$Quality_tag[1]), " quality tag"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Median_quality_auto_trim[1]), sm$Median_quality_auto_trim[1], "NA")), " median quality in auto trim"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Q20_auto_trim_percent[1]), paste0(sm$Q20_auto_trim_percent[1], "%"), "NA")), " Q≥20 in auto trim"),
          div(class = "peak-flag-pill", tags$strong(ifelse(is.finite(sm$Primary_position_difference_percent[1]), paste0(sm$Primary_position_difference_percent[1], "%"), "NA")), " PLOC vs legacy positions differ")),
      div(class = "compact-hint",
          "alpha.6 retains the active legacy auto trim and a same-length PCON-only comparison. The comparison is observational and does not change the processed sequence.")
    )
  })

  output$ab1_trim_comparison_table <- renderDT({
    d <- pitax_trim_window_comparison(selected_result())
    names(d) <- c(
      "Window", "Status", "Start", "End", "Length", "Quality coverage (%)",
      "Median quality", "Q≥20 (%)", "Q≥30 (%)", "Start shift"
    )
    datatable(d, rownames = FALSE, selection = "none", options = list(dom = "t", scrollX = TRUE))
  })

  output$ab1_evidence_detail_table <- renderDT({
    key <- selected_sample_key()
    r <- selected_result()
    ev <- r$ab1_evidence
    if (is.null(ev) || !is.data.frame(ev$detail) || !nrow(ev$detail)) {
      return(datatable(data.frame(Message = "No per-base Stage 1 evidence available for this sample."),
                       rownames = FALSE, selection = "none", options = list(dom = "t")))
    }
    d <- pitax_add_trim_membership(ev$detail, r)
    d$Sample_ID <- pitax_result_sample_id(r, key)
    d$Final_Name <- qc_display_name(key)
    keep <- c(
      "Sample_ID", "Final_Name", "Position", "Base", "In_auto_trim", "In_quality_proposed_window", "Basecaller_quality",
      "Legacy_primary_peak_pos", "Raw_ABIF_primary_peak_pos", "Primary_peak_pos_delta",
      "Legacy_called_signal", "Legacy_called_is_max",
      "Canonical_called_signal", "Canonical_best_alt_signal", "Canonical_called_to_alt_ratio",
      "Canonical_best_channel", "Canonical_called_is_max"
    )
    d <- d[, intersect(keep, names(d)), drop = FALSE]
    numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    for (nm in numeric_cols) d[[nm]] <- round(d[[nm]], 3)
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE, dom = "tip"))
  })

  output$download_ab1_evidence_run <- downloadHandler(
    filename = function() paste0("PITAX_v3_stage1_AB1_run_audit_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      df <- ab1_evidence_run_summary(rv$results)
      df$Final_Name <- vapply(df$Sample, qc_display_name, character(1))
      df <- df[, c("Sample", "Final_Name", setdiff(names(df), c("Sample", "Final_Name"))), drop = FALSE]
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )

  output$download_ab1_evidence_detail <- downloadHandler(
    filename = function() {
      key <- selected_sample_key()
      r <- selected_result()
      pitax_assert_export_identity(r, key)
      paste0(clean_fasta_name(qc_display_name(key)), "_PITAX_v3_stage1_AB1_base_audit.csv")
    },
    content = function(file) {
      key <- selected_sample_key()
      r <- selected_result()
      d <- pitax_evidence_detail_export(r, key)
      d$Final_Name <- qc_display_name(key)
      d <- d[, c("Sample_ID", "Final_Name", setdiff(names(d), c("Sample_ID", "Final_Name"))), drop = FALSE]
      utils::write.csv(d, file, row.names = FALSE, na = "")
    }
  )
  output$primer_match_table <- renderDT({
    req(isTRUE(input$enable_primer_mapping))
    datatable(primer_match_table(selected_result(), selected_processing_settings()), rownames = FALSE,
              options = list(dom = "t", scrollX = TRUE))
  })
  output$primer_alignment <- renderText({
    req(isTRUE(input$enable_primer_mapping))
    primer_alignment_text(selected_result(), selected_processing_settings())
  })
  output$amplicon_overview <- renderPlot({
    draw_amplicon_overview(selected_result(), selected_processing_settings())
  })
  output$expected_amplicon_note <- renderUI({
    settings <- selected_processing_settings()
    span(class = "peak-flag-pill", title = "Assay metadata only; not a fixed coordinate on the Sanger read.",
         paste0("Expected amplicon: ", settings$expected_amplicon_len, " bp · Locus: ", settings$target))
  })

  all_current_peak_flags <- reactive({
    req(selected_result())
    scope <- if (!is.null(input$peak_flag_scope) && input$peak_flag_scope %in% c("trimmed", "raw")) input$peak_flag_scope else "trimmed"
    df <- ambiguous_peak_flags(selected_result(), scope = scope, params = ambiguous_peak_params_from_settings(selected_processing_settings()))
    if (!nrow(df)) return(df)
    r <- ensure_curation_state(selected_result())
    reviewed <- as.integer(r$curation$reviewed_positions)
    df$Review_status <- ifelse(df$Position %in% reviewed, "Reviewed", "Active")
    df
  })

  current_peak_flags <- reactive({
    df <- all_current_peak_flags()
    if (!nrow(df)) return(df)
    if (!isTRUE(input$show_reviewed_flags)) df <- df[df$Review_status != "Reviewed", , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  center_chromatogram_at <- function(pos, half_window = 20) {
    n <- length(selected_result()$peak_pos)
    pos <- suppressWarnings(as.numeric(pos))
    if (!is.finite(pos) || n < 1) return(invisible(NULL))
    range <- c(max(1, pos - half_window), min(n, pos + half_window))
    proxy <- plotly::plotlyProxy("chromatogram_plot", session)
    plotly::plotlyProxyInvoke(proxy, "relayout", list(xaxis.range = range))
    invisible(NULL)
  }

  output$ambiguous_peak_summary <- renderUI({
    all_df <- all_current_peak_flags()
    active <- if (nrow(all_df)) all_df[all_df$Review_status == "Active", , drop = FALSE] else all_df
    total <- nrow(active)
    strong_count <- if (total) sum(active$Severity == "Strong", na.rm = TRUE) else 0L
    moderate_count <- if (total) sum(active$Severity == "Moderate", na.rm = TRUE) else 0L
    auto_count <- if (total && "Auto_correct_candidate" %in% names(active)) sum(active$Auto_correct_candidate %in% TRUE, na.rm = TRUE) else 0L
    reviewed_count <- if (nrow(all_df)) sum(all_df$Review_status == "Reviewed", na.rm = TRUE) else 0L
    r <- ensure_curation_state(selected_result())
    div(class = "peak-flag-summary",
        div(class = "peak-flag-pill", tags$strong(total), " active"),
        div(class = "peak-flag-pill", tags$strong(strong_count), " strong"),
        div(class = "peak-flag-pill", tags$strong(moderate_count), " moderate"),
        if (auto_count > 0) div(class = "peak-flag-pill", tags$strong(auto_count), " high-confidence proposals") else NULL,
        if (reviewed_count > 0) div(class = "peak-flag-pill", tags$strong(reviewed_count), " reviewed") else NULL,
        if (r$curation$revision > 0) div(class = "manual-edit-badge", paste0("Curation revision ", r$curation$revision)) else NULL)
  })

  output$ambiguous_peak_table <- renderDT({
    df <- current_peak_flags()
    if (!nrow(df)) {
      return(datatable(
        data.frame(Message = "No active ambiguous channel-competition positions in the selected scope."),
        rownames = FALSE, selection = "none", options = list(dom = "t")
      ))
    }
    show <- df[, c("Position", "Call", "Competing_channel", "Peak_ratio", "Competitor_percent",
                   "Competitor_peak_offset", "Severity", "Flag", "Auto_correct_candidate", "Review_status"), drop = FALSE]
    show$Auto_correct_candidate <- ifelse(show$Auto_correct_candidate %in% TRUE, "Yes", "")
    names(show) <- c("Position", "Call", "Competing channel", "Peak ratio", "Competitor / called (%)",
                     "Peak offset", "Severity", "Flag", "Auto", "Status")
    datatable(
      show,
      rownames = FALSE,
      selection = list(mode = "single", selected = NULL, target = "row"),
      filter = "none",
      callback = JS("
        table.on('contextmenu', 'tbody tr', function(e) {
          e.preventDefault();
          var row = table.row(this).data();
          if (!row || row.length < 1) return;
          Shiny.setInputValue('peak_flag_context', {position: row[0], nonce: Math.random()}, {priority: 'event'});
        });
      "),
      options = list(pageLength = 10, lengthChange = FALSE, scrollX = TRUE, dom = "tip", order = list(list(0, "asc")))
    )
  })

  observeEvent(input$ambiguous_peak_table_rows_selected, {
    idx <- input$ambiguous_peak_table_rows_selected
    df <- current_peak_flags()
    if (!length(idx) || !nrow(df) || idx[1] < 1 || idx[1] > nrow(df)) return()
    center_chromatogram_at(df$Position[idx[1]])
  }, ignoreInit = TRUE)

  flag_evidence_text <- function(flag_row) {
    if (is.null(flag_row) || !nrow(flag_row)) return("")
    paste0(
      "Call ", flag_row$Call[1], "; competitor ", flag_row$Competing_channel[1],
      "; called signal ", flag_row$Called_signal[1],
      "; competitor signal ", flag_row$Competitor_signal[1],
      "; peak ratio ", flag_row$Peak_ratio[1],
      "; competitor/called ", flag_row$Competitor_percent[1], "%",
      if ("Competitor_peak_offset" %in% names(flag_row) && is.finite(flag_row$Competitor_peak_offset[1])) paste0("; peak offset ", flag_row$Competitor_peak_offset[1]) else "",
      "; flag: ", flag_row$Flag[1]
    )
  }

  show_peak_curation_menu <- function(flag_row) {
    if (is.null(flag_row) || !nrow(flag_row)) return()
    rv$context_peak_flag <- flag_row[1, , drop = FALSE]
    pos <- as.integer(flag_row$Position[1])
    call <- as.character(flag_row$Call[1])
    comp <- as.character(flag_row$Competing_channel[1])
    r <- ensure_curation_state(selected_result())
    st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
    in_trim <- is.finite(st) && is.finite(en) && pos >= st && pos <= en
    trim_len <- if (in_trim) en - st + 1L else NA_integer_
    edge_zone <- if (in_trim) max(8L, min(40L, ceiling(trim_len * 0.15))) else NA_integer_
    recommendation <- ""
    if (in_trim) {
      dl <- pos - st; dr <- en - pos
      if (dl <= edge_zone || dr <= edge_zone) recommendation <- if (dl <= dr) "Recommended: trim the left edge through this position." else "Recommended: trim the right edge from this position."
    }
    comp_candidates <- unique(strsplit(gsub("[^ACGT/]", "", comp), "/", fixed = TRUE)[[1]])
    comp_candidates <- comp_candidates[comp_candidates %in% c("A","C","G","T")]
    base_choices <- unique(c(comp_candidates, c("A","C","G","T","N")))
    base_choices <- setNames(base_choices, base_choices)
    pair_code <- if (call %in% c("A","C","G","T") && length(comp_candidates) && comp_candidates[1] != call) iupac_for_pair(call, comp_candidates[1]) else NA_character_

    showModal(modalDialog(
      title = paste0("Review chromatogram position ", pos),
      easyClose = TRUE, size = "m",
      div(class = "about-callout",
          tags$strong(paste0(call, "  ↔  ", comp)), tags$br(),
          span(flag_evidence_text(flag_row))),
      if (nzchar(recommendation)) div(class = "status-note", recommendation) else NULL,
      selectInput("manual_base_choice", "Change current base to", choices = base_choices, selected = if (length(comp_candidates)) comp_candidates[1] else "N"),
      div(class = "curation-actions",
          actionButton("context_change_base_preview", "Review base change", class = "btn-primary"),
          if (!is.na(pair_code)) actionButton("context_iupac_preview", paste0("Set ambiguity ", pair_code)) else NULL,
          actionButton("context_mark_reviewed", "Keep call / mark reviewed")),
      if (in_trim) tagList(
        tags$hr(),
        div(class = "curation-actions",
            actionButton("context_trim_left", paste0("Trim left through ", pos)),
            actionButton("context_trim_right", paste0("Trim right from ", pos)))
      ) else NULL,
      footer = modalButton("Close")
    ))
  }

  observeEvent(input$peak_flag_context, {
    pos <- suppressWarnings(as.integer(input$peak_flag_context$position))
    df <- all_current_peak_flags()
    row <- df[df$Position == pos, , drop = FALSE]
    if (!nrow(row)) return()
    center_chromatogram_at(pos)
    show_peak_curation_menu(row[1, , drop = FALSE])
  }, ignoreInit = TRUE)

  show_curation_confirmation <- function(type, position, value = NULL) {
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(position)
    flag_row <- rv$context_peak_flag
    evidence <- flag_evidence_text(flag_row)
    if (type == "base") {
      calls <- curated_raw_calls(r); before <- if (pos >= 1 && pos <= length(calls)) calls[pos] else ""
      after <- toupper(as.character(value))
      title <- "Confirm base edit"
      body <- tagList(
        h4(paste0("Position ", pos, ": ", before, " → ", after)),
        p(evidence),
        div(class = "status-note", "The raw AB1 trace and original base call are preserved. This edit changes only the curated sequence and will be recorded in the audit log."))
    } else {
      st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
      side <- as.character(value)
      new_len <- if (side == "left") en - pos else pos - st
      title <- "Confirm manual trimming"
      body <- tagList(
        h4(if (side == "left") paste0("Remove positions ", st, "–", pos) else paste0("Remove positions ", pos, "–", en)),
        p(paste0("Current length: ", en - st + 1L, " bp · New length: ", max(0, new_len), " bp")),
        p(evidence),
        div(class = "status-note", "The automatic trimming boundaries remain stored separately. This manual boundary change is reversible and logged."))
    }
    rv$pending_curation <- list(type = type, position = pos, value = value, evidence = evidence)
    showModal(modalDialog(
      title = title, body,
      footer = tagList(modalButton("Cancel"), actionButton("confirm_curation_action", "Confirm change", class = "btn-danger")),
      easyClose = FALSE
    ))
  }

  observeEvent(input$context_change_base_preview, {
    req(rv$context_peak_flag, input$manual_base_choice)
    show_curation_confirmation("base", rv$context_peak_flag$Position[1], input$manual_base_choice)
  }, ignoreInit = TRUE)

  observeEvent(input$context_iupac_preview, {
    req(rv$context_peak_flag)
    call <- as.character(rv$context_peak_flag$Call[1])
    comp <- strsplit(as.character(rv$context_peak_flag$Competing_channel[1]), "/", fixed = TRUE)[[1]][1]
    code <- iupac_for_pair(call, comp)
    if (is.na(code)) { showNotification("No two-base IUPAC ambiguity code is available for this pair.", type = "warning"); return() }
    show_curation_confirmation("base", rv$context_peak_flag$Position[1], code)
  }, ignoreInit = TRUE)

  observeEvent(input$context_trim_left, {
    req(rv$context_peak_flag)
    show_curation_confirmation("trim", rv$context_peak_flag$Position[1], "left")
  }, ignoreInit = TRUE)

  observeEvent(input$context_trim_right, {
    req(rv$context_peak_flag)
    show_curation_confirmation("trim", rv$context_peak_flag$Position[1], "right")
  }, ignoreInit = TRUE)

  observeEvent(input$context_mark_reviewed, {
    req(rv$context_peak_flag, input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(rv$context_peak_flag$Position[1])
    snap <- curation_review_snapshot(r, pos)
    row <- data.frame(Action = "Mark reviewed / keep call", Position = pos,
                      Before = as.character(rv$context_peak_flag$Call[1]), After = as.character(rv$context_peak_flag$Call[1]),
                      Method = "Manual review", Evidence = flag_evidence_text(rv$context_peak_flag), Details = "No sequence change", stringsAsFactors = FALSE)
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), paste0("Reviewed position ", pos))
    commit_curated_result(input$inspect_sample, r2, paste0("Reviewed position ", pos))
    removeModal()
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_curation_action, {
    req(rv$pending_curation, input$inspect_sample)
    pnd <- rv$pending_curation
    r <- ensure_curation_state(selected_result())
    pos <- as.integer(pnd$position)
    calls <- curated_raw_calls(r)
    before_call <- if (pos >= 1 && pos <= length(calls)) calls[pos] else ""
    if (identical(pnd$type, "base")) {
      new_base <- toupper(as.character(pnd$value))
      if (identical(new_base, before_call)) {
        snap <- curation_review_snapshot(r, pos)
        row <- data.frame(Action = "Mark reviewed / keep call", Position = pos, Before = before_call, After = before_call,
                          Method = "Manual review", Evidence = pnd$evidence, Details = "Selected base equals current call; no sequence change", stringsAsFactors = FALSE)
        label <- paste0("Reviewed position ", pos)
        r2 <- curation_commit(r, snap, row, selected_processing_settings(), label)
        commit_curated_result(input$inspect_sample, r2, label)
        rv$pending_curation <- NULL
        removeModal(); showNotification(paste0(label, "."), type = "message")
        return()
      }
      snap <- curation_set_base_snapshot(r, pos, new_base)
      if (is.null(snap)) { showNotification("Could not apply the requested base edit.", type = "error"); return() }
      method <- if (new_base %in% c("R","Y","S","W","K","M")) "Manual IUPAC ambiguity" else "Manual base correction"
      row <- data.frame(Action = "Base edit", Position = pos, Before = before_call, After = new_base,
                        Method = method, Evidence = pnd$evidence, Details = "Confirmed by user", stringsAsFactors = FALSE)
      label <- paste0("Base edit at ", pos, ": ", before_call, "→", new_base)
    } else {
      side <- as.character(pnd$value)
      st <- as.integer(r$summary$trim_start[1]); en <- as.integer(r$summary$trim_end[1])
      snap <- curation_trim_snapshot(r, pos, side)
      if (is.null(snap)) { showNotification("The requested trimming boundary is not valid.", type = "error"); return() }
      row <- data.frame(Action = if (side == "left") "Manual trim left" else "Manual trim right", Position = pos,
                        Before = paste0(st, "-", en), After = paste0(snap$trim_start, "-", snap$trim_end),
                        Method = "Manual chromatogram curation", Evidence = pnd$evidence,
                        Details = "Flagged position excluded from retained sequence", stringsAsFactors = FALSE)
      label <- paste0(if (side == "left") "Trim left through " else "Trim right from ", pos)
    }
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), label)
    commit_curated_result(input$inspect_sample, r2, label)
    rv$pending_curation <- NULL
    removeModal()
    showNotification(paste0(label, "."), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$auto_correct_settings, {
    p <- ambiguous_peak_params_from_settings(selected_processing_settings())
    showModal(modalDialog(
      title = "Auto-correction criteria",
      div(class = "compact-help-note",
          span(icon("info-circle"), " These values control only automatic correction proposals. Flags and manual curation remain available even when a position does not meet these criteria.")),
      fluidRow(
        column(6,
          numericInput("auto_cfg_alt_called", "Alternative / current signal ≥", value = p$auto_min_alt_to_called, min = 1.01, max = 10, step = 0.05),
          numericInput("auto_cfg_alt_third", "Alternative / third channel ≥", value = p$auto_min_alt_to_third, min = 1.01, max = 10, step = 0.05)
        ),
        column(6,
          numericInput("auto_cfg_peak_offset", "Maximum peak offset (trace samples)", value = p$auto_max_peak_offset, min = 0, max = 10, step = 1),
          numericInput("auto_cfg_relative_signal", "Alternative signal / retained median ≥", value = p$auto_min_relative_signal, min = 0.05, max = 3, step = 0.05)
        )
      ),
      footer = tagList(
        actionButton("reset_auto_correct_criteria", "Reset defaults", icon = icon("refresh")),
        modalButton("Cancel"),
        actionButton("apply_auto_correct_criteria", "Apply criteria", class = "btn-primary")
      ),
      size = "m", easyClose = TRUE
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$reset_auto_correct_criteria, {
    d <- ambiguous_peak_defaults()
    updateNumericInput(session, "auto_cfg_alt_called", value = d$auto_min_alt_to_called)
    updateNumericInput(session, "auto_cfg_alt_third", value = d$auto_min_alt_to_third)
    updateNumericInput(session, "auto_cfg_peak_offset", value = d$auto_max_peak_offset)
    updateNumericInput(session, "auto_cfg_relative_signal", value = d$auto_min_relative_signal)
  }, ignoreInit = TRUE)

  observeEvent(input$apply_auto_correct_criteria, {
    vals <- c(
      alt_called = suppressWarnings(as.numeric(input$auto_cfg_alt_called)),
      alt_third = suppressWarnings(as.numeric(input$auto_cfg_alt_third)),
      peak_offset = suppressWarnings(as.numeric(input$auto_cfg_peak_offset)),
      relative_signal = suppressWarnings(as.numeric(input$auto_cfg_relative_signal))
    )
    if (any(!is.finite(vals)) || vals[["alt_called"]] <= 1 || vals[["alt_third"]] <= 1 ||
        vals[["peak_offset"]] < 0 || vals[["peak_offset"]] > 10 ||
        vals[["relative_signal"]] <= 0 || vals[["relative_signal"]] > 3) {
      showNotification("Check the auto-correction criteria values.", type = "error")
      return()
    }
    settings <- rv$settings
    settings$auto_correct_min_alt_to_called <- vals[["alt_called"]]
    settings$auto_correct_min_alt_to_third <- vals[["alt_third"]]
    settings$auto_correct_max_peak_offset <- as.integer(round(vals[["peak_offset"]]))
    settings$auto_correct_min_relative_signal <- vals[["relative_signal"]]
    rv$settings <- settings
    for (nm in names(rv$results)) {
      if (!is.list(rv$results[[nm]]$processing_settings)) rv$results[[nm]]$processing_settings <- settings
      rv$results[[nm]]$processing_settings$auto_correct_min_alt_to_called <- vals[["alt_called"]]
      rv$results[[nm]]$processing_settings$auto_correct_min_alt_to_third <- vals[["alt_third"]]
      rv$results[[nm]]$processing_settings$auto_correct_max_peak_offset <- as.integer(round(vals[["peak_offset"]]))
      rv$results[[nm]]$processing_settings$auto_correct_min_relative_signal <- vals[["relative_signal"]]
    }
    rv$auto_correct_preview_df <- data.frame()
    removeModal()
    showNotification("Auto-correction criteria updated.", type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$auto_correct_preview, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    df <- high_confidence_autocorrections(r, selected_processing_settings())
    if (nrow(df) && length(r$curation$reviewed_positions)) df <- df[!df$Position %in% r$curation$reviewed_positions, , drop = FALSE]
    if (!nrow(df)) {
      showNotification("No positions meet the current auto-correction criteria. Use Criteria to adjust the thresholds.", type = "message")
      return()
    }
    df$Proposed_base <- df$Competing_channel
    rv$auto_correct_preview_df <- df
    showModal(modalDialog(
      title = "Preview high-confidence base corrections",
      p(paste0(nrow(df), " position(s) meet all conservative auto-correction criteria. Nothing changes until you confirm.")),
      tableOutput("auto_correct_preview_table"),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_auto_correct", paste0("Apply ", nrow(df), " correction(s)"), class = "btn-danger")),
      size = "l", easyClose = FALSE
    ))
  }, ignoreInit = TRUE)

  output$auto_correct_preview_table <- renderTable({
    df <- rv$auto_correct_preview_df
    if (!is.data.frame(df) || !nrow(df)) return(NULL)
    out <- df[, c("Position","Call","Proposed_base","Alternative_to_called_ratio","Alternative_to_third_ratio","Competitor_peak_offset","Competitor_signal"), drop = FALSE]
    names(out) <- c("Position","Current","Proposed","Alternative / current","Alternative / third","Peak offset","Signal")
    out
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  observeEvent(input$confirm_auto_correct, {
    req(input$inspect_sample)
    df <- rv$auto_correct_preview_df
    if (!is.data.frame(df) || !nrow(df)) return()
    r <- ensure_curation_state(selected_result())
    temp <- r
    rows <- list()
    for (i in seq_len(nrow(df))) {
      pos <- as.integer(df$Position[i]); new_base <- as.character(df$Proposed_base[i])
      calls <- curated_raw_calls(temp); before <- calls[pos]
      snap_i <- curation_set_base_snapshot(temp, pos, new_base)
      if (is.null(snap_i)) next
      temp <- curation_restore_snapshot(temp, snap_i, selected_processing_settings())
      ev <- paste0("alternative/current=", df$Alternative_to_called_ratio[i],
                   "; alternative/third=", df$Alternative_to_third_ratio[i],
                   "; peak offset=", df$Competitor_peak_offset[i],
                   "; alternative signal=", df$Competitor_signal[i])
      rows[[length(rows)+1L]] <- data.frame(Action="Base edit", Position=pos, Before=before, After=new_base,
                                            Method="Auto high-confidence correction", Evidence=ev,
                                            Details="Applied after preview and user confirmation", stringsAsFactors=FALSE)
    }
    if (!length(rows)) { removeModal(); showNotification("No valid corrections remained to apply.", type="warning"); return() }
    action_rows <- do.call(rbind, rows)
    final_snap <- curation_snapshot(temp)
    label <- paste0("Auto-corrected ", nrow(action_rows), " high-confidence position(s)")
    r2 <- curation_commit(r, final_snap, action_rows, selected_processing_settings(), label)
    commit_curated_result(input$inspect_sample, r2, label)
    rv$auto_correct_preview_df <- data.frame()
    removeModal()
    showNotification(label, type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$curation_undo, {
    req(input$inspect_sample)
    ans <- curation_undo(selected_result(), selected_processing_settings())
    if (!isTRUE(ans$changed)) { showNotification("Nothing to undo for this sample.", type = "message"); return() }
    commit_curated_result(input$inspect_sample, ans$result, paste0("Undo: ", ans$label))
    showNotification(paste0("Undid: ", ans$label), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$curation_redo, {
    req(input$inspect_sample)
    ans <- curation_redo(selected_result(), selected_processing_settings())
    if (!isTRUE(ans$changed)) { showNotification("Nothing to redo for this sample.", type = "message"); return() }
    commit_curated_result(input$inspect_sample, ans$result, paste0("Redo: ", ans$label))
    showNotification(paste0("Redid: ", ans$label), type = "message")
  }, ignoreInit = TRUE)

  output$curation_history_table <- renderDT({
    r <- ensure_curation_state(selected_result())
    lg <- r$curation$audit_log
    if (!is.data.frame(lg) || !nrow(lg)) return(datatable(data.frame(Message="No manual curation actions for this sample."), rownames=FALSE, options=list(dom="t")))
    datatable(lg, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE, order=list(list(0,"desc"))))
  })

  observeEvent(input$curation_history, {
    req(input$inspect_sample)
    showModal(modalDialog(title = paste0("Curation history — ", input$inspect_sample), DTOutput("curation_history_table"), size="l", easyClose=TRUE, footer=modalButton("Close")))
  }, ignoreInit = TRUE)

  observeEvent(input$curation_reset, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    if (!length(r$curation$undo_stack) && !nrow(r$curation$base_edits) && identical(as.integer(r$curation$trim_start), as.integer(r$curation$auto_trim_start)) && identical(as.integer(r$curation$trim_end), as.integer(r$curation$auto_trim_end))) {
      showNotification("This sample is already at the automatic trimming/base-call state.", type="message"); return()
    }
    showModal(modalDialog(
      title="Reset sample to automatic processing?",
      p("This will restore the original automatic trim boundaries and remove all manual base edits/review marks from the active curated sequence. The reset itself is logged and can be undone."),
      footer=tagList(modalButton("Cancel"), actionButton("confirm_curation_reset", "Reset sample", class="btn-danger")), easyClose=FALSE
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_curation_reset, {
    req(input$inspect_sample)
    r <- ensure_curation_state(selected_result())
    snap <- curation_reset_snapshot(r)
    row <- data.frame(Action="Reset to automatic processing", Position=NA_integer_, Before=paste0(r$summary$trim_start[1], "-", r$summary$trim_end[1]),
                      After=paste0(r$curation$auto_trim_start, "-", r$curation$auto_trim_end), Method="Manual reset", Evidence="",
                      Details="Manual trims, base edits and review marks cleared from active curation state", stringsAsFactors=FALSE)
    r2 <- curation_commit(r, snap, row, selected_processing_settings(), "Reset to automatic processing")
    commit_curated_result(input$inspect_sample, r2, "Reset to automatic processing")
    removeModal(); showNotification("Sample reset to the automatic processing state.", type="message")
  }, ignoreInit = TRUE)

  output$chromatogram_plot <- plotly::renderPlotly({
    plot_result <- selected_result()
    plot_result$display_name <- qc_display_name(selected_sample_key())
    make_chromatogram_plotly(
      plot_result, selected_processing_settings(),
      flags = current_peak_flags(),
      show_flags = isTRUE(input$show_peak_flags)
    )
  })

  output$qc_plot <- renderPlot({
    draw_qc_metrics(selected_result(), selected_processing_settings())
  })

  # ---------------- Rename key ----------------
  rename_key_data <- reactive({
    req(input$rename_key_file)
    f <- input$rename_key_file; ext <- tolower(tools::file_ext(f$name))
    key <- if (ext=="xlsx") readxl::read_excel(f$datapath) else if (ext=="csv") read.csv(f$datapath,stringsAsFactors=FALSE,check.names=FALSE) else stop("Assignment key must be XLSX or CSV.")
    key <- as.data.frame(key,stringsAsFactors=FALSE)
    normalized <- tolower(gsub("[^A-Za-z0-9]+", "_", names(key)))
    names(key) <- normalized
    if ("gene" %in% names(key) && !"locus" %in% names(key)) names(key)[names(key) == "gene"] <- "locus"
    required <- c("old_id", "isolate", "locus", "direction")
    if (!all(required %in% names(key))) stop("Assignment key must contain old_id, isolate, locus (or gene), and direction columns.")
    key <- key[, required, drop = FALSE]
    for (nm in names(key)) key[[nm]] <- trimws(as.character(key[[nm]]))
    key$direction <- vapply(key$direction, stage2_normalize_direction, character(1))
    key[key$old_id != "", , drop = FALSE]
  })

  output$rename_key_status <- renderUI({
    if (is.null(input$rename_key_file)) return(div(class = "compact-hint", "Key columns: old_id, isolate, locus or gene, direction."))
    div(class = "compact-hint", paste("Selected:", input$rename_key_file$name))
  })
  output$download_assignment_key_template <- downloadHandler(
    filename = function() "PITAX_assignment_key_template.csv",
    content = function(file) file.copy("assignment_key_template.csv", file, overwrite = TRUE)
  )

  observeEvent(input$apply_rename_key, {
    req(nrow(rv$read_assignments))
    key <- tryCatch(rename_key_data(), error=function(e){showNotification(conditionMessage(e),type="error");NULL})
    if (is.null(key)) return()
    matched <- 0L
    for(i in seq_len(nrow(rv$read_assignments))) {
      original <- rv$read_assignments$Source_ID[i]
      idx <- which(key$old_id==original)
      if(!length(idx)) idx <- which(startsWith(original,key$old_id))
      if(length(idx)==1) {
        rv$read_assignments$Isolate[i] <- key$isolate[idx]
        rv$read_assignments$Locus[i] <- key$locus[idx]
        rv$read_assignments$Direction[i] <- key$direction[idx]
        matched <- matched+1L
      }
    }
    sync_assignment_state()
    showNotification(paste(matched,"of",nrow(rv$read_assignments),"reads matched; identity fields and generated names were updated."),type="message")
  })

  selected_assignment_rows <- function() {
    rows <- unique(c(input$assignment_upload_table_rows_selected, input$assignment_review_table_rows_selected))
    rows <- rows[!is.na(rows) & rows >= 1L & rows <= nrow(rv$read_assignments)]
    if (length(rows)) rows else seq_len(nrow(rv$read_assignments))
  }

  observeEvent(input$apply_assignment_batch, {
    req(nrow(rv$read_assignments))
    rows <- selected_assignment_rows()
    isolate_values <- as.character(rv$read_assignments$Isolate[rows])
    if (nzchar(input$batch_isolate_find)) isolate_values <- gsub(input$batch_isolate_find, input$batch_isolate_replace, isolate_values, fixed = TRUE)
    if (nzchar(input$batch_isolate_prefix)) isolate_values <- paste0(input$batch_isolate_prefix, isolate_values)
    if (nzchar(input$batch_isolate_suffix)) isolate_values <- paste0(isolate_values, input$batch_isolate_suffix)
    rv$read_assignments$Isolate[rows] <- isolate_values
    if (nzchar(trimws(input$batch_locus))) rv$read_assignments$Locus[rows] <- trimws(input$batch_locus)
    if (nzchar(input$batch_direction)) rv$read_assignments$Direction[rows] <- input$batch_direction
    sync_assignment_state()
    showNotification(paste("Updated", length(rows), "read assignment(s)."), type = "message")
  })

  observeEvent(input$reset_assignments, {
    req(nrow(rv$read_assignments))
    rows <- selected_assignment_rows()
    rv$read_assignments$Isolate[rows] <- ""
    rv$read_assignments$Locus[rows] <- ""
    rv$read_assignments$Direction[rows] <- "Unknown"
    rv$read_assignments$Primer[rows] <- ""
    sync_assignment_state()
  })

  assignment_table_widget <- function() {
    req(nrow(rv$read_assignments))
    df <- rv$read_assignments[, c("Source_ID", "Isolate", "Locus", "Direction", "Final_Name"), drop = FALSE]
    names(df) <- c("Upload barcode / source", "Isolate", "Gene / locus", "Direction", "Final read / FASTA name")
    datatable(df, rownames = FALSE, selection = "multiple",
              editable = list(target = "cell", disable = list(columns = c(0, 4))),
              options = list(pageLength = 20, scrollX = TRUE, dom = "tip"))
  }
  output$assignment_upload_table <- renderDT(assignment_table_widget())
  output$assignment_review_table <- renderDT(assignment_table_widget())

  apply_assignment_cell_edit <- function(info) {
    if (is.null(info) || !info$col %in% 1:3) return(invisible(NULL))
    field <- c("Isolate", "Locus", "Direction")[[info$col]]
    value <- as.character(info$value)
    if (field == "Direction") value <- stage2_normalize_direction(value)
    rv$read_assignments[info$row, field] <- value
    sync_assignment_state()
  }
  observeEvent(input$assignment_upload_table_cell_edit, apply_assignment_cell_edit(input$assignment_upload_table_cell_edit))
  observeEvent(input$assignment_review_table_cell_edit, apply_assignment_cell_edit(input$assignment_review_table_cell_edit))

  rename_error <- reactive({
    stage2_identity_error(rv$read_assignments)
  })
  assignment_validation_ui <- function() {
    e <- stage2_identity_error(rv$read_assignments)
    if(is.null(e)) div(class="status-ok","✓ Explicit identity fields are complete; final read names were generated by PITAX.") else div(class="status-error",paste0("⚠ ",e))
  }
  output$rename_validation <- renderUI(assignment_validation_ui())
  output$upload_assignment_validation <- renderUI(assignment_validation_ui())

  # ---------------- Export records ----------------
  export_records <- reactive({
    req(rv$results,rv$rename)
    records <- list()
    for(original_name in names(rv$results)) {
      r <- rv$results[[original_name]]
      idx <- match(original_name,rv$rename$Original_name)
      r$final_name <- if(!is.na(idx)) rv$rename$New_name[idx] else original_name
      records[[original_name]] <- r
    }
    records
  })
  export_summary_df <- reactive({
    req(rv$summary,rv$rename)
    df <- rv$summary
    df$final_name <- rv$rename$New_name[match(df$sample_id,rv$rename$Original_name)]
    df
  })

  observeEvent(input$to_export, {
    e <- rename_error(); if(!is.null(e)){showNotification(e,type="error");return()}
    updateTabsetPanel(session,"pipeline_step",selected="export")
  })

  project_export_stem <- function() {
    if (is.list(rv$architecture) && is.data.frame(rv$architecture$loci) && nrow(rv$architecture$loci) > 1L) return("PITAX_multi_locus")
    if (!is.null(rv$settings)) return(clean_fasta_name(stage2_scalar_text(rv$settings$target, "PITAX_project")))
    "PITAX_project"
  }

  output$export_summary <- renderUI({
    req(rv$summary,rv$settings)
    loci <- unique(as.character(rv$summary$target))
    architecture_summary <- stage2_architecture_summary(rv$architecture)
    tagList(
      p(strong("Locus/loci: "), paste(loci, collapse = ", ")),
      p(strong("Architecture: "), paste0(architecture_summary$Isolates, " isolate(s), ", architecture_summary$Loci, " locus/loci, ", architecture_summary$Reads, " read(s)")),
      p(strong("Forward primer: "),ifelse(rv$settings$forward_primer=="","Not specified",rv$settings$forward_primer)),
      p(strong("Reverse primer: "),ifelse(rv$settings$reverse_primer=="","Not specified",rv$settings$reverse_primer)),
      p(strong("Samples processed: "),nrow(rv$summary)),
      p(strong("Sequences available for export: "),sum(rv$summary$trimmed_length>0,na.rm=TRUE))
    )
  })

  # ---------------- Checkpoints ----------------
  output$download_trim_checkpoint <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_checkpoint_B_qc.zip"),
    content=function(file){
      write_checkpoint_zip(file,"qc",export_records(),export_summary_df(),rv$settings,rv$rename,results=rv$results,
                           read_assignments=rv$read_assignments,architecture=rv$architecture)
    })
  output$download_rename_checkpoint <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_checkpoint_A_rename.zip"),
    content=function(file) {
      req(is.null(stage2_identity_error(rv$read_assignments)), is.list(rv$architecture), is.list(rv$settings))
      write_assignment_checkpoint_zip(file, rv$rename, rv$read_assignments, rv$architecture, rv$settings)
    })

  output$download_blast_fasta <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_BLAST.fasta"),
    content=function(file) writeLines(make_fasta(export_records(),FALSE),file))
  output$download_full_fasta <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_trimmed_sequences.fasta"),
    content=function(file) writeLines(make_fasta(export_records(),TRUE,export_summary_df()),file))
  output$download_summary_csv <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_trim_summary.csv"),
    content=function(file) write.csv(export_summary_df(),file,row.names=FALSE,fileEncoding="UTF-8"))
  output$download_all_zip <- downloadHandler(
    filename=function() paste0(project_export_stem(),"_Sanger_pipeline_results.zip"),
    content=function(file) write_checkpoint_zip(file,"final",export_records(),export_summary_df(),rv$settings,rv$rename,results=rv$results,
                                                read_assignments=rv$read_assignments,architecture=rv$architecture))

  # ---------------- BLAST workspace ----------------
  resolve_blast_original_name <- function(value, records = export_records()) {
    if (is.null(value) || !length(value) || !length(records)) return(NULL)
    value <- as.character(value)[1]
    if (!nzchar(value)) return(NULL)
    if (value %in% names(records)) return(value)
    final_names <- vapply(records, function(x) as.character(x$final_name), character(1))
    idx <- which(final_names == value)
    if (length(idx)) return(names(records)[idx[1]])
    NULL
  }

  # Keep the BLAST sequence selector synchronized with the current processed/renamed
  # record set without resetting the user's selection every time the tab is opened.
  observe({
    if (!length(rv$results) || is.null(rv$rename)) {
      updateSelectInput(session, "blast_sample", choices = character(), selected = character())
      return()
    }
    rec <- export_records()
    if (!length(rec)) {
      updateSelectInput(session, "blast_sample", choices = character(), selected = character())
      return()
    }
    final_names <- vapply(rec, function(x) as.character(x$final_name), character(1))
    choices <- setNames(names(rec), final_names)
    current <- isolate(input$blast_sample)
    selected <- resolve_blast_original_name(current, rec)
    if (is.null(selected) && length(rec)) selected <- names(rec)[1]
    updateSelectInput(session, "blast_sample", choices=choices, selected=if (is.null(selected)) character() else selected)
  })

  observeEvent(input$to_blast, {
    updateTabsetPanel(session, "pipeline_step", selected="blast")
  })

  blast_selected <- reactive({
    rec <- export_records()
    original_name <- resolve_blast_original_name(input$blast_sample, rec)
    req(!is.null(original_name))
    rec[[original_name]]
  })

  output$blast_sequence_preview_ui <- renderUI({
    r <- blast_selected()
    seq_text <- if (is.null(r$seq)) "" else as.character(r$seq)
    if (!nzchar(seq_text)) {
      return(tags$pre(class="sequence-code-viewer sequence-code-empty", "No processed sequence is available for this selection."))
    }
    tags$pre(class="sequence-code-viewer", wrap_sequence(seq_text, 80))
  })

  observeEvent(input$copy_blast_sequence,
               session$sendCustomMessage("copyText", list(text=blast_selected()$seq)))
  observeEvent(input$open_ncbi_blast,
               session$sendCustomMessage("openUrl", list(url="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn")))

  output$download_selected_blast <- downloadHandler(
    filename=function() paste0(clean_fasta_name(blast_selected()$final_name), ".fasta"),
    content=function(file) writeLines(
      paste0(">", clean_fasta_name(blast_selected()$final_name), "\n", wrap_sequence(blast_selected()$seq, 80)),
      file
    )
  )

  # NCBI BLAST is a shared service. Keep all automated BLAST contacts at least
  # 10 seconds apart, in addition to the per-RID one-minute polling rule.
  wait_for_ncbi_contact_slot <- function(min_seconds = 10) {
    last <- rv$ncbi_last_contact
    if (!is.null(last) && length(last) && !is.na(last)) {
      elapsed <- as.numeric(difftime(Sys.time(), last, units="secs"))
      if (is.finite(elapsed) && elapsed < min_seconds) {
        Sys.sleep(min_seconds - elapsed)
      }
    }
    rv$ncbi_last_contact <- Sys.time()
  }

  parse_submit_response <- function(txt) {
    rid_match <- regmatches(txt, regexpr("RID = [A-Z0-9-]+", txt))
    rtoe_match <- regmatches(txt, regexpr("RTOE = [0-9]+", txt))
    rid <- if (length(rid_match) && nzchar(rid_match)) sub("RID = ", "", rid_match, fixed=TRUE) else ""
    rtoe <- if (length(rtoe_match) && nzchar(rtoe_match)) sub("RTOE = ", "", rtoe_match, fixed=TRUE) else NA_character_
    list(rid=rid, rtoe=rtoe)
  }

  latest_job_index_for_sample <- function(original_name) {
    idx <- which(rv$blast_jobs$original_name == original_name)
    if (!length(idx)) return(NA_integer_)
    tail(idx, 1)
  }

  matching_active_blast_job_index <- function(original_name, database, hitlist_size) {
    if (!is.data.frame(rv$blast_jobs) || !nrow(rv$blast_jobs)) return(NA_integer_)
    idx <- which(
      rv$blast_jobs$original_name == original_name &
      rv$blast_jobs$status %in% c("SUBMITTED", "WAITING", "READY") &
      nzchar(as.character(rv$blast_jobs$database)) &
      as.character(rv$blast_jobs$database) == as.character(database) &
      suppressWarnings(as.integer(rv$blast_jobs$hitlist_size)) == as.integer(hitlist_size)
    )
    if (!length(idx)) return(NA_integer_)
    tail(idx, 1)
  }

  submit_one_blast <- function(original_name) {
    records <- export_records()
    if (!original_name %in% names(records)) {
      return(list(ok=FALSE, message="Sequence is not available in the processed record set."))
    }
    r <- records[[original_name]]
    if (is.null(r) || !nzchar(r$seq)) {
      return(list(ok=FALSE, message="Processed sequence is empty."))
    }

    hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
    wait_for_ncbi_contact_slot(10)

    req_obj <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
      httr2::req_body_form(
        CMD="Put", PROGRAM="blastn", DATABASE=input$blast_database,
        QUERY=paste0(">", clean_fasta_name(r$final_name), "\n", r$seq),
        HITLIST_SIZE=hitlist
      ) |>
      httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")

    txt <- tryCatch(
      httr2::resp_body_string(httr2::req_perform(req_obj)),
      error=function(e) structure(list(error=conditionMessage(e)), class="blast_submit_error")
    )
    if (inherits(txt, "blast_submit_error")) {
      return(list(ok=FALSE, message=txt$error))
    }

    parsed <- parse_submit_response(txt)
    if (!nzchar(parsed$rid)) {
      return(list(ok=FALSE, message="NCBI did not return a BLAST RID."))
    }

    rv$blast_jobs <- rbind(rv$blast_jobs, data.frame(
      final_name=r$final_name,
      original_name=original_name,
      rid=parsed$rid,
      rtoe=parsed$rtoe,
      database=as.character(input$blast_database),
      hitlist_size=hitlist,
      status="SUBMITTED",
      submitted_at=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      last_checked_at="",
      stringsAsFactors=FALSE
    ))

    list(ok=TRUE, rid=parsed$rid, rtoe=parsed$rtoe, final_name=r$final_name)
  }

  observeEvent(input$submit_ncbi_blast, {
    rec <- export_records()
    original_name <- resolve_blast_original_name(input$blast_sample, rec)
    req(!is.null(original_name))
    requested_hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
    existing_idx <- matching_active_blast_job_index(original_name, input$blast_database, requested_hitlist)
    if (!is.na(existing_idx)) {
      existing <- rv$blast_jobs[existing_idx, , drop=FALSE]
      showNotification(
        paste0("Matching BLAST job already exists (RID ", existing$rid[1], ", ",
               existing$database[1], ", ", existing$hitlist_size[1], " hits; status ",
               existing$status[1], "). No new request was submitted."),
        type="message", duration=8
      )
      return()
    }
    ans <- submit_one_blast(original_name)
    if (!isTRUE(ans$ok)) {
      showNotification(ans$message, type="error")
      return()
    }
    showNotification(paste("Submitted to NCBI. RID:", ans$rid), type="message")
  })

  observeEvent(input$submit_all_ncbi_blast, {
    records <- export_records()
    sample_names <- names(records)[vapply(records, function(x) !is.null(x$seq) && nzchar(x$seq), logical(1))]
    if (!length(sample_names)) {
      showNotification("No processed sequences are available for BLAST submission.", type="warning")
      return()
    }

    submitted <- 0L
    skipped <- 0L
    failed <- 0L
    failures <- character()
    rv$blast_batch_status_text <- paste0("Submitting ", length(sample_names), " processed sequence(s) to NCBI...")

    withProgress(message="Submitting all sequences to NCBI BLAST", value=0, {
      for (i in seq_along(sample_names)) {
        nm <- sample_names[i]
        requested_hitlist <- min(100L, max(1L, as.integer(input$blast_hitlist)))
        idx <- matching_active_blast_job_index(nm, input$blast_database, requested_hitlist)
        if (!is.na(idx)) {
          skipped <- skipped + 1L
          incProgress(
            1/length(sample_names),
            detail=paste("Skipping matching active job:", records[[nm]]$final_name,
                         paste0("(", input$blast_database, ", ", requested_hitlist, " hits)"))
          )
          next
        }

        incProgress(0, detail=paste("Submitting", records[[nm]]$final_name))
        ans <- submit_one_blast(nm)
        if (isTRUE(ans$ok)) {
          submitted <- submitted + 1L
        } else {
          failed <- failed + 1L
          failures <- c(failures, paste0(records[[nm]]$final_name, ": ", ans$message))
        }
        incProgress(1/length(sample_names))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Batch submission complete · submitted: ", submitted,
      " · skipped matching active jobs: ", skipped,
      " · failed: ", failed,
      if (length(failures)) paste0(" · ", paste(failures, collapse=" | ")) else ""
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  store_retrieved_hits <- function(hits, r, original_name, rid) {
    # v2.9 BLAST parser: one biological/database hit per accession. Multiple
    # HSPs (local alignment segments) must not receive independent weight in
    # taxonomy or sequence-evidence calculations.
    hits <- normalize_blast_hits_unique_accession(hits)
    hits$final_name <- r$final_name
    hits$original_name <- original_name
    hits$rid <- rid

    hit_cols <- c(
      "final_name","original_name","rid","rank","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","alignment_length","hsp_count","match_support"
    )
    for (nm in setdiff(hit_cols, names(hits))) hits[[nm]] <- NA
    hits <- hits[, hit_cols, drop=FALSE]

    if (nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits)) {
      rv$blast_hits <- rv$blast_hits[rv$blast_hits$rid != rid, , drop=FALSE]
    }
    rv$blast_hits <- rbind(rv$blast_hits, hits)

    top <- hits[1, , drop=FALSE]
    top_cols <- c(
      "final_name","original_name","rid","organism","record_title","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","match_support"
    )

    # Preliminary identification is one current top-hit row per sequence.
    if (nrow(rv$blast_ids) && "original_name" %in% names(rv$blast_ids)) {
      rv$blast_ids <- rv$blast_ids[rv$blast_ids$original_name != original_name, , drop=FALSE]
    }
    rv$blast_ids <- rbind(rv$blast_ids, top[, top_cols, drop=FALSE])
  }

  retrieve_blast_job <- function(job_idx) {
    if (is.na(job_idx) || job_idx < 1 || job_idx > nrow(rv$blast_jobs)) {
      return(list(status="ERROR", contacted=FALSE, message="BLAST job index is invalid."))
    }
    if (identical(as.character(rv$blast_jobs$status[job_idx]), "STALE")) {
      return(list(status="STALE", contacted=FALSE, message="This RID belongs to a sequence version that was changed during manual curation. Submit the current curated sequence as a new BLAST job."))
    }

    original_name <- rv$blast_jobs$original_name[job_idx]
    records <- export_records()
    if (!original_name %in% names(records)) {
      return(list(status="ERROR", contacted=FALSE, message="Processed sequence is no longer available."))
    }
    r <- records[[original_name]]
    rid <- rv$blast_jobs$rid[job_idx]

    # Do not poll a single RID more often than once per minute.
    reference_time <- rv$blast_jobs$last_checked_at[job_idx]
    if (!nzchar(reference_time)) reference_time <- rv$blast_jobs$submitted_at[job_idx]
    ref <- suppressWarnings(as.POSIXct(reference_time, format="%Y-%m-%d %H:%M:%S"))
    if (!is.na(ref)) {
      elapsed <- as.numeric(difftime(Sys.time(), ref, units="secs"))
      if (is.finite(elapsed) && elapsed < 60) {
        return(list(
          status="TOO_SOON", contacted=FALSE,
          message=paste0("Wait ", ceiling(60-elapsed), " more seconds before checking RID ", rid, ".")
        ))
      }
    }

    wait_for_ncbi_contact_slot(10)
    rv$blast_jobs$last_checked_at[job_idx] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    req_obj <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
      httr2::req_url_query(CMD="Get", RID=rid, FORMAT_TYPE="XML2") |>
      httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")

    txt <- tryCatch(
      httr2::resp_body_string(httr2::req_perform(req_obj)),
      error=function(e) structure(list(error=conditionMessage(e)), class="blast_retrieve_error")
    )
    if (inherits(txt, "blast_retrieve_error")) {
      return(list(status="ERROR", contacted=TRUE, message=txt$error))
    }

    if (grepl("Status=WAITING", txt)) {
      rv$blast_jobs$status[job_idx] <- "WAITING"
      return(list(status="WAITING", contacted=TRUE, message="NCBI BLAST is still running."))
    }
    if (grepl("Status=FAILED|Status=UNKNOWN", txt)) {
      rv$blast_jobs$status[job_idx] <- "FAILED/UNKNOWN"
      return(list(status="FAILED", contacted=TRUE, message="NCBI reports failed or unknown job."))
    }

    rv$blast_jobs$status[job_idx] <- "READY"
    rv$blast_raw[[rid]] <- txt

    requested_hits <- suppressWarnings(as.integer(rv$blast_jobs$hitlist_size[job_idx]))
    if (!is.finite(requested_hits) || requested_hits < 1) requested_hits <- as.integer(input$blast_hitlist)
    requested_hits <- min(100L, max(1L, requested_hits))

    hits <- parse_blast_xml2_hits(
      txt,
      query_length=nchar(r$seq),
      max_hits=requested_hits
    )

    if (!nrow(hits)) {
      # Defensive fallback to NCBI's headerless 12-field tabular CSV.
      wait_for_ncbi_contact_slot(10)
      csv_req <- httr2::request("https://blast.ncbi.nlm.nih.gov/Blast.cgi") |>
        httr2::req_url_query(CMD="Get", RID=rid, FORMAT_TYPE="CSV", ALIGNMENT_VIEW="Tabular") |>
        httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI BLAST client")
      csv_txt <- tryCatch(
        httr2::resp_body_string(httr2::req_perform(csv_req)),
        error=function(e) NULL
      )
      if (!is.null(csv_txt)) {
        hits <- parse_blast_csv_hits_fallback(csv_txt, query_length=nchar(r$seq))
        if (nrow(hits)) {
          hits <- head(hits, requested_hits)
          rv$blast_raw[[rid]] <- csv_txt
        }
      }
    }

    hits <- normalize_blast_hits_unique_accession(hits)
    if (nrow(hits) > requested_hits) hits <- hits[seq_len(requested_hits), , drop = FALSE]

    if (!nrow(hits)) {
      return(list(status="NO_HITS", contacted=TRUE, message="BLAST is ready, but no unique accession-level hits could be parsed."))
    }

    # Enrich every unique accession whose XML/tabular result lacks display metadata.
    needs_meta <- which(
      (!nzchar(ifelse(is.na(hits$organism), "", hits$organism))) |
      (!nzchar(ifelse(is.na(hits$record_title), "", hits$record_title)))
    )

    if (length(needs_meta)) {
      accessions_needed <- unique(hits$accession[needs_meta])
      meta_df <- tryCatch(
        fetch_ncbi_nucleotide_metadata_batch(accessions_needed),
        error=function(e) data.frame()
      )

      if (nrow(meta_df)) {
        strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))
        for (i in needs_meta) {
          acc <- hits$accession[i]
          idx <- which(meta_df$accession == acc)
          if (!length(idx) && "primary_accession" %in% names(meta_df)) idx <- which(meta_df$primary_accession == strip_version(acc))
          if (!length(idx)) idx <- which(strip_version(meta_df$accession) == strip_version(acc))
          if (length(idx)) {
            m <- meta_df[idx[1], , drop=FALSE]
            if ((!nzchar(hits$organism[i]) || is.na(hits$organism[i])) && nzchar(m$organism[1])) hits$organism[i] <- m$organism[1]
            if ((!nzchar(hits$record_title[i]) || is.na(hits$record_title[i])) && nzchar(m$title[1])) hits$record_title[i] <- m$title[1]
            if ((is.na(hits$taxid[i]) || !is.finite(hits$taxid[i])) && !is.na(m$taxid[1])) hits$taxid[i] <- m$taxid[1]
          }
        }
      }
    }

    store_retrieved_hits(hits, r, original_name, rid)
    list(status="READY", contacted=TRUE, hits=nrow(hits), message=paste0(nrow(hits), " hit(s) retrieved."))
  }

  observeEvent(input$retrieve_ncbi_blast, {
    selected_rows <- input$blast_jobs_table_rows_selected
    selected_rows <- suppressWarnings(as.integer(selected_rows))
    selected_rows <- selected_rows[is.finite(selected_rows) & selected_rows >= 1L & selected_rows <= nrow(rv$blast_jobs)]

    # Backward-compatible fallback: if no job-table rows are selected, retrieve
    # the newest job for the sequence currently shown in Query workspace.
    if (!length(selected_rows)) {
      rec <- export_records()
      original_name <- resolve_blast_original_name(input$blast_sample, rec)
      if (is.null(original_name)) {
        showNotification("Select one or more BLAST jobs, or choose a sequence in Query workspace.", type="warning")
        return()
      }
      job_idx <- latest_job_index_for_sample(original_name)
      if (is.na(job_idx)) {
        showNotification("No BLAST job exists for the selected sequence.", type="warning")
        return()
      }
      selected_rows <- job_idx
    }

    # Preserve table order and avoid contacting the same RID twice if the UI ever
    # reports duplicate row selections.
    selected_rows <- unique(selected_rows)
    ready <- 0L; waiting <- 0L; too_soon <- 0L; failed <- 0L; no_hits <- 0L; stale <- 0L

    withProgress(message="Retrieving selected NCBI BLAST jobs", value=0, {
      for (k in seq_along(selected_rows)) {
        idx <- selected_rows[k]
        rid <- rv$blast_jobs$rid[idx]
        final_name <- rv$blast_jobs$final_name[idx]
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(as.character(rv$blast_jobs$status[idx]), "READY") && already_stored) {
          ans <- list(status="READY", message="Already retrieved.")
        } else {
          incProgress(0, detail=paste("Checking", final_name))
          ans <- retrieve_blast_job(idx)
        }

        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status == "WAITING") waiting <- waiting + 1L
        else if (ans$status == "TOO_SOON") too_soon <- too_soon + 1L
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(1/length(selected_rows))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Selected retrieval · jobs: ", length(selected_rows),
      " · ready: ", ready,
      " · still running: ", waiting,
      " · too soon to poll: ", too_soon,
      " · no parsed hits: ", no_hits,
      " · stale: ", stale,
      " · failed: ", failed
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  observeEvent(input$retrieve_all_ncbi_blast, {
    if (!nrow(rv$blast_jobs)) {
      showNotification("No BLAST jobs have been submitted yet.", type="warning")
      return()
    }

    # Only the newest RID for each sequence is relevant to the current workspace.
    latest_indices <- unname(vapply(unique(rv$blast_jobs$original_name), latest_job_index_for_sample, integer(1)))
    latest_indices <- latest_indices[is.finite(latest_indices)]
    if (!length(latest_indices)) return()

    ready <- 0L; waiting <- 0L; too_soon <- 0L; failed <- 0L; no_hits <- 0L; stale <- 0L
    rv$blast_batch_status_text <- paste0("Checking ", length(latest_indices), " BLAST job(s)...")

    withProgress(message="Retrieving submitted NCBI BLAST jobs", value=0, {
      for (k in seq_along(latest_indices)) {
        idx <- latest_indices[k]
        final_name <- rv$blast_jobs$final_name[idx]
        rid <- rv$blast_jobs$rid[idx]

        # A READY RID that is already present in the hit store does not need another server request.
        already_stored <- nrow(rv$blast_hits) && "rid" %in% names(rv$blast_hits) && rid %in% rv$blast_hits$rid
        if (identical(rv$blast_jobs$status[idx], "READY") && already_stored) {
          ready <- ready + 1L
          incProgress(1/length(latest_indices), detail=paste(final_name, "already retrieved"))
          next
        }

        incProgress(0, detail=paste("Checking", final_name))
        ans <- retrieve_blast_job(idx)
        if (ans$status == "READY") ready <- ready + 1L
        else if (ans$status == "WAITING") waiting <- waiting + 1L
        else if (ans$status == "TOO_SOON") too_soon <- too_soon + 1L
        else if (ans$status == "NO_HITS") no_hits <- no_hits + 1L
        else if (ans$status == "STALE") stale <- stale + 1L
        else failed <- failed + 1L
        incProgress(1/length(latest_indices))
      }
    })

    rv$blast_batch_status_text <- paste0(
      "Batch retrieval · ready: ", ready,
      " · still running: ", waiting,
      " · too soon to poll: ", too_soon,
      " · no parsed hits: ", no_hits,
      " · stale: ", stale,
      " · failed: ", failed
    )
    showNotification(rv$blast_batch_status_text, type=if(failed) "warning" else "message", duration=10)
  })

  output$blast_batch_status <- renderUI({
    div(class="tax-note", strong("Batch status: "), rv$blast_batch_status_text)
  })

  output$blast_job_status <- renderUI({
    req(input$blast_sample)
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
    if (!nrow(jobs)) return(p(class="settings-note", "No NCBI submission yet for this sequence."))
    j <- jobs[nrow(jobs), ]
    div(
      class=if(j$status == "READY") "status-ok" else "status-warning",
      paste("Selected sequence · RID:", j$rid, "| Status:", j$status,
            "| Database:", if (nzchar(as.character(j$database))) j$database else "legacy/unknown",
            "| Requested hits:", j$hitlist_size,
            "| Estimated wait (RTOE):", j$rtoe, "seconds")
    )
  })

  output$blast_jobs_table <- renderDT({
    df <- rv$blast_jobs
    if (!nrow(df)) return(datatable(data.frame(Message="No BLAST jobs submitted yet."), rownames=FALSE, options=list(dom="t")))
    keep <- c("final_name","original_name","rid","database","hitlist_size","rtoe","status","submitted_at","last_checked_at")
    df <- df[, intersect(keep, names(df)), drop=FALSE]
    friendly <- c(
      final_name="Sample", original_name="Original sample", rid="RID", database="Database", hitlist_size="Hits requested",
      rtoe="Estimated wait (s)", status="Status", submitted_at="Submitted at", last_checked_at="Last checked at"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(df, rownames=FALSE, selection=list(mode="multiple", target="row"), options=list(pageLength=25, scrollX=TRUE))
  })

  output$blast_identification_table <- renderDT({
    # Deliberately show all retrieved sequences here: one current top hit per
    # sequence. The selected-sequence filter belongs to the detailed hit table.
    df <- rv$blast_ids
    if (!nrow(df)) {
      return(datatable(
        data.frame(Message="No parsed identification results yet. Submit and retrieve one or more BLAST jobs."),
        rownames=FALSE, options=list(dom="t")
      ))
    }

    df <- df[order(tolower(as.character(df$final_name))), , drop=FALSE]
    keep <- intersect(c(
      "final_name","organism","accession","identity_percent","query_coverage_percent",
      "evalue","bit_score","match_support","record_title","taxid","rid"
    ), names(df))
    df <- df[, keep, drop=FALSE]
    friendly <- c(
      final_name="Sample", organism="Organism / taxon", accession="Accession",
      identity_percent="Identity (%)", query_coverage_percent="Query coverage (%)",
      evalue="E-value", bit_score="Bit score", match_support="Match support",
      record_title="NCBI hit title", taxid="NCBI TaxID", rid="RID"
    )
    names(df) <- unname(friendly[names(df)])

    datatable(
      df, rownames=FALSE, filter="top",
      options=list(
        pageLength=25, scrollX=TRUE, autoWidth=TRUE,
        columnDefs=list(
          list(targets=1, className="dt-organism", width="190px"),
          list(targets=8, className="dt-hit-title", width="440px"),
          list(targets=c(0,2,3,4,5,6,7,9,10), className="dt-nowrap")
        )
      )
    )
  })

  output$blast_hits_table <- renderDT({
    df <- rv$blast_hits
    if (nrow(df) && !is.null(input$blast_sample)) {
      jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
      if (nrow(jobs)) {
        latest_rid <- jobs$rid[nrow(jobs)]
        df <- df[df$rid == latest_rid, , drop=FALSE]
      } else {
        df <- df[0, , drop=FALSE]
      }
    }
    if (!nrow(df)) return(datatable(data.frame(Message="No retrieved hits yet for the selected sequence."), rownames=FALSE, options=list(dom="t")))

    keep <- intersect(c(
      "rank","organism","accession","identity_percent","query_coverage_percent",
      "evalue","bit_score","hsp_count","match_support","record_title","taxid"
    ), names(df))
    df <- df[, keep, drop=FALSE]
    friendly <- c(
      rank="Rank", organism="Organism / taxon", accession="Accession",
      identity_percent="Identity (%)", query_coverage_percent="Query coverage (%)",
      evalue="E-value", bit_score="Bit score", hsp_count="HSP count", match_support="Match support",
      record_title="NCBI hit title", taxid="NCBI TaxID"
    )
    names(df) <- unname(friendly[names(df)])

    datatable(
      df, rownames=FALSE, filter="top",
      options=list(
        pageLength=25, scrollX=TRUE, autoWidth=TRUE,
        columnDefs=list(
          list(targets=1, className="dt-organism", width="190px"),
          list(targets=9, className="dt-hit-title", width="440px"),
          list(targets=c(0,2,3,4,5,6,7,8,10), className="dt-nowrap")
        )
      )
    )
  })

  output$blast_raw_preview <- renderText({
    req(input$blast_sample)
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == input$blast_sample, , drop=FALSE]
    if (!nrow(jobs)) return("")
    rid <- jobs$rid[nrow(jobs)]
    txt <- rv$blast_raw[[rid]]
    if (is.null(txt)) return("")
    substr(txt, 1, min(5000, nchar(txt)))
  })

  output$download_blast_jobs <- downloadHandler(
    filename=function() "NCBI_BLAST_jobs_and_identifications.csv",
    content=function(file) {
      jobs <- rv$blast_jobs
      if (nrow(rv$blast_ids)) jobs <- merge(jobs, rv$blast_ids, by=c("final_name","original_name","rid"), all.x=TRUE)
      write.csv(jobs, file, row.names=FALSE, fileEncoding="UTF-8")
    }
  )

  output$download_blast_hits <- downloadHandler(
    filename=function() "NCBI_BLAST_all_hits.csv",
    content=function(file) write.csv(rv$blast_hits, file, row.names=FALSE, fileEncoding="UTF-8")
  )

  # ---------------- Taxonomic interpretation ----------------
  latest_blast_rid_for_sample <- function(original_name) {
    jobs <- rv$blast_jobs[rv$blast_jobs$original_name == original_name,,drop=FALSE]
    if (!nrow(jobs)) return("")
    jobs$rid[nrow(jobs)]
  }

  blast_hits_for_sample <- function(original_name) {
    rid <- latest_blast_rid_for_sample(original_name)
    if (!nzchar(rid) || !nrow(rv$blast_hits)) return(rv$blast_hits[0,,drop=FALSE])
    rv$blast_hits[rv$blast_hits$rid == rid,,drop=FALSE]
  }

  update_tax_sample_choices <- function(selected = NULL) {
    if (!nrow(rv$blast_hits)) {
      updateSelectInput(session, "tax_sample", choices = character())
      return(invisible(NULL))
    }
    pairs <- unique(rv$blast_hits[,c("original_name","final_name"),drop=FALSE])
    pairs <- pairs[nzchar(as.character(pairs$original_name)),,drop=FALSE]
    choices <- setNames(as.character(pairs$original_name), as.character(pairs$final_name))
    if (is.null(selected) || !selected %in% choices) selected <- if (length(choices)) choices[1] else character()
    updateSelectInput(session, "tax_sample", choices = choices, selected = selected)
  }

  observeEvent(input$pipeline_step, {
    if (identical(input$pipeline_step, "taxonomy") && nrow(rv$blast_hits)) {
      selected <- if (!is.null(input$tax_sample) && nzchar(input$tax_sample)) input$tax_sample else if (!is.null(input$blast_sample) && nzchar(input$blast_sample)) input$blast_sample else NULL
      update_tax_sample_choices(selected)
    }
  }, ignoreInit = TRUE)

  # Keep the taxonomy selector synchronized as each sequence result arrives,
  # including results retrieved through the batch BLAST workflow.
  observe({
    rv$blast_hits
    if (nrow(rv$blast_hits)) {
      selected <- isolate(input$tax_sample)
      update_tax_sample_choices(selected)
    }
  })


  observeEvent(input$to_taxonomy, {
    if (!nrow(rv$blast_hits)) {
      showNotification("Retrieve at least one BLAST result before taxonomic interpretation.", type="warning")
      return()
    }
    selected <- if (!is.null(input$blast_sample) && nzchar(input$blast_sample)) input$blast_sample else NULL
    update_tax_sample_choices(selected)
    updateTabsetPanel(session, "pipeline_step", selected="taxonomy")
  })

  analyze_taxonomy_sample <- function(sample_key, quiet = FALSE) {
    src <- normalize_blast_hits_unique_accession(blast_hits_for_sample(sample_key))
    if (!nrow(src)) return(list(ok=FALSE, message="No retrieved unique accession-level BLAST hits are available for this sequence."))

    rid <- unique(as.character(src$rid))[1]
    final_name <- unique(as.character(src$final_name))[1]
    top_n <- nrow(src)
    sample_settings <- settings_for_result(rv$results[[sample_key]])
    sample_target <- if (is.list(sample_settings)) stage2_scalar_text(sample_settings$target, "Other") else "Other"

    enriched <- NULL
    result <- NULL
    err <- NULL
    enriched <- tryCatch(
      enrich_hits_with_taxonomy(src),
      error=function(e) { err <<- conditionMessage(e); NULL }
    )
    if (!is.null(enriched)) {
      result <- tryCatch(
        build_taxonomic_consensus(enriched, target=sample_target, top_n=top_n),
        error=function(e) { err <<- conditionMessage(e); NULL }
      )
    }
    if (is.null(result) || !nrow(result$summary)) {
      return(list(ok=FALSE, message=if (!is.null(err)) paste("Taxonomic analysis failed:", err) else "Taxonomic analysis did not produce a summary."))
    }

    tax_error <- attr(enriched, "taxonomy_error")
    sm <- result$summary
    sm$final_name <- final_name
    sm$original_name <- sample_key
    sm$rid <- rid
    sm$target <- sample_target
    sm$max_hits_inspected <- top_n
    sm$analyzed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    sm$app_version <- APP_VERSION
    order_cols <- c("final_name","original_name","rid","target","recommended_identification","recommended_level","confidence")
    sm <- sm[,c(order_cols, setdiff(names(sm), order_cols)),drop=FALSE]

    th <- result$hits
    th$final_name <- final_name
    th$original_name <- sample_key
    th$rid <- rid

    tc <- result$counts
    tc$final_name <- final_name
    tc$original_name <- sample_key
    tc$rid <- rid

    if (nrow(rv$taxonomy_summary) && "rid" %in% names(rv$taxonomy_summary)) rv$taxonomy_summary <- rv$taxonomy_summary[rv$taxonomy_summary$rid != rid,,drop=FALSE]
    if (nrow(rv$taxonomy_hits) && "rid" %in% names(rv$taxonomy_hits)) rv$taxonomy_hits <- rv$taxonomy_hits[rv$taxonomy_hits$rid != rid,,drop=FALSE]
    if (nrow(rv$taxonomy_counts) && "rid" %in% names(rv$taxonomy_counts)) rv$taxonomy_counts <- rv$taxonomy_counts[rv$taxonomy_counts$rid != rid,,drop=FALSE]

    rv$taxonomy_summary <- rbind_fill(rv$taxonomy_summary, sm)
    rv$taxonomy_hits <- rbind_fill(rv$taxonomy_hits, th)
    rv$taxonomy_counts <- rbind_fill(rv$taxonomy_counts, tc)

    msg <- paste0(
      "Analysis complete for ", final_name, ": ",
      sm$recommended_identification[1], " · ", sm$recommended_level[1],
      " · confidence ", sm$confidence[1], "."
    )
    if (!is.null(tax_error) && nzchar(tax_error)) msg <- paste0(msg, " Some NCBI lineage lookups reported: ", tax_error)
    list(ok=TRUE, message=msg, final_name=final_name)
  }

  observeEvent(input$run_taxonomy, {
    req(input$tax_sample)
    rv$taxonomy_status_text <- paste0("Resolving NCBI taxonomy for ", input$tax_sample, "...")
    ans <- NULL
    withProgress(message="Taxonomic interpretation", value=0, {
      incProgress(0.2, detail="Resolving NCBI taxonomy and lineages")
      ans <- analyze_taxonomy_sample(input$tax_sample)
      incProgress(0.8, detail="Preparing consensus")
    })
    rv$taxonomy_status_text <- ans$message
    if (isTRUE(ans$ok)) showNotification("Taxonomic consensus completed.", type="message") else showNotification(ans$message, type="error")
  })

  observeEvent(input$run_taxonomy_all, {
    if (!nrow(rv$blast_hits)) {
      showNotification("Retrieve BLAST results before running batch taxonomy.", type="warning")
      return()
    }
    pairs <- unique(rv$blast_hits[,c("original_name","final_name"),drop=FALSE])
    samples <- as.character(pairs$original_name)
    samples <- samples[nzchar(samples)]
    if (!length(samples)) return()

    ok_n <- 0L; fail_n <- 0L; failures <- character()
    withProgress(message="Analyzing all retrieved sequences", value=0, {
      for (i in seq_along(samples)) {
        incProgress(1/length(samples), detail=paste0(i, "/", length(samples), " · ", pairs$final_name[match(samples[i], pairs$original_name)]))
        ans <- analyze_taxonomy_sample(samples[i], quiet=TRUE)
        if (isTRUE(ans$ok)) ok_n <- ok_n + 1L else { fail_n <- fail_n + 1L; failures <- c(failures, paste0(samples[i], ": ", ans$message)) }
        if (i < length(samples)) Sys.sleep(0.35)
      }
    })
    rv$taxonomy_batch_status_text <- paste0("Batch taxonomy complete: ", ok_n, " analyzed, ", fail_n, " failed.")
    rv$taxonomy_status_text <- rv$taxonomy_batch_status_text
    if (fail_n) showNotification(paste0(rv$taxonomy_batch_status_text, " Review failed samples individually."), type="warning") else showNotification(rv$taxonomy_batch_status_text, type="message")
    update_tax_sample_choices(isolate(input$tax_sample))
  })

  selected_tax_summary <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_summary)) return(rv$taxonomy_summary[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_summary[rv$taxonomy_summary$rid == rid,,drop=FALSE]
  })

  selected_tax_hits <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_hits)) return(rv$taxonomy_hits[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_hits[rv$taxonomy_hits$rid == rid,,drop=FALSE]
  })

  selected_tax_counts <- reactive({
    req(input$tax_sample)
    if (!nrow(rv$taxonomy_counts)) return(rv$taxonomy_counts[0,,drop=FALSE])
    rid <- latest_blast_rid_for_sample(input$tax_sample)
    rv$taxonomy_counts[rv$taxonomy_counts$rid == rid,,drop=FALSE]
  })

  output$taxonomy_status <- renderUI({
    div(class="tax-note", strong("Status: "), rv$taxonomy_status_text)
  })

  output$taxonomy_hero <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)

    txt <- function(name, fallback = "—") {
      if (!name %in% names(sm)) return(fallback)
      x <- as.character(sm[[name]][1])
      if (is.na(x) || !nzchar(trimws(x))) fallback else x
    }
    num <- function(name) {
      if (!name %in% names(sm)) return(NA_real_)
      suppressWarnings(as.numeric(sm[[name]][1]))
    }
    fmt_pct <- function(x, digits=2) if (is.finite(x)) paste0(round(x, digits), "%") else "—"

    recommendation <- txt("recommended_identification", "Unresolved")
    level <- txt("recommended_level", "unresolved")
    confidence <- txt("confidence", "Low / review")
    conf_low <- tolower(confidence)
    conf_key <- if (grepl("high", conf_low)) "high" else if (grepl("moderate", conf_low)) "moderate" else "low"
    hero_class <- paste("taxonomy-result-hero", paste0("conf-", conf_key))
    hero_icon <- if (identical(conf_key, "high")) "check-circle" else if (identical(conf_key, "moderate")) "info-circle" else "exclamation-triangle"

    best_match <- txt("best_molecular_match", txt("species_candidate", ""))
    best_acc <- txt("best_match_accession", "")
    best_id <- num("best_match_identity_percent")
    if (!is.finite(best_id)) best_id <- num("candidate_identity_percent")
    best_cov <- num("best_match_query_coverage_percent")
    if (!is.finite(best_cov)) best_cov <- num("candidate_query_coverage_percent")

    alt <- txt("closest_alternative_species", txt("species_best_competitor", ""))
    alt_id <- num("closest_alternative_identity_percent")
    alt_cov <- num("closest_alternative_query_coverage_percent")
    locus_disc <- txt("locus_discrimination", "Not assessed")
    species_conclusion <- txt("species_level_conclusion", "")
    tax_support <- txt("taxonomic_support")
    seq_evidence <- txt("sequence_evidence")
    ref_support <- txt("reference_support")
    best_n <- num("best_species_database_accessions")
    alt_n <- num("closest_species_database_accessions")

    why_title <- if (tolower(level) == "genus" && nzchar(best_match) && best_match != "—") {
      "Why genus and not species?"
    } else if (tolower(level) == "species") {
      "Why this species?"
    } else if (tolower(level) == "unresolved") {
      "Why unresolved?"
    } else {
      "Why this recommendation?"
    }
    why_text <- txt("decision_reason", "Review the best molecular match and close alternatives below.")

    best_line <- if (nzchar(best_match) && best_match != "—") {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("flask")),
          div(strong("Best molecular match: "), tags$em(best_match),
              if (nzchar(best_acc) && best_acc != "—") paste0(" · ", best_acc) else "",
              paste0(" · Identity ", fmt_pct(best_id), " · Coverage ", fmt_pct(best_cov, 1))))
    } else NULL

    alt_line <- if (nzchar(alt) && alt != "—") {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("exchange-alt")),
          div(strong("Closest alternative: "), tags$em(alt),
              if (is.finite(alt_id) || is.finite(alt_cov)) paste0(" · Identity ", fmt_pct(alt_id), " · Coverage ", fmt_pct(alt_cov, 1)) else ""))
    } else {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("check")),
          div(strong("Closest alternative: "), "No different resolved species inside the close-match window"))
    }

    db_line <- if (is.finite(best_n) && best_n > 0) {
      div(class = "taxonomy-hero-line",
          span(class = "mini-icon", icon("database")),
          div(strong("Database representation: "), paste0(round(best_n), " accession(s) for the best species"),
              if (is.finite(alt_n) && alt_n > 0 && nzchar(alt) && alt != "—") paste0(" · ", round(alt_n), " for the closest alternative") else "",
              span(style="color:#64748b;", " · context only")))
    } else NULL

    div(class = hero_class,
      div(class = "taxonomy-hero-icon-wrap",
        div(class = "taxonomy-hero-icon", icon(hero_icon))
      ),
      div(class = "taxonomy-hero-main",
        div(class = "taxonomy-hero-title",
          recommendation,
          span(class = "level", paste0(" (", level, ")")),
          span(class = "taxonomy-hero-confidence", paste0("— ", confidence))
        ),
        div(class = "taxonomy-hero-lines",
          best_line,
          alt_line,
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("chart-line")),
              div(strong("Evidence: "), "Taxonomic support ", strong(tax_support), " · Sequence evidence ", strong(seq_evidence))),
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("microscope")),
              div(strong("Locus discrimination: "), locus_disc,
                  if (nzchar(species_conclusion) && species_conclusion != "—") paste0(" · ", species_conclusion) else "")),
          div(class = "taxonomy-hero-line",
              span(class = "mini-icon", icon("book")),
              div(strong("Reference support: "), ref_support)),
          db_line
        )
      ),
      div(class = "taxonomy-why",
        div(class = "taxonomy-why-title", span(why_title), info_tip("The recommendation is based on the best molecular match and whether other named taxa are nearly indistinguishable by Identity and query coverage.")),
        div(class = "taxonomy-why-text", why_text)
      )
    )
  })

  output$taxonomy_summary_cards <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)

    sval <- function(name, fallback = "—") {
      if (!name %in% names(sm)) return(fallback)
      x <- as.character(sm[[name]][1])
      if (is.na(x) || !nzchar(trimws(x))) fallback else x
    }
    nval <- function(name) {
      if (!name %in% names(sm)) return(NA_real_)
      suppressWarnings(as.numeric(sm[[name]][1]))
    }
    pct <- function(x, digits=1) if (is.finite(x)) paste0(round(x, digits), "%") else "—"

    ident <- nval("best_match_identity_percent")
    if (!is.finite(ident)) ident <- nval("candidate_identity_percent")
    cov <- nval("best_match_query_coverage_percent")
    if (!is.finite(cov)) cov <- nval("candidate_query_coverage_percent")
    close_n <- nval("close_species_count")
    locus_disc <- sval("locus_discrimination", "Not assessed")
    seq_ev <- sval("sequence_evidence")
    overall <- sval("confidence")
    level <- sval("recommended_level", "unresolved")
    overall_sub <- paste0("Recommendation level: ", level)

    div(class = "taxonomy-metric-grid",
      div(class = "tax-metric metric-purple",
        div(class = "tax-metric-icon", icon("dna")),
        div(class = "tax-metric-label", "Best-match Identity"),
        div(class = "tax-metric-value metric-colored", pct(ident, 2)),
        div(class = "tax-metric-sub", "Best Identity within the near-best coverage band")
      ),
      div(class = "tax-metric metric-blue",
        div(class = "tax-metric-icon", icon("arrows-alt-h")),
        div(class = "tax-metric-label", "Query coverage"),
        div(class = "tax-metric-value metric-colored", pct(cov, 1)),
        div(class = "tax-metric-sub", "Coverage of the selected best molecular match")
      ),
      div(class = "tax-metric metric-green",
        div(class = "tax-metric-icon", icon("clone")),
        div(class = "tax-metric-label", "Close species alternatives"),
        div(class = "tax-metric-value", if (is.finite(close_n)) as.character(round(close_n)) else "—"),
        div(class = "tax-metric-sub", "Nearly indistinguishable named species")
      ),
      div(class = "tax-metric metric-cyan",
        div(class = "tax-metric-icon", icon("microscope")),
        div(class = "tax-metric-label", "Locus discrimination"),
        div(class = "tax-metric-value metric-colored", locus_disc),
        div(class = "tax-metric-sub", paste0("Sequence evidence: ", seq_ev))
      ),
      div(class = "tax-metric metric-amber",
        div(class = "tax-metric-icon", icon("shield")),
        div(class = "tax-metric-label", "Confidence in recommendation"),
        div(class = "tax-metric-value metric-colored", overall),
        div(class = "tax-metric-sub", overall_sub)
      )
    )
  })

  output$taxonomy_summary_table <- renderDT({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(datatable(data.frame(Message="Run multi-hit taxonomic analysis for this sequence."), rownames=FALSE, options=list(dom="t")))

    keep <- intersect(c(
      "algorithm_version","recommended_identification","recommended_level","confidence",
      "best_molecular_match","species_level_conclusion","locus_discrimination"
    ), names(sm))
    df <- sm[, keep, drop=FALSE]
    friendly <- c(
      algorithm_version="Decision engine",
      recommended_identification="Recommended identification",
      recommended_level="Recommended level",
      confidence="Confidence",
      best_molecular_match="Best molecular match",
      species_level_conclusion="Species-level conclusion",
      locus_discrimination="Locus discrimination"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(
      df,
      rownames=FALSE,
      class="compact stripe tax-decision-table",
      options=list(dom="t", autoWidth=FALSE, ordering=FALSE,
                   columnDefs=list(list(targets="_all", className="dt-wrap")))
    )
  })

  output$taxonomy_locus_note <- renderUI({
    sm <- selected_tax_summary()
    if (!nrow(sm)) return(NULL)
    decision <- if ("decision_reason" %in% names(sm) && nzchar(as.character(sm$decision_reason[1]))) as.character(sm$decision_reason[1]) else "Review the best molecular match and close alternatives."
    locus_flag <- if ("locus_flag" %in% names(sm)) as.character(sm$locus_flag[1]) else ""
    tagList(
      div(class="tax-callout",
        div(class="tax-callout-icon", icon("info-circle")),
        div(strong("Decision explanation"), decision)
      ),
      if (nzchar(locus_flag)) div(class="tax-callout",
        div(class="tax-callout-icon", icon("flag")),
        div(strong("Locus flag"), locus_flag)
      ) else NULL
    )
  })

  output$taxonomy_counts_table <- renderDT({
    df <- selected_tax_counts()
    if (!nrow(df)) return(datatable(data.frame(Message="No resolved species evidence profile yet."), rownames=FALSE, options=list(dom="t")))
    keep <- intersect(c(
      "interpretation","taxon","best_identity_percent","best_query_coverage_percent",
      "delta_identity_pp","delta_coverage_pp","accession_count"
    ), names(df))
    df <- df[,keep,drop=FALSE]
    friendly <- c(
      interpretation="Interpretation",
      taxon="Species",
      best_identity_percent="Best Identity (%)",
      best_query_coverage_percent="Best coverage (%)",
      delta_identity_pp="ΔIdentity (pp)",
      delta_coverage_pp="ΔCoverage (pp)",
      accession_count="Accessions"
    )
    names(df) <- unname(friendly[names(df)])
    datatable(
      df,
      rownames=FALSE,
      class="compact stripe species-profile-table",
      options=list(dom="t", autoWidth=FALSE, paging=FALSE, ordering=FALSE,
                   columnDefs=list(
                     list(targets=0, width="115px"),
                     list(targets=1, width="155px", className="dt-taxon"),
                     list(targets=c(2,3,4,5,6), width="74px", className="dt-compact-number")
                   ))
    )
  })

  output$taxonomy_score_plot <- plotly::renderPlotly({
    df <- selected_tax_hits()
    if (!nrow(df)) return(plotly::plot_ly() |> plotly::layout(title="Run taxonomic analysis to display the BLAST score landscape."))

    df$analysis_rank <- if ("analysis_rank" %in% names(df)) suppressWarnings(as.numeric(df$analysis_rank)) else seq_len(nrow(df))
    df$bit_score_num <- suppressWarnings(as.numeric(df$bit_score))
    sp <- clean_taxon_text(df$species)
    gn <- clean_taxon_text(df$genus)
    df$display_taxon <- ifelse(nzchar(sp), sp, ifelse(nzchar(gn), gn, clean_taxon_text(df$organism)))
    resolved_genera <- unique(gn[nzchar(gn)])
    resolved_species <- unique(sp[nzchar(sp)])
    if (length(resolved_genera) > 1) {
      df$color_taxon <- ifelse(nzchar(gn), gn, "Unresolved genus")
      legend_title <- "Genus"
    } else if (length(resolved_species) > 1) {
      df$color_taxon <- ifelse(nzchar(sp), sp, "Unresolved species")
      legend_title <- "Species"
    } else {
      df$color_taxon <- ifelse(nzchar(df$display_taxon), df$display_taxon, "Unresolved")
      legend_title <- "Taxon"
    }

    refq <- if ("reference_quality" %in% names(df)) as.character(df$reference_quality) else ""
    acc <- if ("accession" %in% names(df)) as.character(df$accession) else ""
    ident <- if ("identity_percent" %in% names(df)) suppressWarnings(as.numeric(df$identity_percent)) else NA_real_
    cov <- if ("query_coverage_percent" %in% names(df)) suppressWarnings(as.numeric(df$query_coverage_percent)) else NA_real_
    status <- if ("molecular_evidence_status" %in% names(df)) as.character(df$molecular_evidence_status) else "Database match"
    is_best <- if ("is_best_molecular_match" %in% names(df)) as.logical(df$is_best_molecular_match) else rep(FALSE, nrow(df))
    df$marker_size <- ifelse(is_best, 12, 7)

    df$hover <- paste0(
      "<b>", df$display_taxon, "</b>",
      "<br>Accession: ", acc,
      "<br>Identity: ", round(ident,2), "%",
      "<br>Query coverage: ", round(cov,1), "%",
      "<br>Bit score: ", round(df$bit_score_num,2),
      "<br>Evidence: ", status,
      "<br>Reference: ", refq
    )

    p <- plotly::plot_ly(df, x=~analysis_rank, y=~bit_score_num, type="scatter", mode="lines",
                         line=list(width=1.5, color="#cbd5e1"), hoverinfo="skip", showlegend=FALSE)
    p <- plotly::add_markers(
      p, data=df, x=~analysis_rank, y=~bit_score_num,
      color=~color_taxon, size=~marker_size, sizes=c(7,12),
      text=~hover, hoverinfo="text",
      marker=list(line=list(width=0.5, color="#ffffff")),
      showlegend=TRUE
    )
    p <- plotly::layout(
      p,
      xaxis=list(title="BLAST hit rank", dtick=1, gridcolor="#edf1f6", zeroline=FALSE, titlefont=list(color="#52647a"), tickfont=list(color="#64748b")),
      yaxis=list(title="Bit score", gridcolor="#edf1f6", zeroline=FALSE, titlefont=list(color="#52647a"), tickfont=list(color="#64748b")),
      legend=list(title=list(text=legend_title), orientation="h", x=0, y=1.14, font=list(color="#475569")),
      hovermode="closest",
      paper_bgcolor="rgba(0,0,0,0)",
      plot_bgcolor="rgba(0,0,0,0)",
      margin=list(l=62, r=20, t=35, b=55)
    )
    plotly::config(
      p,
      displaylogo=FALSE,
      modeBarButtonsToRemove=c("lasso2d", "select2d", "hoverCompareCartesian")
    )
  })

  output$taxonomy_hits_table <- renderDT({
    df <- selected_tax_hits()
    if (!nrow(df)) return(datatable(data.frame(Message="No taxonomy-enriched hits yet."), rownames=FALSE, options=list(dom="t")))
    keep <- intersect(c(
      "analysis_rank","organism","genus","species","family","accession","taxid",
      "identity_percent","query_coverage_percent","evalue","bit_score","delta_bit_from_best",
      "reference_quality","hsp_count"
    ), names(df))
    df <- df[,keep,drop=FALSE]
    friendly <- c(
      analysis_rank="Rank",
      organism="NCBI organism",
      genus="Genus",
      species="Species",
      family="Family",
      accession="Accession",
      taxid="TaxID",
      identity_percent="Identity (%)",
      query_coverage_percent="Coverage (%)",
      evalue="E-value",
      bit_score="Bit score",
      delta_bit_from_best="ΔBit",
      reference_quality="Reference",
      hsp_count="HSP"
    )
    names(df) <- unname(friendly[names(df)])

    idx <- function(nm) { w <- which(names(df) == nm); if (length(w)) w[1]-1L else integer() }
    defs <- list(
      list(targets=idx("Rank"), width="42px", className="dt-compact-number"),
      list(targets=idx("NCBI organism"), width="105px", className="dt-organism-compact"),
      list(targets=idx("Genus"), width="85px", className="dt-taxon-compact"),
      list(targets=idx("Species"), width="125px", className="dt-taxon-compact"),
      list(targets=idx("Family"), width="90px", className="dt-wrap-compact"),
      list(targets=idx("Accession"), width="88px", className="dt-nowrap"),
      list(targets=idx("TaxID"), width="62px", className="dt-compact-number"),
      list(targets=unlist(lapply(c("Identity (%)","Coverage (%)","E-value","Bit score","ΔBit"), idx)), width="68px", className="dt-compact-number"),
      list(targets=idx("Reference"), width="92px", className="dt-wrap-compact"),
      list(targets=idx("HSP"), width="40px", className="dt-compact-number")
    )
    defs <- Filter(function(x) length(x$targets) > 0, defs)

    datatable(
      df, rownames=FALSE, filter="top", class="compact stripe taxonomy-hits-compact",
      options=list(
        pageLength=25, scrollX=FALSE, autoWidth=FALSE,
        columnDefs=defs
      )
    )
  })


  all_qc_peak_flags <- reactive({
    if (is.null(rv$results) || !length(rv$results)) return(data.frame())
    collect_ambiguous_peak_flags(rv$results, scope = "trimmed", settings = rv$settings)
  })

  all_curation_log <- reactive({
    if (is.null(rv$results) || !length(rv$results)) return(data.frame())
    collect_curation_log(rv$results)
  })

  team_summary_df <- reactive({
    base <- if (!is.null(rv$summary) && nrow(rv$summary)) export_summary_df() else data.frame()
    if (!nrow(base)) return(data.frame())

    out <- data.frame(
      Sample = as.character(base$final_name),
      Original_sample = as.character(base$sample_id),
      Target = if ("target" %in% names(base)) as.character(base$target) else if (!is.null(rv$settings)) rv$settings$target else "",
      Raw_length = if ("raw_length" %in% names(base)) base$raw_length else NA,
      Trimmed_length = if ("trimmed_length" %in% names(base)) base$trimmed_length else NA,
      QC_status = if ("status" %in% names(base)) as.character(base$status) else "",
      Ambiguous_peak_flags = 0L,
      Strong_peak_flags = 0L,
      Manual_curation_actions = 0L,
      Manual_base_edits = 0L,
      Curation_revision = 0L,
      Identification = "Not analyzed",
      Identification_level = "",
      Overall_confidence = "",
      Taxonomic_support = "",
      Sequence_evidence = "",
      Reference_support = "",
      Best_molecular_match = "",
      Best_match_identity_percent = NA_real_,
      Best_match_query_coverage_percent = NA_real_,
      Closest_alternative_species = "",
      Closest_alternative_identity_percent = NA_real_,
      Closest_alternative_query_coverage_percent = NA_real_,
      Close_species_count = NA_integer_,
      Locus_discrimination = "",
      Species_level_conclusion = "",
      Best_species_database_accessions = NA_integer_,
      Species_candidate = "",
      Species_candidate_confidence = "",
      Candidate_identity_percent = NA_real_,
      Candidate_query_coverage_percent = NA_real_,
      Genus_best_competitor = "",
      Genus_delta_bit = NA_real_,
      Species_best_competitor = "",
      Species_delta_bit = NA_real_,
      BLAST_hits_used = NA_integer_,
      Top_accession = "",
      Top_NCBI_organism = "",
      RID = "",
      Comment = "",
      Application_version = APP_VERSION,
      stringsAsFactors = FALSE
    )

    qcf <- all_qc_peak_flags()
    if (nrow(qcf)) {
      for (i in seq_len(nrow(out))) {
        sf <- qcf[qcf$Sample == out$Original_sample[i], , drop = FALSE]
        if ("Review_status" %in% names(sf)) sf <- sf[sf$Review_status == "Active", , drop = FALSE]
        out$Ambiguous_peak_flags[i] <- nrow(sf)
        out$Strong_peak_flags[i] <- if (nrow(sf)) sum(sf$Severity == "Strong", na.rm = TRUE) else 0L
      }
    }
    clog <- all_curation_log()
    for (i in seq_len(nrow(out))) {
      nm <- out$Original_sample[i]
      r <- rv$results[[nm]]
      if (!is.null(r)) {
        r <- ensure_curation_state(r)
        base_edit_count <- if (is.data.frame(r$curation$base_edits)) nrow(r$curation$base_edits) else 0L
        trim_left_changed <- is.finite(as.numeric(r$curation$trim_start)) && is.finite(as.numeric(r$curation$auto_trim_start)) &&
          as.integer(r$curation$trim_start) != as.integer(r$curation$auto_trim_start)
        trim_right_changed <- is.finite(as.numeric(r$curation$trim_end)) && is.finite(as.numeric(r$curation$auto_trim_end)) &&
          as.integer(r$curation$trim_end) != as.integer(r$curation$auto_trim_end)
        out$Manual_curation_actions[i] <- as.integer(base_edit_count + trim_left_changed + trim_right_changed)
        out$Manual_base_edits[i] <- base_edit_count
        out$Curation_revision[i] <- as.integer(r$curation$revision)
      }
    }

    if (nrow(rv$taxonomy_summary)) {
      for (i in seq_len(nrow(out))) {
        sm <- rv$taxonomy_summary[rv$taxonomy_summary$original_name == out$Original_sample[i],,drop=FALSE]
        if (!nrow(sm)) next
        sm <- sm[nrow(sm),,drop=FALSE]
        out$Identification[i] <- as.character(sm$recommended_identification[1])
        out$Identification_level[i] <- as.character(sm$recommended_level[1])
        out$Overall_confidence[i] <- as.character(sm$confidence[1])
        if ("taxonomic_support" %in% names(sm)) out$Taxonomic_support[i] <- as.character(sm$taxonomic_support[1])
        if ("sequence_evidence" %in% names(sm)) out$Sequence_evidence[i] <- as.character(sm$sequence_evidence[1])
        if ("reference_support" %in% names(sm)) out$Reference_support[i] <- as.character(sm$reference_support[1])
        if ("best_molecular_match" %in% names(sm)) out$Best_molecular_match[i] <- as.character(sm$best_molecular_match[1])
        if ("best_match_identity_percent" %in% names(sm)) out$Best_match_identity_percent[i] <- sm$best_match_identity_percent[1]
        if ("best_match_query_coverage_percent" %in% names(sm)) out$Best_match_query_coverage_percent[i] <- sm$best_match_query_coverage_percent[1]
        if ("closest_alternative_species" %in% names(sm)) out$Closest_alternative_species[i] <- as.character(sm$closest_alternative_species[1])
        if ("closest_alternative_identity_percent" %in% names(sm)) out$Closest_alternative_identity_percent[i] <- sm$closest_alternative_identity_percent[1]
        if ("closest_alternative_query_coverage_percent" %in% names(sm)) out$Closest_alternative_query_coverage_percent[i] <- sm$closest_alternative_query_coverage_percent[1]
        if ("close_species_count" %in% names(sm)) out$Close_species_count[i] <- sm$close_species_count[1]
        if ("locus_discrimination" %in% names(sm)) out$Locus_discrimination[i] <- as.character(sm$locus_discrimination[1])
        if ("species_level_conclusion" %in% names(sm)) out$Species_level_conclusion[i] <- as.character(sm$species_level_conclusion[1])
        if ("best_species_database_accessions" %in% names(sm)) out$Best_species_database_accessions[i] <- sm$best_species_database_accessions[1]
        if ("species_candidate" %in% names(sm)) out$Species_candidate[i] <- as.character(sm$species_candidate[1])
        if ("species_candidate_confidence" %in% names(sm)) out$Species_candidate_confidence[i] <- as.character(sm$species_candidate_confidence[1])
        if ("candidate_identity_percent" %in% names(sm)) out$Candidate_identity_percent[i] <- sm$candidate_identity_percent[1]
        if ("candidate_query_coverage_percent" %in% names(sm)) out$Candidate_query_coverage_percent[i] <- sm$candidate_query_coverage_percent[1]
        if ("genus_best_competitor" %in% names(sm)) out$Genus_best_competitor[i] <- as.character(sm$genus_best_competitor[1])
        if ("genus_delta_bit" %in% names(sm)) out$Genus_delta_bit[i] <- sm$genus_delta_bit[1]
        if ("species_best_competitor" %in% names(sm)) out$Species_best_competitor[i] <- as.character(sm$species_best_competitor[1])
        if ("species_delta_bit" %in% names(sm)) out$Species_delta_bit[i] <- sm$species_delta_bit[1]
        if ("hits_used" %in% names(sm)) out$BLAST_hits_used[i] <- sm$hits_used[1]
        out$RID[i] <- as.character(sm$rid[1])
        if ("decision_reason" %in% names(sm)) out$Comment[i] <- as.character(sm$decision_reason[1])

        bh <- rv$blast_hits[rv$blast_hits$rid == out$RID[i],,drop=FALSE]
        if (nrow(bh)) {
          if ("rank" %in% names(bh)) bh <- bh[order(suppressWarnings(as.numeric(bh$rank))),,drop=FALSE]
          if ("accession" %in% names(bh)) out$Top_accession[i] <- as.character(bh$accession[1])
          if ("organism" %in% names(bh)) out$Top_NCBI_organism[i] <- as.character(bh$organism[1])
        }
      }
    }
    out
  })

  output$team_summary_table <- renderDT({
    df <- team_summary_df()
    if (!nrow(df)) return(datatable(data.frame(Message="No processed sequences yet."), rownames=FALSE, options=list(dom="t")))
    show <- intersect(c(
      "Sample","Target","Trimmed_length","QC_status","Ambiguous_peak_flags","Manual_curation_actions",
      "Identification","Identification_level","Overall_confidence","Best_molecular_match",
      "Best_match_identity_percent","Best_match_query_coverage_percent","Closest_alternative_species",
      "Close_species_count","Locus_discrimination","Reference_support","Comment"
    ), names(df))
    comment_idx <- which(show == "Comment") - 1L
    taxon_idx <- which(show %in% c("Identification","Best_molecular_match","Closest_alternative_species")) - 1L
    defs <- list()
    if (length(comment_idx)) defs[[length(defs)+1]] <- list(targets=comment_idx, width="420px", className="dt-comment")
    if (length(taxon_idx)) defs[[length(defs)+1]] <- list(targets=taxon_idx, width="145px", className="dt-taxon-compact")
    team_view <- df[,show,drop=FALSE]
    if ("Manual_curation_actions" %in% names(team_view)) names(team_view)[names(team_view) == "Manual_curation_actions"] <- "Active_curation_changes"
    datatable(team_view, rownames=FALSE, filter="top", class="compact stripe",
              options=list(pageLength=25, scrollX=TRUE, autoWidth=FALSE, columnDefs=defs))
  })

  output$download_team_summary_csv <- downloadHandler(
    filename=function() paste0(project_export_stem(), "_team_identification_summary.csv"),
    content=function(file) write.csv(team_summary_df(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_team_summary_xlsx <- downloadHandler(
    filename=function() paste0(project_export_stem(), "_team_identification_summary.xlsx"),
    content=function(file) {
      req(nrow(team_summary_df()))
      wb <- openxlsx::createWorkbook(creator="Sanger Sequence Pipeline")
      header_style <- openxlsx::createStyle(fgFill="#1F4E78", fontColour="#FFFFFF", textDecoration="bold", halign="center", valign="center")
      wrap_style <- openxlsx::createStyle(wrapText=TRUE, valign="top")

      add_sheet <- function(name, df, freeze=TRUE) {
        if (is.null(df) || !ncol(df)) df <- data.frame(Message="No data available", stringsAsFactors=FALSE)
        openxlsx::addWorksheet(wb, name)
        openxlsx::writeData(wb, name, df, withFilter=nrow(df)>0, headerStyle=header_style)
        if (freeze) openxlsx::freezePane(wb, name, firstRow=TRUE)
        if (nrow(df)) openxlsx::addStyle(wb, name, wrap_style, rows=2:(nrow(df)+1), cols=seq_len(ncol(df)), gridExpand=TRUE, stack=TRUE)
        openxlsx::setColWidths(wb, name, cols=seq_len(ncol(df)), widths="auto")
        long_cols <- which(names(df) %in% c("Comment","decision_reason","record_title","hit_title","NCBI hit title"))
        if (length(long_cols)) openxlsx::setColWidths(wb, name, cols=long_cols, widths=45)
      }

      add_sheet("Summary", team_summary_df())
      add_sheet("BLAST Hits", rv$blast_hits)
      add_sheet("Taxonomy Details", rv$taxonomy_summary)
      add_sheet("Species Evidence", rv$taxonomy_counts)
      add_sheet("QC", export_summary_df())
      add_sheet("QC Flags", all_qc_peak_flags())
      add_sheet("Manual Curation", all_curation_log())
      if (!is.null(rv$rename)) add_sheet("Rename Map", rv$rename)
      if (is.data.frame(rv$read_assignments) && nrow(rv$read_assignments)) add_sheet("Read Assignments", rv$read_assignments)
      if (is.list(rv$architecture)) {
        add_sheet("Isolates", rv$architecture$isolates)
        add_sheet("Loci", rv$architecture$loci)
        add_sheet("Reads", rv$architecture$reads)
      }
      settings_df <- data.frame(
        Field=c("Application version","Exported at", if (!is.null(rv$settings)) names(rv$settings) else character()),
        Value=c(APP_VERSION, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), if (!is.null(rv$settings)) unlist(rv$settings, use.names=FALSE) else character()),
        stringsAsFactors=FALSE
      )
      add_sheet("Run Settings", settings_df, freeze=FALSE)
      openxlsx::saveWorkbook(wb, file, overwrite=TRUE)
    }
  )

  output$changelog_text <- renderText({
    path <- "CHANGELOG.md"
    if (!file.exists(path)) return("No changelog file found.")
    paste(readLines(path, warn=FALSE), collapse="\n")
  })

  output$download_taxonomy_summary <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_taxonomic_summary.csv"),
    content=function(file) write.csv(selected_tax_summary(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_taxonomy_hits <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_taxonomy_enriched_hits.csv"),
    content=function(file) write.csv(selected_tax_hits(), file, row.names=FALSE, fileEncoding="UTF-8")
  )

  output$download_taxonomy_checkpoint <- downloadHandler(
    filename=function() paste0(clean_fasta_name(ifelse(is.null(input$tax_sample), "sample", input$tax_sample)), "_checkpoint_D_taxonomy.zip"),
    content=function(file) {
      req(nrow(selected_tax_summary()))
      source_hits <- blast_hits_for_sample(input$tax_sample)
      make_taxonomy_checkpoint_zip(file, selected_tax_summary(), selected_tax_hits(), selected_tax_counts(), source_hits)
    }
  )

  # ---------------- Reset ----------------
  reset_pipeline_state <- function() {
    rv$results <- list(); rv$summary <- NULL; rv$rename <- NULL; rv$settings <- NULL
    rv$read_assignments <- stage2_empty_assignments(); rv$architecture <- NULL; rv$project_migration_log <- ""
    rv$blast_jobs <- rv$blast_jobs[0,]; rv$blast_raw <- list(); rv$blast_ids <- data.frame(); rv$blast_hits <- data.frame()
    rv$ncbi_last_contact <- as.POSIXct(NA); rv$blast_batch_status_text <- "No batch operation has been run yet."
    rv$taxonomy_summary <- data.frame(); rv$taxonomy_hits <- data.frame(); rv$taxonomy_counts <- data.frame()
    rv$taxonomy_status_text <- "No taxonomic analysis has been run yet."
    rv$taxonomy_batch_status_text <- "No batch taxonomic analysis has been run yet."
    rv$project_status_text <- "Current session has not been saved as a project."
    rv$project_loaded_name <- ""
    updateTabsetPanel(session,"pipeline_step",selected="upload")
  }
  observeEvent(input$reset_pipeline, reset_pipeline_state())
  observeEvent(input$reset_pipeline_tax, reset_pipeline_state())
}

shinyApp(ui=ui,server=server)
