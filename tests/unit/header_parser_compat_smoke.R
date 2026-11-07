# Test group 3: Sequence Library fasta_header.py compatibility.
stopifnot(file.exists("VERSION.txt"))
source("R/export/header_templates.R", local = FALSE)

check <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

parser <- "C:/SequenceLibrary/Application/app/parsers/fasta_header.py"
parser_available <- file.exists(parser)
python <- Sys.which("python")
if (!nzchar(python)) python <- Sys.which("python3")

parse_with_python <- function(header) {
  if (!parser_available || !nzchar(python)) return(NULL)
  tmp_py <- tempfile(fileext = ".py")
  tmp_hdr <- tempfile(fileext = ".txt")
  writeLines(header, tmp_hdr, useBytes = TRUE)
  writeLines(c(
    "import json, sys",
    "sys.path.insert(0, r'C:/SequenceLibrary/Application')",
    "from app.parsers.fasta_header import parse_fasta_header",
    "h = open(sys.argv[1], encoding='utf-8').read().strip()",
    "p = parse_fasta_header(h)",
    "print(json.dumps({'name': p.name, 'tags': p.tags, 'extras': p.extras, 'notes': p.notes, 'warnings': p.warnings}))"
  ), tmp_py)
  out <- system2(python, c(tmp_py, tmp_hdr), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    stop(paste("Python parser failed:", paste(out, collapse = "\n")), call. = FALSE)
  }
  jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE)
}

rec <- list(final_name = "FB83_TEF_R", sequence = "ACGT", isolate = "FB83",
            locus = "TEF1", direction = "Reverse", missing = FALSE)

# Build a profile covering every canonical tag.
fields <- c("organism", "isolate", "sample", "locus", "assay", "role", "molecule",
            "accession", "direction", "notes", "pitax_id", "experiment", "seq_date",
            "treatment", "fungicide", "location", "lab_batch")
cfg <- pitax_header_default_config()
cfg$draft$ordered_fields <- setdiff(fields, "lab_batch")
cfg <- pitax_header_create_custom_field(cfg, "Lab batch", tag = "lab_batch", activate = TRUE)$config
cfg$draft$ordered_fields <- fields
cfg$draft$tag_style <- "bare"
cfg <- pitax_header_refresh_derived(cfg)

key <- pitax_header_record_key("qc_trimmed_sequences", "S1")
pairs <- list(
  organism = "Erysiphe necator",
  sample = "leaf-12",
  role = "query",
  accession = "ACC123",
  notes = "manual note",
  experiment = "Fungicide resistance X",
  seq_date = "2026-10-28",
  treatment = "sprayed",
  fungicide = "QoI",
  location = "Golan",
  lab_batch = "B42",
  assay = "TEF1-assay"
)
for (nm in names(pairs)) cfg <- pitax_header_set_value(cfg, key, nm, pairs[[nm]])

header <- pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg)
check(startsWith(header, "FB83_TEF_R "), paste("primary token:", header))
check(grepl("organism=\"Erysiphe necator\"", header, fixed = TRUE) ||
        grepl("organism=Erysiphe", header, fixed = TRUE), "organism present")
check(grepl("isolate=FB83", header, fixed = TRUE), "isolate")
check(grepl("locus=TEF1", header, fixed = TRUE), "locus")
check(grepl("molecule=DNA", header, fixed = TRUE), "molecule default DNA")
check(grepl("direction=R", header, fixed = TRUE), "direction compact")
check(grepl("experiment=", header, fixed = TRUE), "experiment present when selected")
check(grepl("seq_date=2026-10-28", header, fixed = TRUE), "canonical seq_date not date=")
check(grepl("lab_batch=B42", header, fixed = TRUE), "custom extras tag")
check(!grepl("fungicide=$", header) && !grepl("fungicide= ", header), "no empty tags")

if (parser_available && nzchar(python) && requireNamespace("jsonlite", quietly = TRUE)) {
  parsed <- parse_with_python(header)
  check(!is.null(parsed), "python parsed")
  check(identical(parsed$name, "FB83_TEF_R"), "python name")
  check(identical(parsed$tags$organism, "Erysiphe necator"), "python organism")
  check(identical(parsed$tags$isolate, "FB83"), "python isolate")
  check(identical(parsed$tags$sample, "leaf-12"), "python sample")
  check(identical(parsed$tags$locus, "TEF1"), "python locus")
  check(identical(parsed$tags$assay, "TEF1-assay"), "python assay")
  check(identical(parsed$tags$molecule, "DNA"), "python molecule")
  check(identical(parsed$tags$accession, "ACC123"), "python accession")
  check(identical(parsed$tags$direction, "R"), "python direction")
  check(identical(parsed$tags$pitax_id, "S1"), "python pitax_id")
  check(identical(parsed$tags$experiment, "Fungicide resistance X"), "python experiment")
  check(identical(parsed$tags$seq_date, "2026-10-28"), "python seq_date")
  check(identical(parsed$tags$treatment, "sprayed"), "python treatment")
  check(identical(parsed$tags$fungicide, "QoI"), "python fungicide")
  check(identical(parsed$tags$location, "Golan"), "python location")
  check(identical(parsed$extras$lab_batch, "B42"), "python extras")
  # Notes folding includes experiment/date-like fields.
  check(grepl("experiment=", parsed$notes %||% "", fixed = TRUE), "notes fold experiment")
  check(grepl("seq_date=", parsed$notes %||% "", fixed = TRUE), "notes fold seq_date")
  message("PASS: python Sequence Library parser round-trip")
} else {
  message("SKIP: python parser round-trip (parser or python unavailable)")
}

# Bracket style
cfg_b <- cfg
cfg_b$draft$tag_style <- "bracket"
cfg_b <- pitax_header_refresh_derived(cfg_b)
hb <- pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg_b)
check(grepl("[isolate=FB83]", hb, fixed = TRUE), "bracket style")
if (parser_available && nzchar(python) && requireNamespace("jsonlite", quietly = TRUE)) {
  pb <- parse_with_python(hb)
  check(identical(pb$tags$isolate, "FB83"), "python bracket isolate")
}

# Alias acceptance is a parser concern; generator uses canonical tags only.
# Name override
cfg_n <- pitax_header_default_config()
cfg_n$draft$primary_identifier <- "name"
cfg_n$draft$ordered_fields <- c("isolate")
cfg_n <- pitax_header_set_value(cfg_n, key, "name", "OVERRIDE1")
cfg_n <- pitax_header_refresh_derived(cfg_n)
hn <- pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg_n)
check(startsWith(hn, "OVERRIDE1 "), paste("name override:", hn))

# Assay without locus value fails
rec2 <- list(final_name = "X", sequence = "A", isolate = "I", locus = "", direction = "F")
cfg_a <- pitax_header_default_config()
cfg_a$draft$ordered_fields <- c("assay", "locus")
cfg_a <- pitax_header_set_value(cfg_a, pitax_header_record_key("qc_trimmed_sequences", "Sx"),
                               "assay", "A1")
check(inherits(try(pitax_header_render(rec2, "qc_trimmed_sequences", "Sx", cfg_a), silent = TRUE),
               "try-error"),
      "assay without locus value fails")

message("PASS: header_parser_compat_smoke")
