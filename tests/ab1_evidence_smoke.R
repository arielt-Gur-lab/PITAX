# PITAX v3.0.0-alpha.2 - Stage 1 AB1 evidence helper tests.
# Synthetic tests for corrected raw-ABIF interpretation and export identity.

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "ab1_evidence.R"))

if (!methods::isClass("PitaxTestAbif")) methods::setClass("PitaxTestAbif", slots = c(data = "list"))
abif <- methods::new(
  "PitaxTestAbif",
  data = list(
    "PLOC.2" = c(9, 19, 29, 39, 49),
    "PCON.2" = c(5, 25, 40, 35, 10),
    "FWO.1" = "GATC"
  )
)

# 1. Raw ABIF PLOC.2 is converted to the same 1-based coordinate system used
# by sangerseqR's primary peak positions.
pos <- pitax_abif_primary_positions(abif, 5)
stopifnot(identical(pos$tag, "ABIF PLOC.2 + 1"))
stopifnot(identical(as.numeric(pos$values), c(10,20,30,40,50)))

# 2. PCON.2 is retained by base index without shifting.
q <- pitax_abif_quality(abif, 5)
stopifnot(identical(q$tag, "PCON.2"))
stopifnot(identical(as.numeric(q$values), c(5,25,40,35,10)))

# 3. sangerseqR traceMatrix is treated as canonical A/C/G/T.
trace <- matrix(1, nrow = 60, ncol = 4)
colnames(trace) <- c("A","C","G","T")
map <- pitax_canonical_channel_map(trace)
stopifnot(identical(unname(as.integer(map[c("A","C","G","T")])), 1:4))
stopifnot(pitax_maps_identical(map, setNames(1:4, c("A","C","G","T"))))

# 4. Evidence is evaluated at the primary ABIF base-call positions, not by
# treating peakPosMatrix columns as A/C/G/T.
seq_string <- "ACGTN"
trace[10,] <- c(100,10,8,5)
trace[20,] <- c(8,110,9,7)
trace[30,] <- c(6,8,120,9)
trace[40,] <- c(7,6,8,130)
legacy_pos <- c(10,20,30,40,50)
legacy_map <- setNames(1:4, c("A","C","G","T"))
ev <- build_ab1_evidence(NULL, seq_string, trace, legacy_pos, legacy_map, abif)
stopifnot(identical(ev$schema, "ab1-evidence-audit-v2"))
stopifnot(all(ev$detail$Primary_peak_pos_delta == 0))
stopifnot(all(ev$detail$Canonical_called_is_max[1:4] %in% TRUE))
stopifnot(identical(ev$detail$Canonical_best_channel[1:4], c("A","C","G","T")))
stopifnot(ev$summary$channel_maps_match[1] %in% TRUE)
stopifnot(ev$summary$different_primary_position_percent[1] == 0)

# 5. Run summary computes quality metrics specifically inside auto trim.
res <- list(
  sample_id = "S001",
  ab1_evidence = ev,
  curation = list(auto_trim_start = 2L, auto_trim_end = 4L),
  summary = data.frame(trim_start = 2L, trim_end = 4L)
)
sm <- ab1_evidence_result_summary(res)
stopifnot(sm$Sample[1] == "S001")
stopifnot(sm$Median_quality_auto_trim[1] == 35)
stopifnot(sm$Q20_auto_trim_percent[1] == 100)
stopifnot(sm$Q30_auto_trim_percent[1] == round(100 * 2/3, 2))
stopifnot(sm$Canonical_called_is_max_auto_trim_percent[1] == 100)

# 6. Selected-base export contains its Sample_ID and refuses a selection/result
# mismatch instead of silently writing a misleading filename/content pair.
d <- pitax_evidence_detail_export(res, "S001")
stopifnot("Sample_ID" %in% names(d))
stopifnot(all(d$Sample_ID == "S001"))
blocked <- tryCatch({ pitax_evidence_detail_export(res, "S999"); FALSE }, error = function(e) TRUE)
stopifnot(blocked)

# 7. Run-table row selection resolves back to the actual result key.
results <- list(S001 = res, S002 = modifyList(res, list(sample_id = "S002")))
stopifnot(pitax_result_key_for_sample(results, "S002") == "S002")

cat("v3.0.0-alpha.2 AB1 evidence helper tests passed.\n")
