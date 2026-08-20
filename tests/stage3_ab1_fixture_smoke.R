# PITAX v3.0.0-alpha.8 - controlled clean and conflict AB1 pair validation.

source(file.path("R", "stage3_consensus.R"))

fixture_dir <- file.path("TEST", "Stage3_synthetic_pair")
forward_path <- file.path(fixture_dir, "TESTPAIR001_forward.ab1")
reverse_path <- file.path(fixture_dir, "TESTPAIR001_reverse.ab1")
stopifnot(file.exists(forward_path), file.exists(reverse_path))

forward_abif <- sangerseqR::read.abif(forward_path)
reverse_abif <- sangerseqR::read.abif(reverse_path)
forward_read <- sangerseqR::readsangerseq(forward_path)
reverse_read <- sangerseqR::readsangerseq(reverse_path)

forward_calls <- as.character(sangerseqR::primarySeq(forward_read))
reverse_calls <- as.character(sangerseqR::primarySeq(reverse_read))
stopifnot(identical(stage3_reverse_complement(reverse_calls), forward_calls))
stopifnot(identical(rev(as.integer(forward_abif@data$PCON.2)), as.integer(reverse_abif@data$PCON.2)))

trace_length <- length(forward_abif@data$DATA.9)
expected_reverse_positions <- trace_length - 1L - rev(as.integer(forward_abif@data$PLOC.2))
stopifnot(identical(expected_reverse_positions, as.integer(reverse_abif@data$PLOC.2)))
stopifnot(all(diff(as.integer(reverse_abif@data$PLOC.2)) >= 0L))

# FWO_.1 is GATC in this fixture: DATA.9=G, DATA.10=A, DATA.11=T,
# DATA.12=C. Reverse-complementing therefore swaps A/T and C/G while reversing
# the trace-sample order.
stopifnot(identical(rev(as.integer(forward_abif@data$DATA.11)), as.integer(reverse_abif@data$DATA.10)))
stopifnot(identical(rev(as.integer(forward_abif@data$DATA.10)), as.integer(reverse_abif@data$DATA.11)))
stopifnot(identical(rev(as.integer(forward_abif@data$DATA.9)), as.integer(reverse_abif@data$DATA.12)))
stopifnot(identical(rev(as.integer(forward_abif@data$DATA.12)), as.integer(reverse_abif@data$DATA.9)))

conflict_dir <- file.path("TEST", "Stage3_conflict_pair")
conflict_forward <- as.character(sangerseqR::primarySeq(sangerseqR::readsangerseq(file.path(conflict_dir, "TESTCONFLICT001_forward.ab1"))))
conflict_reverse <- as.character(sangerseqR::primarySeq(sangerseqR::readsangerseq(file.path(conflict_dir, "TESTCONFLICT001_reverse.ab1"))))
conflict_oriented <- stage3_reverse_complement(conflict_reverse)
mismatch_positions <- which(strsplit(conflict_forward, "", fixed = TRUE)[[1]] != strsplit(conflict_oriented, "", fixed = TRUE)[[1]])
stopifnot(identical(mismatch_positions, 400L))
stopifnot(substr(conflict_forward, 400, 400) == "C")
stopifnot(substr(conflict_oriented, 400, 400) == "A")
stopifnot(stage3_iupac_pair("C", "A") == "M")

cat("v3.0.0-alpha.8 controlled AB1 pair tests passed.\n")
