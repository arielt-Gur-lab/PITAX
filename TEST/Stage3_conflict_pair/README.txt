PITAX Stage 3 controlled conflict pair
======================================

Upload TESTCONFLICT001_forward.ab1 and TESTCONFLICT001_reverse.ab1 in the
Paired Forward/Reverse project model, then import assignment_key.csv.

This fixture is derived from the same permissively licensed public test trace
documented in ../Stage3_synthetic_pair/README.txt. The Reverse file is a
trace-aware reverse complement with one controlled, trace-consistent call
conflict at Forward-oriented called-base position 400:

  Forward call: C
  Reverse-oriented call: A
  Automatic ambiguity call: M

The pair is intended to exercise Stage 3 Review & Curation, including linked
chromatograms, acceptance of Forward / Reverse / IUPAC, revision history and
Undo / Redo. It is deterministic software-validation material, not an
independently sequenced biological pair and therefore does not replace a real
paired-AB1 acceptance set.
