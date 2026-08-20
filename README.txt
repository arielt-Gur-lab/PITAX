PITAX v3.0.0-alpha.8.2 — Stage 3 consensus review
==================================================

Purpose
-------
Stage 2 is closed. Stage 3 creates the auditable sequence used by Export and
BLAST without overwriting either source AB1 read or its curated QC sequence.

Alpha.8–alpha.8.2 changes
------------------------
- Alpha.8.2 fixes Shiny startup/upload navigation by isolating the one-time
  reactive project-mode read and passing mode explicitly to the UI helper.
- Alpha.8.1 fixes application startup by namespace-qualifying the two new
  Plotly consensus chromatogram outputs.
- Simple projects no longer show a redundant Consensus screen. After QC,
  PITAX automatically orients each independent read and continues to Export.
- Paired projects retain Stage 3 and now include Consensus Review & Curation.
- A selected conflict shows focused Forward and Reverse chromatograms at the
  original called-base positions.
- The reviewer may accept the Forward call, Reverse call, or intentionally
  retain the automatic IUPAC ambiguity call.
- Every decision creates a monotonic consensus revision and audit row. Undo and
  Redo are also new audited revisions; they never erase history.
- BLAST jobs retain the exact consensus revision used for each RID.
- Consensus curation audit CSVs are included in checkpoint/result packages.
- Analysis sequence summary widths are fixed and no longer change on sorting.

Project read models
-------------------
- Simple independent reads: one oriented analysis sequence per uploaded read;
  no Forward/Reverse matching and no user-facing consensus step.
- Paired Forward/Reverse: reads pair only from explicit Isolate, Gene/Locus and
  Direction fields. A lone read remains a valid single-read representative.

Scientific contract
-------------------
- One upload/run contains exactly one Gene/Locus. Stage 4 will combine loci.
- Reverse-complementing occurs only in a derived analysis view.
- Weak overlaps are blocked rather than concatenated.
- A mismatch is resolved automatically only with a sufficient quality
  advantage. Otherwise an IUPAC call is retained and Export remains blocked
  until a reviewer records a decision.
- Source sequences, source curation revisions, alignment columns, quality,
  chromatogram positions, automatic calls and manual decisions remain linked.

Controlled test data
--------------------
- `TEST/Stage3_synthetic_pair`: clean trace-aware mirror pair for orientation,
  pairing and overlap mechanics.
- `TEST/Stage3_conflict_pair`: the same controlled material with one
  trace-consistent C/A conflict at Forward-oriented position 400 (IUPAC M).

These fixtures are reproducible software-validation material. They are not an
independently sequenced biological Forward/Reverse pair and do not by themselves
close biological acceptance of Stage 3.

Testing
-------
Run `run_tests.bat` on Windows with R installed.

Expected final line:

  All PITAX v3.0.0-alpha.8.2 tests passed.

Roadmap
-------
- Stage 3: complete biological acceptance with a real paired-AB1 truth set.
- Stage 4: multi-locus isolate profile with separate evidence per locus.
- PITAX 4 preparation: stored alignment, ITS region annotation, then
  phylogenetic trees based on the reviewed alignment.
