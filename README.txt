PITAX v3.0.0-alpha.6 — Stage 2 explicit identity build
=====================================================

Purpose
-------
Corrected Stage 2 alpha. It introduces Project -> Isolate -> Locus -> Read as a
backward-compatible metadata architecture while keeping the upload barcode as
the immutable source identifier. Biological identity is entered explicitly and
is never inferred by parsing the uploaded filename.

What changed
------------
- Project schema is now version 3.
- The visible order is Upload -> Assay -> Rename review -> Trim & QC.
- Upload includes manual row editing, expanded batch assignment and XLSX/CSV
  assignment-key import.
- The assignment source of truth is three separate fields: Isolate, Gene/Locus
  and Direction.
- PITAX generates the final read/FASTA name as <Isolate>_<Locus>_<F/R>; the name
  is an output, not a field from which identity is reconstructed.
- Assay opens Rename review immediately; trimming starts only with Start
  trimming.
- A single read is valid. A pair is counted only when the same isolate/locus has
  two distinct source reads assigned Forward and Reverse.
- Paired reads remain separate and auditable; no consensus is created in Stage 2.
- Schema-1 and schema-2 projects migrate in memory without rewriting QC,
  curation, BLAST or taxonomy evidence. Save to persist schema 3.
- Project files and checkpoint ZIPs retain source IDs, explicit assignments,
  generated names and isolate/locus/read architecture tables.

Assignment key
--------------
Use XLSX or CSV columns:

  old_id, isolate, locus, direction

`gene` may be used instead of `locus`. old_id matches the upload barcode/source
name, with unique prefix matching supported.

See `assignment_key_template.csv` for a ready-to-edit example.

Testing
-------
Run `run_tests.bat`.

Expected final line:

  All PITAX v3.0.0-alpha.6 tests passed.

Then follow `V3_STAGE2.md` for the focused manual checks.

Important boundary
------------------
Stage 2 stores identity and pairing metadata only. Reverse-complement alignment,
overlap calculation, conflict resolution and Forward/Reverse consensus belong to
Stage 3 and are intentionally absent from this build.

Branching
---------
Keep this alpha on the v3 development branch. Stable 2.14.2 remains the
production baseline until the gated 3.0 process is complete.

Online deployment
-----------------
If you intentionally deploy this branch separately, regenerate `manifest.json`
locally with:

  source("prepare_connect_cloud.R")

See DEPLOYMENT.md for the stable deployment workflow.
