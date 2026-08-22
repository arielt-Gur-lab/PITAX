PITAX v3.0.0-alpha.10.3 — Alpha 10 schema foundation
=========================================================

Purpose
-------
Stage 3 is accepted for continued development using the controlled clean and
single-conflict AB1 fixtures. A real independently sequenced Forward/Reverse
pair remains a deferred biological validation item.

Alpha 10.3 retains the visual review and multi-locus layers, adds the schema-6
assay foundation, and makes runtime source encoding independent of the Windows
or Connect locale. Raw non-ASCII UI symbols are represented by Unicode escapes
and protected by a source contract.

The authoritative guide for continuing development is:

  docs/PITAX_MASTER.md

Stage 4 workflow
----------------
1. Process and save each Gene/Locus as its own PITAX project.
2. Open `9 · Multi-locus`.
3. Include the current session and/or add completed `.sangerproject` files.
4. Build the profile. Projects are joined only by the explicit Isolate field.
5. Select each isolate in the visual workspace and compare its locus cards,
   combined status and Identity/query-coverage plot.
6. Review the retained profile and per-locus audit tables.
7. Export CSV, FASTA or Checkpoint F ZIP.

Scientific contract
-------------------
- One source project must contain exactly one Gene/Locus.
- Every Isolate/Locus combination must be unique across imported projects.
- Every locus retains its sequence, consensus revision, BLAST RID, taxonomic
  interpretation, reference context and source-project fingerprint.
- Concordant loci may support the same rank.
- A genus or species conflict is retained and no combined call is reported.
- Missing locus-level taxonomy is shown as partial evidence.
- Taxon-specific marker recommendations are not automated in Alpha 10.3; that
  layer will be added only after literature review and curated scenario tests.

Stage 3 retained behavior
-------------------------
- Simple projects orient each read and bypass the Consensus screen.
- Paired projects create an auditable Forward/Reverse consensus.
- Weak overlaps are blocked rather than concatenated.
- Unresolved mismatches retain an International Union of Pure and Applied
  Chemistry (IUPAC) ambiguity call until reviewed.
- Manual consensus decisions, Undo/Redo and BLAST queries remain revision-bound.

Controlled Stage 3 test data
----------------------------
- `tests/fixtures/Stage3_synthetic_pair`: clean trace-aware mirror pair.
- `tests/fixtures/Stage3_conflict_pair`: one trace-consistent C/A conflict at
  Forward-oriented position 400 (IUPAC M).

Testing
-------
Run `run_tests.bat` on Windows with R installed.

Expected final line:

  All PITAX v3.0.0-alpha.10.3 tests passed.

The suite now has 15 groups. Alpha 10 tests cover controlled loci, assay-linked reads,
schema-5 migration and multi-assay architecture. Stage 4 tests prove duplicate Isolate/Locus
blocking, provenance retention, missing-evidence reporting, correct separation
of multiple isolates in the visual selector and that a 2:1 locus majority
cannot vote away a cross-genus conflict. Group 13 guards horizontal
DataTables header/body alignment throughout the application. Group 14 protects
the organized file/module boundaries and rejects raw non-ASCII runtime source.

Release update workflow
-----------------------
- Keep the existing Git working folder and its hidden `.git` directory.
- Copy release files over that folder; release ZIP files do not contain `.git`.
- Run `run_tests.bat`, then `Git.BAT`.
- `Git.BAT` fetches `origin/main` and stops before commit if the remote branch
  contains changes that are not present locally.

Roadmap
-------
- Alpha 10 next: schema 6, controlled locus vocabulary and multiple Assay
  profiles inside one project, with direct multi-locus analysis by Isolate.
- Stage 4 later: literature-backed marker profiles and curated fungal scenarios.
- PITAX 4 preparation: stored alignment, ITS region annotation, then
  phylogenetic trees based on the reviewed alignment.
