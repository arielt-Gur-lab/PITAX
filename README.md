PITAX v3.6.5 - Browser-safe custom FASTA ZIP download
=========================================================

Purpose
-------
Stage 3 is accepted for continued development using the controlled clean and
single-conflict AB1 fixtures. A real independently sequenced Forward/Reverse
pair remains a deferred biological validation item.

Alpha 10.4 keeps the visual review and multi-locus layers and replaces the
scattered result-download handlers with one shared Export dialog. Packages are
built from an explicit component catalog, a frozen project snapshot and a
RUN_INFO.txt manifest. Save project and the blank assignment-key template remain
separate from result export.

The authoritative guide for continuing development is:

  docs/PITAX_MASTER.md

Custom FASTA Header Builder (Export)
------------------------------------
1. Open Export from any stage after read assignments/sequence processing.
2. Under FASTA headers, select "Custom FASTA headers".
3. Edit the one-line token template, e.g.:
   {Final_Name} organism={Organism} experiment="{Experiment}" date={Seq_Date}
4. Double-click metadata cells for individual records, or select rows and
   use bulk field/value updates; use All records for global experiment values.
5. Add custom fields (e.g. Fungicide or Location) and use {Fungicide}
   in the template. Built-ins are Final_Name, Isolate, Locus, Direction,
   Source_ID and Record_ID.
6. Review the live preview of the FASTA components actually selected.
   Only fields referenced in the saved active template are required. Removing
   experiment="{Experiment}" from the template removes the experiment field
   from the header and its validation requirements, even if the metadata
   editor no longer shows the unused Experiment column for that template.\n   Other named templates may continue to use Experiment.
   Resolve missing values/duplicate identifiers before downloading.
   On error, the app now explains the failing record/field in the export
   dialog instead of navigating to an HTTP 500 page.
7. After pressing **1. Prepare ZIP**, wait until the green **2. Download ready ZIP**
   button appears and click it. This is an actual Shiny browser download link;
   the old invisible automatic download click is no longer used.
7. Download the ZIP: per-read filenames remain 315_TEF_R.fasta while headers
   contain the chosen metadata. FASTA_header_manifest.csv lists the outputs.

Saved templates and metadata values persist in the .sangerproject and shared
snapshots. Template CSV export/import transfers reusable templates and their
field names between projects; per-record values do not transfer automatically.
Sequencing dates are entered by the user, not guessed from the export date.
These headers never change the authoritative read/analysis identity or BLAST.

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
- A project may contain multiple Gene/Locus assays; every Isolate/Locus combination
  must remain unique across evidence sources.
- Every locus retains its sequence, consensus revision, BLAST RID, taxonomic
  interpretation, reference context and source-project fingerprint.
- Concordant loci may support the same rank.
- A genus or species conflict is retained and no combined call is reported.
- Missing locus-level taxonomy is shown as partial evidence.
- Taxon-specific marker recommendations are not automated in Alpha 10; that
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

  PITAX 3.6.5 includes a native Shiny download integration test.

The suite now has 22 groups. Alpha 10 covers controlled loci, multi-assay profiles,
schema-5 migration, gated Steps UI navigation and same-project multi-locus integration.
Stage 4 tests prove duplicate Isolate/Locus blocking, provenance retention, missing-evidence
reporting, correct separation of multiple isolates in the visual selector and that a 2:1
locus majority cannot vote away a cross-genus conflict. Group 13 guards horizontal
DataTables header/body alignment throughout the application. Group 14 protects
the organized file/module boundaries and rejects raw non-ASCII runtime source.
Group 16 covers reflection follow-up correctness (invalidation breadth, Assay_ID
signature, READY-only taxonomy RID binding and BLAST revision checks).
Group 20 covers the shared export engine (catalog availability, package planning,
sequence roles, BLAST partial states and RUN_INFO.txt contents).

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
