PITAX v3.0.0-alpha.3 — Stage 1
================================

Purpose
-------
Development branch build for the first PITAX 3.0 gate: validating the AB1 evidence model without changing the established v2.14.2 trimming, curation, BLAST or taxonomy decisions.

Normal workflow remains:

AB1 upload -> Assay & trimming -> QC -> Rename -> Export -> NCBI BLAST -> Taxonomic summary

Stage 1 adds an observational AB1 evidence audit inside QC.

alpha.3 changes
----------------
- Keeps the corrected alpha.2 ABIF evidence model and export synchronization.
- Adds explicit auto-trim start, end and length to the Stage 1 run audit.
- Adds `In_auto_trim` to every exported per-base audit row.
- Adds a same-length PCON-only quality-window proposal for comparison. It is ranked by Q>=20, then Q>=30, median quality, mean quality and coverage.
- Shows legacy auto trim and the PCON-only proposal side by side for the selected sample.
- Adds `In_quality_proposed_window` to the per-base table and CSV.
- The proposed quality window is observational only: it does not change trimming, FASTA, curation, BLAST or taxonomy.
- Simplifies the Plotly title to the sample ID, avoiding the Windows byte-marker rendering of the em dash.
- Moves the chromatogram legend outside the data area and separates called-base labels from QC markers.

Testing
-------
Run `run_tests.bat`.

Expected final line:

  All PITAX v3.0.0-alpha.3 tests passed.

Then follow `V3_STAGE1.md`. The key manual regression is to process `265DMAA001_customer`, confirm its active auto trim remains unchanged, and compare it with the observational PCON-only window.

Branching
---------
Keep this alpha on the v3 development branch. Do not replace the stable 2.14.2 `main` deployment while Stage 1 is under validation.

Online deployment
-----------------
The stable app can remain on `main` in Posit Connect Cloud. If you intentionally deploy this branch separately, regenerate `manifest.json` locally with:

  source("prepare_connect_cloud.R")

See DEPLOYMENT.md for the stable deployment workflow.
