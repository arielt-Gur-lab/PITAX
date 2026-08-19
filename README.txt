PITAX v3.0.0-alpha.4 — Stage 1 gate-closing build
================================

Purpose
-------
Development branch build for the first PITAX 3.0 gate: validating the AB1 evidence model without changing the established v2.14.2 trimming, curation, BLAST or taxonomy decisions.

User-facing workflow:

AB1 upload -> Assay & trimming -> Rename -> QC -> Export -> NCBI BLAST -> Taxonomic summary

Stage 1 adds an observational AB1 evidence audit inside QC.

alpha.4 changes
----------------
- Moves Rename before QC so chromatograms and QC tables use resolved sample names.
- Trimming now continues directly to Rename; validated names continue to QC; QC continues to Export.
- QC selectors, tables, chromatogram titles and audit displays show the resolved name while preserving the original sample ID in audit exports.
- Reorders checkpoints: A = renamed sequences; B = renamed and curated QC state.
- Cleans the chromatogram navigator by clipping base letters and QC markers from its compressed overview.
- Keeps the validated alpha.3 AB1 evidence model and same-length PCON comparison unchanged and observational only.
- Adds an automated workflow-order regression test.

Testing
-------
Run `run_tests.bat`.

Expected final line:

  All PITAX v3.0.0-alpha.4 tests passed.

Then follow `V3_STAGE1.md`. The final manual regression is the Rename -> QC transition and resolved-name display. Once that is green, Stage 1 is closed.

Branching
---------
Keep this alpha on the v3 development branch. Do not replace the stable 2.14.2 `main` deployment while Stage 1 is under validation.

Online deployment
-----------------
The stable app can remain on `main` in Posit Connect Cloud. If you intentionally deploy this branch separately, regenerate `manifest.json` locally with:

  source("prepare_connect_cloud.R")

See DEPLOYMENT.md for the stable deployment workflow.
