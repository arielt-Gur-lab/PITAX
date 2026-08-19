PITAX v3.0.0-alpha.2 — Stage 1
================================

Purpose
-------
Development branch build for the first PITAX 3.0 gate: validating the AB1 evidence model without changing the established v2.14.2 trimming, curation, BLAST or taxonomy decisions.

Normal workflow remains:

AB1 upload -> Assay & trimming -> QC -> Rename -> Export -> NCBI BLAST -> Taxonomic summary

Stage 1 adds an observational AB1 evidence audit inside QC.

alpha.2 changes
----------------
- Corrected the alpha.1 interpretation of raw sangerseqR peak matrices.
- Raw ABIF primary positions are now read from PLOC.2 (+1 to match R/sangerseqR coordinates).
- PCON.2/PCON.1 basecaller quality is recorded when available.
- traceMatrix is treated as canonical A,C,G,T evidence.
- Existing PITAX inferred mapping and primary positions are compared against those raw/canonical sources without changing trimming.
- Added auto-trim quality summaries (median quality, Q>=20, Q>=30, called-base dominance and called/alternative signal ratio).
- Fixed selected-base audit export: selection, result content and filename now use one sample identity; every exported row includes Sample_ID; mismatches are blocked.
- Clicking a run-audit row selects that same sample in the QC workspace.

Testing
-------
Run `run_tests.bat`.

Expected final line:

  All PITAX v3.0.0-alpha.2 tests passed.

Then follow `V3_STAGE1.md`. The key manual regression is to select `265DMAA001_customer`, download the selected-base audit CSV and confirm both filename and every Sample_ID row identify 001.

Branching
---------
Keep this alpha on the v3 development branch. Do not replace the stable 2.14.2 `main` deployment while Stage 1 is under validation.

Online deployment
-----------------
The stable app can remain on `main` in Posit Connect Cloud. If you intentionally deploy this branch separately, regenerate `manifest.json` locally with:

  source("prepare_connect_cloud.R")

See DEPLOYMENT.md for the stable deployment workflow.
