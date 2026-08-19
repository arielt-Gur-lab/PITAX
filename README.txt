Sanger Sequence Pipeline v2.14.0
================================

Purpose
-------
Local Shiny application for Sanger AB1 processing and sequence identification:

AB1 upload -> Assay & trimming -> QC -> Rename -> Export -> NCBI BLAST -> Taxonomic summary


Key v2.14 identification update
--------------------------------
- Taxonomic interpretation now starts from the best molecular match rather than a primary Bit-score cluster.
- Near-full-length matches are preferred over short partial alignments. Within a comparable coverage tier, Identity is ranked first, then query coverage; Bit score and E-value are tie-breakers.
- A different species is treated as a close alternative when its best comparable match is within 0.5 Identity percentage points and no more than 2 query-coverage points below the best match. These are application review heuristics, not universal species thresholds.
- If close alternatives remain within one genus, the species stays unresolved while the genus can remain a strong recommendation. If close alternatives cross genera, genus confidence is withheld.
- Species evidence profile reports one row per resolved species, including its best accession, Identity, coverage, reference quality and accession count. Counts describe database representation only; they do not vote on the identification.
- Sequence evidence and locus discrimination are separate. A high-quality read can therefore be flagged as having poor species discrimination when several species are nearly indistinguishable with the selected locus.
- The BLAST score landscape colors points by genus when several genera occur and by species when uncertainty is confined to one genus.
- Team and taxonomy tables were compacted; the Team Comment field now receives a wider display column.

Interface and QC retained from v2.13
-----------------------------------
- Aptos is the preferred UI font when available.
- Manual chromatogram curation, editable high-confidence auto-correction criteria, Undo/Redo, audit logging, stale-BLAST protection, batch BLAST retrieval and project Save/Load remain unchanged.
- The v2.14 release changes taxonomic interpretation and presentation; trimming, QC and manual-curation logic are not changed.

Project files
-------------
Use the Project bar near the top of the app:

  Save project  -> downloads a .sangerproject file
  Load project  -> restores a previously saved project

Projects retain processed/curated sequences, chromatogram data, automatic trim state, manual curation history, QC, rename mapping, RIDs, retrieved BLAST hits and taxonomic analyses. v2.14.0 remains backward-compatible with earlier v2 project state; older projects simply begin with an empty curation history.

Run
---
Double-click run_app.bat, or from R:

  shiny::runApp('.', launch.browser = TRUE)

First run
---------
Missing R packages are installed automatically.

Scientific interpretation
-------------------------
BLAST similarity and database agreement support an identification hypothesis; they are not formal taxonomic validation. Species-level resolution depends on locus, clade and reference-database quality. The app exposes each decision component so the user can audit why a genus/species recommendation was made.

Scientific basis highlighted in Help / About
-------------------------------------------
- Schoch et al. (2012): ITS as the universal fungal barcode marker.
- Vu et al. (2018): broad fungal ITS identity benchmarks (~99.6% species, ~94.3% genus).
- NCBI RefSeq Targeted Loci: curated fungal ITS reference resources, largely type-derived.
- NCBI BLAST documentation: Bit scores and HSP/accession-level alignment structure.

Roadmap
-------
v3.0: isolate-level evidence integration, including Forward/Reverse consensus and multi-locus identification.


TESTING
-------
Regression tests can be run from any working directory. Easiest option on Windows: double-click run_tests.bat. It runs taxonomy logic, ambiguous-peak detection, and manual curation/undo-redo smoke tests. These tests are not required for normal app startup.


v2.11 manual curation
---------------------
The QC & Chromatogram stage now separates three layers: immutable raw AB1 evidence, the automatic trim/base calls, and the current curated sequence. Left-click a flag to inspect it and right-click to open curation actions. Sequence-changing actions are confirmed, reversible, logged, exported, and automatically invalidate downstream BLAST/taxonomy results that were generated from an older sequence version.
