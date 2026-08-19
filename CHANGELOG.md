# Changelog

This file records changes that affect the way the pipeline behaves, presents results, or stores project data. Small internal refactors are omitted unless they change something visible or scientifically relevant.

## 3.0.0-alpha.2
**Stage 1 evidence correction and export synchronization**
- Corrected the Stage 1 ABIF model: fresh `sangerseqR` `peakPosMatrix[,1]` is treated as the primary ABIF base-call coordinate rather than an A-channel peak column.
- Raw primary positions are audited against `PLOC.2 + 1`; the invalid alpha.1 A/C/G/T interpretation of raw `peakPosMatrix`/`peakAmpMatrix` was removed.
- Added auto-trim PCON summaries (median, Q>=20, Q>=30) and canonical trace competition metrics.
- Fixed selected-base audit downloads so sample selection, filename and file contents share one identity; CSV rows now include `Sample_ID` and mismatches are blocked.
- Run-audit row selection now updates the QC sample selector.
- No change to the v2.14.2 trimming, manual curation, BLAST or taxonomy decision paths.

## 3.0.0-alpha.1
**Stage 1 · AB1 evidence foundation**
- Starts the PITAX 3.0 migration without changing the established v2.14.2 trimming, curation, BLAST or taxonomy decision paths.
- Captures ABIF basecaller quality values from `PCON.2`/`PCON.1` when present.
- Captures the full `peakPosMatrix` and `peakAmpMatrix` supplied by `sangerseqR`.
- Adds an observational A/C/G/T evidence model that selects the peak-position column corresponding to the called base and compares it with the v2 legacy first-column position model.
- Adds a collapsed **Stage 1 · AB1 evidence audit** panel in QC with run-level and per-base diagnostics plus CSV downloads.
- Evidence-audit failures are isolated from the production read path: an audit error is reported but does not invalidate a trimming result that the v2.14.2 engine can process.
- Adds a synthetic regression test for documented A/C/G/T channel mapping, called-base peak positions, signal evidence, `peakAmpMatrix` handling and quality-vector alignment.
- Adds the five-stage v3 roadmap and a manual real-AB1 validation checklist.

## 2.14.2
**PITAX branding and deployment preparation**
- Added the PITAX logo to the application header. The source image remains `logo.png` in the project root and is served without exposing the source directory.
- Changed the browser title and main application branding to PITAX.
- Added `prepare_connect_cloud.R` and `DEPLOYMENT.md` for GitHub-based deployment to Posit Connect Cloud.
- No trimming, QC, BLAST or taxonomic decision logic changed in this release.

## 2.14.1
**Loading feedback and match-selection cleanup**
- Added a visible loading overlay when changing workflow steps, with a thin activity bar while Shiny is processing. The QC workspace now shows a loading state after AB1 processing instead of appearing temporarily blank.
- Team summary now reports active sequence curation changes rather than counting every audit-log event. Marking a flag as reviewed without changing the sequence no longer appears as a manual sequence edit.
- Tightened best-molecular-match selection: hits are first restricted to a near-best coverage band (within 2 percentage points of the highest query coverage), then ranked by Identity. This prevents a noticeably shorter hit from being highlighted only because its Identity is marginally higher.

## 2.14.0
**Evidence-first taxonomic interpretation**
- Replaced the primary-score-cluster decision rule with a simpler evidence-first tree: best molecular match first, then close alternatives, then the most conservative supported taxonomic level.
- Best match is selected from comparable-coverage hits by Identity, then query coverage; short partial matches cannot outrank near-full-length hits just because their Identity is 100%.
- Added a Species evidence profile with the best accession, Identity, coverage and accession count for each resolved species. Accession count is shown as database context only and is not used as a majority vote.
- Added locus-discrimination reporting. High-quality reads that match several species almost equally are flagged as poor species-level discrimination rather than poor sequence quality.
- Updated the taxonomy dashboard, Help / About, score-landscape coloring and exports for the new logic.
- Taxonomy-enriched hits no longer show Primary score cluster; the table is compacted to avoid horizontal scrolling on normal desktop widths.
- Team summary gives the Comment column substantially more width.

## 2.13.4
**BLAST workspace fixes**
- The selected processed sequence is shown reliably in Query workspace, including after returning from Taxonomic summary or loading an older project state.
- Submit selected now checks for an existing job with the same sequence, database and requested hit count before contacting NCBI. Matching SUBMITTED, WAITING or READY jobs are reused instead of creating another RID.
- The BLAST jobs table now supports selecting multiple rows. Retrieve selected checks all selected jobs in one operation; if no rows are selected it falls back to the newest job for the sequence shown in Query workspace.

## 2.13.2
**Version history display**
- Rewrote the changelog to be shorter and easier to scan.
- Version history now wraps long lines inside the page instead of requiring horizontal scrolling.

## 2.13.1
**BLAST resubmission**
- An existing BLAST job is reused only when both the database and requested maximum hit count match the current settings.
- Changing either setting creates a new BLAST request and a new RID.
- Matching active jobs are still skipped so the same request is not submitted twice by mistake.
- The BLAST jobs table now keeps the database used for each RID. Older saved projects remain compatible.

## 2.13.0
**UI pass across the full workflow**
- Applied the same visual language to Upload, Assay & Trim, QC & Chromatogram, Rename, Export, NCBI BLAST and Taxonomic summary.
- Switched the main interface font to Aptos when available, with system fallbacks.
- Reorganized each stage around a clear heading, compact controls and wider working areas for tables and plots.
- QC keeps the chromatogram and manual-curation tools prominent; Rename and BLAST use clearer two-column workspaces; Export is grouped by output type.
- No scientific or processing logic changed in this release.

## 2.12.1
**Startup fix**
- Fixed a quoting error in the embedded CSS that prevented v2.12.0 from sourcing correctly.
- No change to the taxonomy engine or processing workflow.

## 2.12.0
**Taxonomy dashboard redesign**
- Rebuilt the Taxonomic summary around a single recommendation panel and compact evidence cards.
- Genus agreement, species agreement, competitor ΔBit, sequence evidence and overall confidence are shown separately.
- Reduced the amount of explanatory text in the working view and moved detail into tooltips and Help / About.
- Simplified the BLAST score-landscape controls and added a short expandable guide for reading the graph.
- The underlying taxonomy decision logic was unchanged.

## 2.11.1
**Consensus and close competitors**
- Fixed a case where one near-tied competitor could erase an otherwise strong primary-cluster consensus.
- A close competitor now lowers confidence instead of automatically forcing an Unresolved result.
- Species promotion still requires strong sequence evidence and adequate separation; otherwise the supported genus is reported and the species remains a candidate for review.
- Added a regression case for a 24/25 Pleurotus cluster with a close Lentinus hit.

## 2.11.0
**Manual sequence curation**
- Added a curated-sequence layer without modifying the original AB1 trace or automatic base calls.
- Right-clicking a flagged chromatogram position allows manual base replacement, a two-base IUPAC ambiguity call, directional trimming, or marking the call as reviewed.
- Sequence-changing actions require confirmation and are written to a Manual Curation audit log.
- Added Undo, Redo, History and Reset to auto for each sample.
- Added a conservative bulk correction tool. It previews only high-confidence proposals and requires confirmation before applying them.
- If a curated sequence changes after BLAST, the old BLAST/taxonomy result is marked stale and the sequence must be submitted again.
- Manual-curation history is included in project exports.

## 2.10.4
**Ambiguous peak review**
- Added flags for positions where another dye channel strongly competes with the called base, the called base is not locally dominant, or the call is N.
- Clicking a flag centers the interactive chromatogram on that position.
- Added optional flag markers on the chromatogram and QC flag counts in exports.
- The default review scope is the retained sequence; the full raw read can also be inspected.

## 2.10.3
**Taxonomy and Help fixes**
- Fixed a Help / About tab configuration that could stop the Shiny app from loading.
- Restored the intended rule that Moderate species evidence keeps the species as a candidate rather than automatically promoting it.
- Clarified the message shown when no resolved alternative species or genus is found.

## 2.10.2
**Genus fallback**
- Fixed a case where an ITS identity benchmark was acting like a hard cutoff and could erase a strong genus consensus.
- Identity and coverage now influence evidence strength and species promotion without automatically removing a supported genus call.
- Removed flat Top-N voting from taxonomic interpretation. Analysis uses the unique accession-level hits retrieved from BLAST.
- Removed the redundant Decision explanation row from the interpretation table and simplified the score-landscape toolbar.

## 2.10.1
**Help / About**
- Split Help / About into separate sections for workflow, trimming/QC, BLAST, taxonomy, references and version history.
- Added clickable scientific and NCBI references.
- Clarified which values come from published literature and which rules are application heuristics.
- Fixed the taxonomy smoke test so it can be run from outside the application working directory.

## 2.10.0
**Score-aware taxonomic interpretation**
- Replaced flat Top-N voting with a primary high-scoring BLAST cluster.
- Added the best alternative species and genus, together with ΔBit from the leading candidate.
- Added rank fallback so uncertain species evidence can still produce a supported genus identification.
- Added reference context for type material, RefSeq, ordinary GenBank records and unresolved/environmental annotations.
- Added the BLAST score landscape and expanded the exported decision evidence.
- Improved accession-level HSP handling so one accession is treated as one database hit.

## 2.9.0
**Save / Load projects and BLAST hit normalization**
- Added portable `.sangerproject` files that preserve processed traces, QC, rename state, BLAST jobs/results and taxonomy analyses.
- Loading a project restores the workflow without repeating completed trimming or BLAST work.
- Multiple HSPs from the same accession are aggregated into one hit before taxonomic analysis.
- Sequence evidence uses competitive unique accessions rather than the weak tail of all retrieved rows.

## 2.8.0
**Reporting and batch taxonomy**
- Added Help / About, visible application versioning and this changelog.
- Added Analyze all retrieved sequences and a consolidated Team identification summary.
- Added CSV and Excel team reports with BLAST, taxonomy, QC, rename and run-setting sheets.
- Exports record the application version and main run settings for provenance.

## 2.7.0
**Taxonomic consensus correction**
- Separated taxonomic agreement from sequence-quality evidence.
- Weak or missing identity/coverage no longer automatically erases a strong genus/species agreement signal.
- Unresolved labels such as `sp.`, uncultured and unclassified are not treated as conflicting named species.

## 2.6.0
**Batch BLAST**
- Added submission and retrieval for all processed sequences while keeping each Sample → RID → Hit set separate.
- Increased the default number of retrieved hits and added support for up to 100.
- Added an overview of retrieved BLAST results across samples.

## 2.5.0
**Multi-hit taxonomy**
- Added NCBI taxonomy lineage resolution.
- Added genus/species agreement and lowest-common-taxon context across multiple hits.
- Added locus-aware interpretation notes.

## 2.4.0
**BLAST metadata and navigation**
- Improved NCBI metadata enrichment and BLAST hit-table readability.
- Added clearer workflow navigation and stage-position indicators.

## 2.3.0
**Interactive chromatogram and BLAST parsing**
- Moved chromatogram inspection to an interactive browser plot.
- Restored QC plots and kept them in exports.
- Improved parsing and storage of multiple BLAST hits.

## 2.2.0
**QC and chromatogram navigation**
- Improved chromatogram navigation and reduced the role of experimental primer-site mapping.
- Improved NCBI result metadata display.

## 2.1.0
**Primer-aware QC experiment**
- Added chromatogram inspection, sequence preview and experimental primer-site mapping.

## 2.0.0
**Modular Shiny application**
- Split the application into separate R source modules.
- Established the staged workflow: Upload → Assay & Trim → QC → Rename → Export → NCBI BLAST → Taxonomic summary.

## 1.x
**Initial workflow**
- Wrapped the original AB1 trimming script in Shiny.
- Added configurable trimming, rename mapping and FASTA/CSV export.
