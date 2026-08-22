# Changelog

This file records changes that affect the way the pipeline behaves, presents results, or stores project data. Small internal refactors are omitted unless they change something visible or scientifically relevant.

## 3.0.0-alpha.10.3
**Encoding and release-process recovery**
- Removed all raw non-ASCII characters from runtime R and JavaScript source; UI punctuation and symbols now use ASCII Unicode escapes so Windows and Connect render the same text independently of locale.
- Removed source-time encoding conversion entirely because runtime source files are now ASCII-only.
- Added a contract that scans every runtime R, JavaScript and CSS source file and fails if a raw non-ASCII byte is introduced.
- Updated `Git.BAT` to fetch and detect a remote-ahead or diverged branch before tests and commit, while still retrying a previously committed local push when there are no new file changes.
- Release archives no longer include `.git`; updates must be copied over the existing Git working folder so repository history and credentials remain local and current.

## 3.0.0-alpha.10.2
**Cross-platform source-encoding startup hotfix**
- Fixed the Connect Cloud startup failure at `assay_profiles.R:35` caused by forcing UTF-8 conversion after the Linux worker fell back to the `C` locale.
- Added a platform-safe source helper: Windows explicitly decodes source files as UTF-8, while Linux/Connect retains the source bytes without locale-dependent conversion.
- Replaced the Greek alpha/beta aliases in the new assay-domain source with ASCII Unicode escapes, preserving their runtime values without requiring non-ASCII parsing during startup.
- Updated the structure contracts to protect both the Windows mojibake fix and the Linux/Connect startup path.

## 3.0.0-alpha.10.1
**Alpha 10 foundation: schema 6, controlled loci and assay provenance**
- Added the schema-6 assay domain with a controlled locus vocabulary; removed the free-text `Other` locus option from the assay screen.
- Split assay-specific fields (locus, primers and length limits) from project-level trimming defaults.
- Added persistent `assay_profiles` and `project_defaults` state and linked every read to an explicit `Assay_ID`.
- Made read locus and direction-specific primer provenance resolve from the assigned assay profile instead of editable locus text.
- Added the direct schema-5 to schema-6 migration: the former run-level settings become one assay and all existing reads are linked to it without changing read, BLAST or taxonomy evidence.
- Blocked speculative loading of project schemas 1–4; those projects must first be opened and resaved with the frozen Alpha 9.1 application.
- Added assay metadata to project architecture and checkpoint exports.
- Fixed Windows mojibake such as `<c2><b7>` and `<e2><80><94>` by loading every application and server module explicitly as UTF-8.
- Expanded the Windows test runner to 15 groups with multi-assay, controlled-vocabulary, dangling-reference and schema-5 migration coverage.

## 3.0.0-alpha.9.2
**Project structure and development continuity**
- Reduced `app.R` to a small application entry point and separated bootstrap, UI, server, browser assets, domain logic, services and exports by responsibility.
- Grouped automated tests into unit, integration and contract suites, and moved controlled AB1 data under `tests/fixtures`.
- Moved documentation, deployment utilities, static assets and reusable templates into dedicated directories.
- Updated all source paths, download paths and the Windows test runner to the organized layout.
- Added `docs/PITAX_MASTER.md` as the authoritative development guide for product intent, scientific invariants, data dependencies, lessons learned and the agreed roadmap.
- Added a fourteenth contract test that protects the composition root, module locations, external browser assets and server-stage load order.
- Kept project schema 5 and all processing, consensus, BLAST, taxonomy and multi-locus behavior unchanged.

## 3.0.0-alpha.9.1
**Multi-locus visual interpretation and global table alignment**
- Added a visual Stage 4 interpretation workspace with isolate selection, project-wide overview counts, a status/conclusion panel, Identity/coverage comparison and one evidence card per locus.
- Kept the complete source, isolate-profile and per-locus audit tables below the new visual layer.
- Kept multiple isolates strictly separated in the visual selector while showing all imported loci for the selected isolate.
- Fixed horizontal DataTables header/body drift by removing the global forced-width override, enabling automatic width calculation for scrolling tables, recalculating visible tables after tab/data changes and synchronizing header/body scrolling.
- Added multi-isolate visual-selection tests and a thirteenth regression group covering the global DataTables alignment contract.
- No project schema or scientific interpretation rules changed in this maintenance release.

## 3.0.0-alpha.9
**Stage 4 multi-locus isolate profile foundation**
- Added a ninth workflow screen that imports separately completed single-locus PITAX projects and joins evidence only by the explicit Isolate field.
- Added project schema 5 with a persistent `pitax-multilocus-profile-v1` snapshot while preserving all Stage 3, BLAST and taxonomy evidence during migration.
- Retained one evidence row per Isolate/Locus with sequence, consensus revision, source-project fingerprint, BLAST RID, taxonomic interpretation, locus limitation and reference context.
- Blocked duplicate Isolate/Locus evidence and repeated source files.
- Added conservative profile states for concordant species, concordant genus, missing taxonomy, partial evidence, species conflict and genus conflict.
- Explicitly prohibited flat locus voting: even a 2:1 majority cannot erase a conflicting genus call.
- Added current-session staleness detection and disabled Stage 4 downloads until a changed current project is rebuilt into the profile.
- Added profile/evidence CSV, multi-locus FASTA and Checkpoint F ZIP exports.
- Expanded the Windows test runner from 10 to 12 groups with Stage 4 logic and Shiny integration contracts.
- Stage 3 controlled-fixture validation is accepted for continued development; independently sequenced paired-AB1 validation remains deferred.

## 3.0.0-alpha.8.2
**Shiny reactive-context startup hotfix**
- Removed the hidden `rv$project_mode` read from the navigation helper and now pass project mode explicitly from reactive observers.
- Wrapped the one-time post-flush project-mode read in `isolate()`, preventing Shiny from rejecting reactive-value access outside a reactive consumer.
- Normalized loaded Simple projects away from the hidden Consensus tab.
- Added a regression contract for the post-flush callback and audited every alpha.8 helper that reads `input` or `rv` for its call context.

## 3.0.0-alpha.8.1
**Plotly namespace startup hotfix**
- Qualified both new consensus chromatogram UI outputs with `plotly::`, preventing application startup failure when Plotly is installed but not attached with `library(plotly)`.
- Added a Stage 3 integration contract that requires the namespace-qualified UI and server bindings.

## 3.0.0-alpha.8
**Consensus review, Simple-mode bypass and stable analysis table**
- Simple projects now materialize correctly oriented independent analysis sequences automatically after QC and continue directly to Export; the user-facing consensus tab is omitted.
- Stabilized the Analysis sequence summary with a compact fixed-width column contract so sorting no longer changes the layout.
- Added a dedicated Consensus Review & Curation workspace for paired projects with Forward, Reverse and intentional IUPAC decisions.
- Linked each selected conflict to focused Forward and Reverse chromatograms at the original called-base positions.
- Added monotonic consensus revisions, persistent audit history and auditable Undo/Redo actions while retaining immutable source-read revision bindings.
- Bound every submitted BLAST RID to the exact consensus revision used as its query.
- Exported consensus curation audit CSVs beside per-column evidence.
- Added a deterministic one-conflict AB1 pair under `TEST/Stage3_conflict_pair` and regression coverage for Forward/IUPAC acceptance, gate opening and Undo/Redo.
- The controlled mirror-derived fixtures validate software mechanics but do not replace an independently sequenced biological Forward/Reverse acceptance set.

## 3.0.0-alpha.7.2
**Wider project-mode choices and controlled Stage 3 AB1 fixture**
- Expanded the Project read model radio group to the full available horizontal width without increasing card height.
- Added a reproducible synthetic Forward/Reverse AB1 pair, assignment key and expected sequence under `TEST/Stage3_synthetic_pair`.
- The Reverse fixture is a trace-aware reverse complement of a permissively licensed public AB1 test trace: calls, quality values, peak positions and A/C/G/T channels are transformed consistently.
- Documented the fixture boundary: it validates mechanics but cannot replace a real independently sequenced Forward/Reverse pair for closing the Stage 3 gate.

## 3.0.0-alpha.7.1
**Rename workspace, project read modes and BLAST depth clarification**
- Removed biological identity editing from Upload and concentrated assignment-key import, batch editing, manual editing and generated-name review in Rename before Trim/QC.
- Replaced the reactive DT identity editor with a stable form table and explicit **Apply table changes** action, preventing full-table redraw after each edit.
- Made Direction a native Forward/Reverse selector in every row.
- Added a saved project read model: Simple independent reads or explicit Forward/Reverse pairing.
- Simple mode keeps each read independent, retains its `<Isolate>_<Locus>_<F/R>` name and orients Reverse reads without attempting a merge.
- Relabeled BLAST hit count as retrieval depth and added audit fields separating retrieved, competitive and noncompetitive tail hits.
- Removed accession-count boosts from confidence; database representation remains audit context and no longer acts as a vote in either direction.
- Added regression coverage proving that many weak tail hits do not lower the recommendation or confidence; a genuinely close late hit may still do so.
- Limited curation invalidation notices to cases where downstream BLAST/taxonomy evidence really existed and included the curation action plus sequence revision transition.

## 3.0.0-alpha.7
**Stage 3 Forward / Reverse consensus foundation**
- Closed the accepted Stage 2 gate and advanced the project schema to version 4.
- Enforced one Gene/Locus per uploaded sequencing run; Stage 4 will combine separately processed locus profiles.
- Added a derived Reverse-complement view and deterministic semi-global overlap alignment without modifying either curated source read.
- Added one isolate-level sequence record per Isolate/Locus, while labeling one-read cases as single-read representatives rather than consensus.
- Added conservative quality-aware mismatch resolution; unresolved contradictions retain an IUPAC ambiguity call and require review.
- Blocked weak/unreliable overlaps instead of concatenating reads.
- Added per-column provenance with Forward/Reverse calls, basecaller quality, raw chromatogram positions, decision and review state.
- Routed downstream FASTA and BLAST through the isolate-level sequence after the Stage 3 gate is green.
- Added consensus staleness and invalidation when a curated source sequence or rebuilt isolate-level sequence changes.
- Added consensus FASTA, summary, alignment and per-column evidence checkpoint exports.
- Added schema-3 → schema-4 migration that preserves established read, QC, curation, BLAST and taxonomy evidence.
- Added Stage 3 algorithm and Shiny integration regressions; manual conflict editing and paired-AB1 truth-set validation remain open.

## 3.0.0-alpha.6
**Stage 2 workflow and identity correction**
- Separated the user-facing stages into Upload → Assay → Rename & Assign → Trim & QC.
- Removed trimming from the Assay-to-Rename transition, so Rename opens immediately and processing starts only with **Start trimming**.
- Added explicit Isolate, Gene/Locus and Direction editing directly on Upload, with expanded batch tools and an XLSX/CSV assignment key.
- Kept the upload barcode immutable and removed filename/name parsing as a source of biological identity.
- Generated `<Isolate>_<Locus>_<F/R>` from the explicit fields and stored both fields and generated name in project schema version 3.
- Placed the generated-name assignment table before the architecture preview on Rename review.
- Added a pre-trim architecture checkpoint and schema-2 → schema-3 migration.
- Added regressions proving one read remains a single read and only distinct Forward plus Reverse reads form a pair.
- Prioritized duplicate source/read-key validation so unsafe duplicate AB1 basenames are reported even while biological assignments are still incomplete.
- Guarded empty pre-trim and all-failed-run selector states so QC/BLAST choice synchronization never calls `setNames()` on `NULL`.
- Replaced free-text Direction cell editing with per-read Forward/Reverse dropdowns in both Upload assignment and Rename review.
- Did not add Forward/Reverse consensus; that remains gated to Stage 3.

## 3.0.0-alpha.5
**Stage 2 isolate / locus / read architecture**
- Started Stage 2 with a Project → Isolate → Locus → Read architecture and project schema version 2.
- Added an editable read-assignment table with conservative filename inference and visible assay-default provenance.
- Added separate per-read assignment and processing-settings provenance so mixed locus/direction metadata is not collapsed into one run-wide label.
- Added in-memory migration for schema-1 projects while preserving established processing, curation, BLAST and taxonomy state.
- Added unique default names for paired reads and retained one-read-per-locus as fully valid.
- Added read-assignment and architecture CSVs to processing checkpoints.
- Added Stage 2 architecture, migration and Shiny integration regression tests.
- Did not add Forward/Reverse consensus; that remains gated to Stage 3.

## 3.0.0-alpha.4
**Stage 1 gate-closing workflow order**
- Moved Rename before QC in the user-facing workflow: Upload → Assay & Trim → Rename → QC → Export.
- Trimming now opens Rename; validated names open QC; QC continues to Export.
- QC selectors, trimming tables, chromatogram titles and audit displays use resolved names while audit CSVs preserve both original `Sample_ID` and `Final_Name`.
- Reordered checkpoints: Checkpoint A stores renamed sequences; Checkpoint B stores the renamed and curated QC state.
- Simplified the Plotly range navigator by clipping compressed base letters and QC markers from its overview.
- Added an automated workflow-order and QC-naming regression test.
- Retained the validated alpha.3 AB1 evidence model and observational PCON comparison without changing trimming, curation, FASTA, BLAST or taxonomy decisions.

## 3.0.0-alpha.3
**Stage 1 trim comparison and chromatogram layout correction**
- Added explicit auto-trim start, end and length to the run-level AB1 audit.
- Added `In_auto_trim` to every per-base audit export row.
- Added a same-length PCON-only comparison window ranked by Q>=20, Q>=30, median quality, mean quality and coverage, with deterministic proximity tie-breaking.
- Added a selected-sample comparison table for the active legacy trim and the observational PCON window.
- Added `In_quality_proposed_window` to the per-base audit table and CSV.
- Kept the quality-window proposal strictly observational; it does not modify sequence processing, FASTA, curation, BLAST or taxonomy.
- Simplified the Plotly title to the sample ID to avoid Windows byte-marker rendering of the em dash.
- Moved the chromatogram legend outside the data area and separated called-base labels from QC markers.

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
