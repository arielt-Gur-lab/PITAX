# PITAX 3.0 development roadmap

PITAX 3.0 will be built as five gated stages. A stage advances only after its tests and manual acceptance checks are green.

## Stage 1 — Read evidence foundation
Goal: validate the raw AB1 evidence model before changing any biological decision.

- Preserve v2.14.2 biological decision behavior; the deliberate Rename-before-QC workflow change is presentation/state flow only.
- Capture basecaller quality, raw ABIF primary-call coordinates and canonical A/C/G/T trace evidence.
- Validate the legacy PITAX primary positions against ABIF PLOC.2 and the inferred channel map against the canonical trace order on real AB1 files.
- Establish regression tests and an audit export.

Gate: real laboratory AB1 files produce interpretable audit evidence without processing failures; original and resolved sample identities remain traceable through Rename → QC; trimming, curation, BLAST and taxonomy decisions remain unchanged.

Status: **Closed in v3.0.0-alpha.4** after the accepted three-sample validation and Rename-before-QC regression.

## Stage 2 — Isolate / locus / read architecture
Goal: make the isolate the primary biological object while preserving single-read workflows.

- Introduce Project → Isolate → Locus → Read data objects.
- Resolve explicit read assignment (isolate, locus, direction, primer) in separate operator-edited fields; preserve the upload barcode and generate the FASTA name without parsing filenames.
- Add project-schema migration from v2 projects.
- Keep one-read-per-locus fully valid.

Gate: old projects load; new projects can represent single and paired reads without losing any v2 evidence or audit trail.

Status: **Gate candidate in v3.0.0-alpha.6**. Upload now supports explicit manual, batch and key-based assignment; PITAX generates `<Isolate>_<Locus>_<F/R>` without parsing the barcode, and single/pair counting has a regression test. Consensus remains reserved for Stage 3.

## Stage 3 — Forward / Reverse consensus
Goal: build an auditable, quality-aware consensus rather than a majority vote.

- Reverse-complement and overlap alignment.
- Per-column provenance from Forward and Reverse reads.
- Automatic resolution only when evidence is sufficiently asymmetric.
- Explicit conflict state for strong contradictory evidence.
- Side-by-side chromatogram review, manual consensus edits, undo/redo and revision binding.

Gate: curated truth-set comparisons and F/R edge-case tests pass.

## Stage 4 — Multi-locus fungal identification
Goal: interpret each locus separately, then integrate evidence at isolate level.

- Per-locus BLAST/taxonomy remains evidence-first.
- Marker profiles for major fungal groups.
- Cross-locus concordance/conflict logic; no flat locus voting.
- Reference provenance layer (type/curated/standard database records).
- Suggested next locus/action when species resolution is insufficient.

Gate: curated Fusarium, Trichoderma, Aspergillus, Penicillium, Pleurotus and general-fungi scenarios behave conservatively and explainably.

## Stage 5 — Production validation and 3.0 stable
Goal: harden the online multi-user system and validate the scientific workflow.

- Global NCBI request scheduler for multiple Shiny sessions.
- Session-isolation and concurrent-user tests.
- Curated-provider adapters where justified (e.g. RefSeq, UNITE, FUSARIUM-ID).
- Performance/project-size review.
- Real-lab truth-set validation and release documentation.

Gate: local and Posit Connect Cloud acceptance suite passes; PITAX 3.0 is released as stable.
