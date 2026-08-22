# PITAX 3.0 — Stage 3 consensus foundation

Retained in version: `3.0.0-alpha.9.1`

## Purpose

Stage 3 creates one downstream sequence per Isolate–Locus without destroying or rewriting either source read. A two-read consensus is treated as an evidence object, not as a text merge.

## Scientific and data contract

- The Forward curated sequence stays in sequencing orientation.
- The Reverse curated sequence is reverse-complemented only in a derived alignment view.
- Pairing comes only from the explicit Isolate, Locus and Direction fields accepted in Stage 2.
- A run contains exactly one locus. Stage 4 will combine separately processed locus results at isolate level.
- A single read is valid but is explicitly reported as a single-read representative.
- A pair must pass minimum overlap and overlap-identity rules.
- Agreement is accepted; a canonical base resolves an `N`; a quality-asymmetric mismatch may be resolved automatically.
- A mismatch without strong asymmetric evidence is retained as an International Union of Pure and Applied Chemistry (IUPAC) ambiguity code and marked for review.
- Internal alignment gaps and weak overlaps block downstream analysis.
- Every call retains source-read and raw chromatogram positions where available.

The default thresholds in alpha.8 are conservative application review settings, not universal biological constants:

| Setting | Default |
| --- | ---: |
| Minimum overlap | 40 bp |
| Minimum overlap identity | 85% |
| Strong basecaller quality | 20 |
| Minimum quality advantage | 10 |

## Output object

Project schema 5 retains the `pitax-consensus-set-v1` object introduced in schema 4, containing:

- one record per Isolate–Locus;
- source read identifiers and exact source sequences used;
- Forward/Reverse oriented alignment;
- isolate-level sequence and status;
- per-column provenance and review state;
- algorithm/settings/timestamp;
- immutable source curation revisions;
- monotonic consensus revision, audit history and Undo/Redo stacks.

If a curated source sequence changes, the prior isolate-level sequence is stale and must be rebuilt. Rebuilding a changed sequence invalidates downstream BLAST and taxonomy evidence.

## Current alpha boundary

Implemented through alpha.8:

- orientation normalization;
- deterministic pairwise overlap alignment;
- conservative automatic call rules;
- single-read representatives;
- weak-overlap and unresolved-conflict gates;
- per-column audit display and checkpoint exports;
- schema-3 to schema-4 migration;
- downstream FASTA/BLAST routing through the isolate-level sequence.
- automatic Simple-mode orientation and direct QC → Export navigation;
- linked Forward/Reverse chromatograms at a selected conflict;
- manual Forward, Reverse or intentional IUPAC decisions;
- monotonic consensus revisions with auditable Undo/Redo;
- BLAST RID binding to the exact consensus revision;
- controlled clean and single-conflict AB1 fixtures.

Deferred biological validation items:

- real independently sequenced paired-AB1 truth-set validation;
- edge-case acceptance for indels, poor overlap boundaries and multiple conflicts.

## Automated test

Run `run_tests.bat` on Windows with R installed.

Expected final line:

```text
All PITAX v3.0.0-alpha.9.1 tests passed.
```

The full alpha.9.1 runner has thirteen groups. The first ten retain taxonomy, QC, curation, Stage 1 evidence, workflow order, Stage 2 architecture, Stage 3 algorithm, AB1 fixtures and Shiny integration coverage; groups 11–12 cover Stage 4 and group 13 protects global table alignment.

## Focused manual acceptance for alpha.8

### A. Single-locus run gate

1. Assign all uploaded reads to `ITS`; confirm trimming can start.
2. Change one row to `TEF1`; confirm Start trimming is blocked with a one-locus message.
3. Restore `ITS` and trim the run.

### B. Single read

1. Process one Forward read and build Stage 3.
2. Confirm status is `SINGLE_READ` and the output name is `<Isolate>_<Locus>`.
3. Repeat with one Reverse read and confirm the displayed/output sequence is reverse-complemented while the QC source sequence is unchanged.

### C. Clean pair

1. Process a known Forward/Reverse pair for the same isolate and ITS locus.
2. Build Stage 3 and inspect oriented alignment, overlap length and identity.
3. Confirm the summary contains one isolate-level row, not two read rows.
4. Confirm the gate is green only if there are no unresolved positions.
5. Download the isolate-level FASTA and verify one record named `<Isolate>_ITS`.

### D. Conflict and weak overlap

1. Upload `tests/fixtures/Stage3_conflict_pair` and import its assignment key in Paired mode.
2. Confirm PITAX retains IUPAC `M`, marks the column for review and blocks Export.
3. Select the conflict and confirm both chromatograms focus on their linked raw positions.
4. Accept Forward, Reverse or IUPAC; confirm Revision increments and Export opens.
5. Undo and Redo; confirm both actions increment Revision and remain in History.
6. Raise the minimum overlap above the available overlap; rebuild and confirm `NO_RELIABLE_OVERLAP` with no emitted sequence.

### E. Audit, staleness and migration

1. Download the Stage 3 checkpoint and confirm it includes summary, FASTA, alignment and per-column evidence files.
2. Change a curated source sequence; confirm Stage 3 must be rebuilt.
3. Load an alpha.6 project; confirm read/QC/BLAST/taxonomy evidence remains present and Stage 3 starts empty.
4. Save and reload; confirm the project now uses schema 5 and retains the Stage 3 object.
