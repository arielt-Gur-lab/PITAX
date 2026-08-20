# PITAX 3.0 — Stage 2 closed gate

Version: `3.0.0-alpha.6`

## Purpose

alpha.6 makes the sequencer filename/barcode an immutable technical identifier and makes biological identity an explicit operator assignment. PITAX does not infer isolate, gene or direction by parsing the uploaded name.

The visible workflow is:

**Upload & assign → Assay → Rename review → Trim & QC → Export**

The Project → Isolate → Locus → Read architecture is additive. It does not merge Forward/Reverse reads, change trimming decisions, or reinterpret accepted Stage 1 evidence.

## Identity contract

For every source read the operator assigns three independent fields:

| Field | Meaning |
| --- | --- |
| Isolate | Laboratory isolate/sample code |
| Gene / locus | Sequenced marker |
| Direction | Forward or Reverse |

PITAX then generates the final read/FASTA name:

```text
<Isolate>_<Locus>_<F/R>
```

For example, explicit fields `FB120`, `ITS`, `Forward` generate `FB120_ITS_F`.

The generated name is an output and display label. The separate fields remain the source of truth, so PITAX never needs to split that name back into biological metadata. The original upload barcode remains stored alongside both.

## Upload assignment interface

The Upload screen supports:

- direct cell editing for Isolate and Gene/Locus;
- a closed Forward/Reverse dropdown for Direction, with a non-valid `Select…` placeholder until assignment;
- multi-row selection;
- batch isolate prefix, suffix, find and replace;
- batch Gene/Locus and Direction assignment;
- XLSX/CSV key import using `old_id`, `isolate`, `locus` (or `gene`) and `direction`;
- automatic generation and validation of the final read/FASTA name.

`assignment_key_template.csv` is included as a ready-to-edit example.

If no rows are selected, a batch action applies to all uploaded reads. Source barcode and generated name columns are read-only.

## Data model

| Object | Meaning | Stable link |
| --- | --- | --- |
| Project | Saved PITAX analysis session | Project schema version 3 |
| Isolate | Explicit biological isolate/sample | `Isolate_ID` |
| Locus | One genetic marker belonging to one isolate | `Locus_ID` → `Isolate_ID` |
| Read | One immutable AB1 source read | `Read_ID` → `Locus_ID`; `Source_ID` preserves the upload barcode |

One read per locus remains fully valid. A pair is counted only when the same isolate/locus contains two distinct source reads explicitly assigned Forward and Reverse. Stage 2 never creates a consensus sequence.

## Workflow behavior

1. **Upload** preserves source barcodes and provides the full assignment interface.
2. **Assay** records primer defaults and trimming parameters. **Continue to Rename** changes screens immediately and performs no trimming.
3. **Rename & Assign** is the final editable review of explicit fields and generated names.
4. **Start trimming** validates the fields, stores the architecture and begins AB1 processing.
5. **Trim & QC** retains the established chromatogram, audit and curation behavior.

Checkpoint A is available after valid assignment and before trimming. It contains the assignment key, explicit read assignments, generated names, architecture tables and run settings.

## Project migration

- Schema 1 projects migrate through the existing evidence-preserving path and are represented as schema 3 in memory.
- Schema 2 projects receive schema-3 fields without parsing or rewriting their existing output names, processing, curation, BLAST or taxonomy evidence.
- A new trim run requires complete explicit Isolate, Gene/Locus and Direction fields.
- Migration never overwrites the source project automatically; saving creates a schema-3 project.

## Automated test

Run `run_tests.bat` on Windows with R installed.

Expected final line:

```text
All PITAX v3.0.0-alpha.6 tests passed.
```

Seven groups cover the accepted taxonomy, QC, curation, Stage 1 evidence, workflow-order, Stage 2 architecture/migration and Shiny integration contracts.

## Focused manual acceptance

### A. Upload and assignment

1. Upload the accepted AB1 batch with generic sequencer/barcode names.
2. Confirm the source column preserves those names unchanged and no isolate is inferred from them.
3. Edit one row manually: Isolate `FB120`, Gene `ITS`, Direction `Forward`.
4. Confirm PITAX generates `FB120_ITS_F`.
5. Select multiple rows and test isolate prefix/suffix, Gene and Direction batch actions.
6. Import a key with `old_id,isolate,locus,direction` and confirm only uniquely matched rows update.

### B. Transition and validation

1. Continue to Assay and click **Continue to Rename**.
2. Confirm Rename opens immediately without AB1 processing.
3. Confirm the generated-name table appears before the architecture preview.
4. Confirm missing isolate, gene, direction or duplicate generated names block **Start trimming**.
5. Click **Start trimming** and confirm processing opens **Trim & QC**.

### C. Single and paired reads

1. Assign one read as `FB121`, `LSU`, `Forward`; confirm one read and zero pairs.
2. Assign two distinct source reads to the same isolate/locus, one Forward and one Reverse; confirm two reads and one pair.
3. Confirm both reads remain separately selectable and no merged or consensus sequence is created.

### D. Audit and migration

1. Download Checkpoint A and confirm it contains `rename_map.csv`, `read_assignments.csv`, `project_architecture/isolates.csv`, `loci.csv`, `reads.csv`, and `run_settings.txt`.
2. After trimming, confirm source barcode, explicit identity and generated name remain traceable.
3. Load a schema-1 or schema-2 project; compare result, rename, BLAST and taxonomy row counts.
4. Save under a new filename, reload it, and confirm it loads directly as schema 3.

## Gate boundary

Stage 2 was accepted and closed after the alpha.6 automated and focused manual checks. Reverse-complement alignment, overlap conflict resolution and Forward/Reverse consensus remain outside this layer and are implemented only by the explicit Stage 3 sequence object.
