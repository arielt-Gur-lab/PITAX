# PITAX 3.0 — Stage 1 validation

Version: `3.0.0-alpha.2`

## Purpose

Stage 1 validates the AB1 evidence model before PITAX changes trimming or builds Forward/Reverse consensus. The established v2.14.2 processing path remains authoritative in this alpha.

alpha.2 follows the first real-lab audit. It corrects the interpretation of `sangerseqR` objects and fixes selected-sample CSV export synchronization. No trimming, curation, BLAST, or taxonomy decision rule is changed.

## What alpha.1 taught us

The first real batch showed that PCON.2 quality was available across the tested AB1 files and that the inferred PITAX channel map agreed with A/C/G/T. Reviewing the `sangerseqR` constructor also clarified an important detail: for freshly read ABIF files, `peakPosMatrix[,1]` is the primary ABIF base-call position (`PLOC.2 + 1`). The remaining raw `peakPosMatrix`/`peakAmpMatrix` columns must not be interpreted as A/C/G/T peak columns unless `makeBaseCalls()` has rebuilt them.

Therefore alpha.2 removes the invalid alpha.1 A/C/G/T peak-matrix comparison.

## What the corrected audit captures

For every newly processed AB1, when available:

- ABIF basecaller quality (`PCON.2`, then `PCON.1` fallback).
- Raw ABIF primary base-call position (`PLOC.2 + 1`).
- The established PITAX primary position (`peakPosMatrix[,1]`).
- Canonical A/C/G/T trace evidence from `traceMatrix`.
- Existing PITAX inferred channel map versus canonical A/C/G/T map.
- Per-base called signal, strongest alternative signal, called/alternative ratio, and whether the called base is the strongest trace channel.
- Auto-trim quality summaries: median quality, Q>=20 %, Q>=30 %, call-is-max %, and median called/alternative ratio.

No Stage 1 evidence value changes the sequence in alpha.2.

## Export fix

The selected sample is now one source of truth for QC and audit export.

- Clicking a row in the run-audit table selects that sample in the QC workspace.
- The selected-sample note shows the sample ID explicitly.
- Selected-base CSV includes a `Sample_ID` column on every row.
- Filename is derived from the same result object used to create the CSV.
- If selected key and result sample ID disagree, PITAX blocks the download rather than producing a misleading file.

## Automated test

Run:

```text
run_tests.bat
```

Expected final line:

```text
All PITAX v3.0.0-alpha.2 tests passed.
```

The Stage 1 test now verifies PLOC.2 conversion, PCON indexing, canonical trace mapping, per-base signal evidence, auto-trim quality summaries, and selected-sample export identity.

## Manual acceptance test

Use the same batch of real AB1 files used for alpha.1 if possible.

### A. Regression

Confirm that:

1. Files process normally.
2. Trim start/end and trimmed sequence are unchanged from v2.14.2/alpha.1.
3. Chromatogram and ambiguous-peak review still work.
4. Manual curation still works.
5. Rename, Export, BLAST and Taxonomy are not blocked by the Stage 1 audit.

### B. Corrected run audit

Open **3 · QC & Chromatogram → Stage 1 · AB1 evidence audit**.

Expected findings for normal ABI files:

- `Quality tag`: usually `PCON.2` if present.
- `Primary position source`: `ABIF PLOC.2 + 1`.
- `Primary positions different (%)`: normally 0% if sangerseqR and raw ABIF positions agree.
- `Maps match`: TRUE when the existing inferred map is A=1,C=2,G=3,T=4.

Quality fields are observations, not pass/fail thresholds yet.

### C. Selection/export regression

This is required for alpha.2:

1. Select `265DMAA001_customer` in the Sample dropdown.
2. Confirm the audit note says `Selected sample: 265DMAA001_customer`.
3. Click **Download selected-base audit CSV**.
4. Filename must start with `265DMAA001_customer`.
5. Open the CSV: every row in `Sample_ID` must be `265DMAA001_customer`.
6. Then click a different row in the run-audit table. The Sample dropdown and selected-sample note must change to that same sample.

### D. Files to send back

After the checks, send:

- the new **Run audit CSV**;
- the **selected-base audit CSV for 265DMAA001_customer**.

That is enough for the final Stage 1 evidence review. A second clean sample is optional because alpha.1 already supplied a clean example.

## Stage 1 alpha.2 is green when

- all four automated test groups pass;
- no v2 workflow regression is observed;
- raw PLOC primary positions agree with the established PITAX positions on the tested AB1 files, or any disagreement is understood;
- selected-sample export identity is correct;
- PCON behavior in the problematic sample is reviewed;
- cloud deployment is not required for the branch gate unless you intentionally deploy the v3 branch.

Do not start Stage 2 until this gate is reviewed.

## Primary implementation reference

The `sangerseqR` ABIF constructor builds `traceMatrix` in A,C,G,T order and initializes `peakPosMatrix` from primary/secondary ABIF base-call positions. `makeBaseCalls()` later rebuilds peak position/amplitude matrices as A,C,G,T per-window peaks, but PITAX does not call it in the established v2 processing path.

This distinction is now encoded explicitly in the Stage 1 audit and its regression tests.
