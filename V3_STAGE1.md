# PITAX 3.0 — Stage 1 validation

Version: `3.0.0-alpha.4`

## Purpose

Stage 1 validates the AB1 evidence model before PITAX changes trimming or builds Forward/Reverse consensus. The established v2.14.2 processing path remains authoritative in this alpha.

alpha.4 is the gate-closing build. It keeps the validated alpha.3 evidence model unchanged, moves Rename before QC, carries resolved names into the QC display, reorders checkpoints, and simplifies the chromatogram navigator. No trimming, curation, BLAST, or taxonomy decision rule is changed.

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

No Stage 1 evidence value changes the sequence in alpha.4.

## Same-length PCON comparison

alpha.4 retains the alpha.3 comparison of the established v2 auto trim with every contiguous PCON-quality window of the same length. Candidate windows require at least 90% available PCON scores and are ranked deterministically by:

1. Q>=20 percentage;
2. Q>=30 percentage;
3. median quality;
4. mean quality;
5. quality coverage;
6. smallest shift from the established auto-trim start.

The selected comparison window is marked **Observational only**. It does not change the active trim, processed sequence, FASTA, curation, BLAST query or taxonomy.

## Export fix

The selected sample is now one source of truth for QC and audit export.

- Clicking a row in the run-audit table selects that sample in the QC workspace.
- The selected-sample note shows the sample ID explicitly.
- Selected-base CSV includes a `Sample_ID` column on every row.
- QC displays the resolved final name; audit exports preserve both `Sample_ID` and `Final_Name`.
- Selected-base CSV includes `In_auto_trim` and `In_quality_proposed_window` on every evidence row.
- Filename is derived from the same result object used to create the CSV.
- If selected key and result sample ID disagree, PITAX blocks the download rather than producing a misleading file.

## Automated test

Run:

```text
run_tests.bat
```

Expected final line:

```text
All PITAX v3.0.0-alpha.4 tests passed.
```

Five test groups verify taxonomy logic, ambiguous-peak logic, manual curation, the AB1 evidence model, and the Rename-before-QC workflow contract.

## Manual acceptance test

The user accepted the three-sample batch (000, 001 and 002) as sufficient for the Stage 1 gate.

### A. Regression

Confirm that:

1. Files process normally.
2. Active trim start/end and trimmed sequence are unchanged from alpha.2.
3. Chromatogram and ambiguous-peak review still work.
4. Manual curation still works.
5. Rename, QC, Export, BLAST and Taxonomy are not blocked by the Stage 1 audit.
6. The PCON-only proposal never replaces the active trim or processed sequence.

### B. Corrected run audit

Open **4 · QC & Chromatogram → Stage 1 · AB1 evidence audit**.

Expected findings for normal ABI files:

- `Quality tag`: usually `PCON.2` if present.
- `Primary position source`: `ABIF PLOC.2 + 1`.
- `Primary positions different (%)`: normally 0% if sangerseqR and raw ABIF positions agree.
- `Maps match`: TRUE when the existing inferred map is A=1,C=2,G=3,T=4.

Quality fields are observations, not pass/fail thresholds yet.

### C. Chromatogram display regression

1. The plot title shows the resolved sample name; it must not contain `<e2><80><94>` or a duplicate `interactive chromatogram` label.
2. The A/C/G/T and QC-marker legend is outside the data area on the right.
3. Called-base labels and QC triangles occupy separate vertical bands and do not overlap.
4. The bottom navigator shows a compact signal overview without compressed base letters or QC markers.

### D. Rename-before-QC regression

This is required for alpha.4:

1. Run trimming; PITAX must continue to **3 · Rename**.
2. Apply a rename key or batch edit and confirm names are valid.
3. Continue to **4 · QC & Chromatogram**.
4. The QC sample selector, trimming table, chromatogram title and audit display must show the resolved name.
5. The audit CSV must preserve the original `Sample_ID` and include the resolved `Final_Name`.
6. Back/Continue navigation must follow Rename → QC → Export.
7. Checkpoint A must be the renamed state; Checkpoint B must include the renamed and curated QC state.

### E. Trim comparison regression

1. Select the resolved name corresponding to `265DMAA001_customer` in the Sample dropdown.
2. Confirm the audit note shows the resolved name and preserves the original ID.
3. Confirm the run audit reports the active auto trim start, end and length. For the previously supplied 001 audit these were 28, 420 and 393.
4. In **Legacy auto trim vs PCON-only comparison**, the legacy row must say `Active output`; the PCON row must say `Observational only` and have the same length.
5. For the previously supplied 001 data, the PCON proposal is expected near bases 265–657 and should improve the displayed Q20/Q30 metrics. This is a regression expectation, not a new trimming rule.
6. Click **Download selected-base audit CSV**. The filename uses the resolved name.
7. Open the CSV: `Sample_ID` must preserve 001, `Final_Name` must contain the resolved name, `In_auto_trim` must be TRUE only for positions 28–420, and `In_quality_proposed_window` must mark positions 265–657.
8. Click a different row in the run-audit table. The Sample dropdown, selected-sample note and comparison table must all change to that sample.

### F. Accepted Stage 1 evidence

The accepted real-data validation consists of:

- a three-sample run audit for 000, 001 and 002;
- a detailed 000 audit validating both 500-base memberships;
- a detailed 001 audit validating active positions 28–420 and proposed positions 265–657;
- successful alpha.3 automated tests;
- chromatogram confirmation that the title and external legend were corrected.

No additional AB1 batch is required for this gate.

## Stage 1 alpha.4 is green when

- all five automated test groups pass;
- no v2 workflow regression is observed;
- raw PLOC primary positions agree with the established PITAX positions on the tested AB1 files, or any disagreement is understood;
- selected-sample export identity is correct;
- auto-trim and proposed-window memberships match their displayed boundaries;
- chromatogram title, legend and QC markers render without overlap or byte markers;
- the PCON proposal remains observational and does not change the established output;
- PCON behavior in 001 and the accepted three-sample batch is reviewed without turning one problematic read into a universal trimming rule;
- Rename precedes QC and resolved names appear consistently throughout QC;
- the chromatogram navigator no longer displays compressed base letters or QC markers;
- cloud deployment is not required for the branch gate unless you intentionally deploy the v3 branch.

After the alpha.4 workflow and navigator checks pass, Stage 1 is closed and Stage 2 may begin.

## Primary implementation reference

The `sangerseqR` ABIF constructor builds `traceMatrix` in A,C,G,T order and initializes `peakPosMatrix` from primary/secondary ABIF base-call positions. `makeBaseCalls()` later rebuilds peak position/amplitude matrices as A,C,G,T per-window peaks, but PITAX does not call it in the established v2 processing path.

This distinction is now encoded explicitly in the Stage 1 audit and its regression tests.
