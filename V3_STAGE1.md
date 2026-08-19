# PITAX 3.0 — Stage 1 validation

Version: `3.0.0-alpha.3`

## Purpose

Stage 1 validates the AB1 evidence model before PITAX changes trimming or builds Forward/Reverse consensus. The established v2.14.2 processing path remains authoritative in this alpha.

alpha.3 follows the review of `265DMAA001_customer`. It keeps the corrected alpha.2 evidence model, fixes the chromatogram title/legend layout, records explicit auto-trim boundaries, and adds a same-length PCON-only comparison window. No trimming, curation, BLAST, or taxonomy decision rule is changed.

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

No Stage 1 evidence value changes the sequence in alpha.3.

## Same-length PCON comparison

alpha.3 compares the established v2 auto trim with every contiguous PCON-quality window of the same length. Candidate windows require at least 90% available PCON scores and are ranked deterministically by:

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
All PITAX v3.0.0-alpha.3 tests passed.
```

The Stage 1 test verifies PLOC.2 conversion, PCON indexing, canonical trace mapping, per-base signal evidence, auto-trim quality summaries, selected-sample export identity, window membership and the non-mutating same-length PCON proposal.

## Manual acceptance test

Use the same batch of 13 real AB1 files used for alpha.1 if possible.

### A. Regression

Confirm that:

1. Files process normally.
2. Active trim start/end and trimmed sequence are unchanged from alpha.2.
3. Chromatogram and ambiguous-peak review still work.
4. Manual curation still works.
5. Rename, Export, BLAST and Taxonomy are not blocked by the Stage 1 audit.
6. The PCON-only proposal never replaces the active trim or processed sequence.

### B. Corrected run audit

Open **3 · QC & Chromatogram → Stage 1 · AB1 evidence audit**.

Expected findings for normal ABI files:

- `Quality tag`: usually `PCON.2` if present.
- `Primary position source`: `ABIF PLOC.2 + 1`.
- `Primary positions different (%)`: normally 0% if sangerseqR and raw ABIF positions agree.
- `Maps match`: TRUE when the existing inferred map is A=1,C=2,G=3,T=4.

Quality fields are observations, not pass/fail thresholds yet.

### C. Chromatogram display regression

1. The plot title shows only the selected sample ID; it must not contain `<e2><80><94>` or a duplicate `interactive chromatogram` label.
2. The A/C/G/T and QC-marker legend is outside the data area on the right.
3. Called-base labels and QC triangles occupy separate vertical bands and do not overlap.

### D. Trim comparison and export regression

This is required for alpha.3:

1. Select `265DMAA001_customer` in the Sample dropdown.
2. Confirm the audit note says `Selected sample: 265DMAA001_customer`.
3. Confirm the run audit reports the active auto trim start, end and length. For the previously supplied 001 audit these were 28, 420 and 393.
4. In **Legacy auto trim vs PCON-only comparison**, the legacy row must say `Active output`; the PCON row must say `Observational only` and have the same length.
5. For the previously supplied 001 data, the PCON proposal is expected near bases 265–657 and should improve the displayed Q20/Q30 metrics. This is a regression expectation, not a new trimming rule.
6. Click **Download selected-base audit CSV**. Filename must start with `265DMAA001_customer`.
7. Open the CSV: every row in `Sample_ID` must identify 001; `In_auto_trim` must be TRUE only for positions 28–420; `In_quality_proposed_window` must mark the displayed PCON window.
8. Click a different row in the run-audit table. The Sample dropdown, selected-sample note and comparison table must all change to that sample.

### E. Files to send back

After the checks, send:

- the new **Run audit CSV**;
- the **selected-base audit CSV for 265DMAA001_customer**.

That is enough for the final Stage 1 evidence review. A second clean sample is optional because alpha.1 already supplied a clean example.

## Stage 1 alpha.3 is green when

- all four automated test groups pass;
- no v2 workflow regression is observed;
- raw PLOC primary positions agree with the established PITAX positions on the tested AB1 files, or any disagreement is understood;
- selected-sample export identity is correct;
- auto-trim and proposed-window memberships match their displayed boundaries;
- chromatogram title, legend and QC markers render without overlap or byte markers;
- the PCON proposal remains observational and does not change the established output;
- PCON behavior in 001 and the rest of the batch is reviewed without turning one problematic read into a universal trimming rule;
- cloud deployment is not required for the branch gate unless you intentionally deploy the v3 branch.

Do not start Stage 2 until this gate is reviewed.

## Primary implementation reference

The `sangerseqR` ABIF constructor builds `traceMatrix` in A,C,G,T order and initializes `peakPosMatrix` from primary/secondary ABIF base-call positions. `makeBaseCalls()` later rebuilds peak position/amplitude matrices as A,C,G,T per-window peaks, but PITAX does not call it in the established v2 processing path.

This distinction is now encoded explicitly in the Stage 1 audit and its regression tests.
