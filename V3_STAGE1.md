# PITAX 3.0 — Stage 1 validation

Version: `3.0.0-alpha.1.1`

## Purpose

Stage 1 validates the AB1 evidence model before PITAX changes trimming or builds Forward/Reverse consensus. The existing v2.14.2 processing path remains authoritative in this alpha.

This means that a successful Stage 1 test should give the same trimming, curation, BLAST and taxonomy behavior you already know, plus a new audit panel in QC.

## What the new audit captures

For every newly processed AB1 file PITAX now stores, when available:

- ABIF basecaller quality values (`PCON.2`, then `PCON.1` as fallback).
- The complete `sangerseqR` `peakPosMatrix`.
- The complete `sangerseqR` `peakAmpMatrix`.
- The established PITAX v2 peak-position vector (`peakPosMatrix[,1]`).
- A called-base-specific peak-position vector: A uses the A column, C the C column, G the G column and T the T column.
- Signal/competition evidence at both position models.
- A direct comparison between the inferred v2 channel map and the canonical A/C/G/T model.

No Stage 1 evidence value is used to change a base or trim boundary in this alpha.

## Automated test

Run:

```text
run_tests.bat
```

Expected final line:

```text
All PITAX v3.0.0-alpha.1.1 tests passed.
```

There are now four test groups; the fourth is the Stage 1 AB1 evidence helper test.

## Manual acceptance test — real AB1 files

Use a small but varied batch, ideally 6–12 files:

- at least two clean reads;
- at least two reads with weak/noisy ends;
- at least one read with visible competing/double peaks;
- if available, files from more than one sequencing run.

### A. Regression check

Process the files with the same assay/trimming settings used in v2.14.2.

Confirm:

1. All files that processed in v2.14.2 still process.
2. Trim start/end and trimmed sequence look unchanged.
3. Chromatogram, ambiguous-peak review and manual curation still work.
4. Existing Pleurotus/taxonomy behavior is unchanged.
5. No new Stage 1 warning prevents continuing to Rename, Export or BLAST.

### B. Open the new audit

In **3 · QC & Chromatogram**, open:

**Stage 1 · AB1 evidence audit → Open evidence audit**

For each sample check:

1. `Evidence` should normally be `Captured`.
2. `Quality tag` should usually be `PCON.2` on ABI files that contain basecaller confidence values. `Unavailable` is allowed and should not fail the read.
3. `Quality coverage (%)` should be inspected for unexpected truncation/misalignment.
4. Compare `Legacy map` with `Documented map`.
5. Inspect `Peak positions different (%)` and `Median |peak Δ|`.
6. Compare:
   - `Legacy call is max (%)`
   - `Called-base position call is max (%)`
   - `PeakAmp call is max (%)`

The goal is not to require a specific percentage in alpha.1. We are collecting evidence to decide which representation should become authoritative in the next Stage 1 iteration.

### C. Inspect several known positions

For one clean read and one problematic read:

1. Find 3–5 visually clean bases in the chromatogram.
2. Find 2–3 ambiguous/competing positions if available.
3. In the per-base audit table compare the base, quality value, legacy peak position, called-base peak position and which channel is strongest.
4. Make sure the row numbering matches the chromatogram base numbering.

### D. Export the audit

Download:

- **Run audit CSV** — required for Stage 1 review.
- **Selected-base audit CSV** for at least one clean and one problematic read — strongly preferred.

Send those CSV files back for analysis. They contain numeric evidence only, not the original AB1 binary.

## Cloud check

After local tests are green, regenerate `manifest.json` because Stage 1 adds a new R source file:

```r
source("prepare_connect_cloud.R")
```

Then commit/push. Confirm the online app opens and can process at least one AB1 through QC. The new audit panel must appear online as well.

## Stage 1 alpha.1 is green when

- all four automated test groups pass;
- the normal v2 workflow shows no regression;
- real AB1 files produce an audit without stopping processing;
- exported run evidence is available for review;
- the Connect Cloud build and one online AB1 run succeed.

Do not start Stage 2 yet. The first decision after this gate is whether Stage 1 evidence validates changing the authoritative peak-position/channel model and whether basecaller quality can be safely incorporated into the future read/consensus model.

## Technical basis for the audit

The Stage 1 audit follows the `sangerseqR` class documentation rather than inferring these structures from observed data:

- `traceMatrix`: four normalized trace-signal columns in A,C,G,T order.
- `peakPosMatrix`: within each base-call window, the position of the maximum A/C/G/T peak, in A,C,G,T column order.
- `peakAmpMatrix`: corresponding maximum A/C/G/T peak amplitudes, in A,C,G,T column order.
- `read.abif()`: exposes the raw ABIF data-field list; available fields vary by instrument/basecaller version. Stage 1 therefore treats PCON as optional evidence, not as a required input.

Primary documentation:

- https://github.com/bioc/sangerseqR/blob/devel/man/sangerseq-class.Rd
- https://github.com/bioc/sangerseqR/blob/devel/man/abif-class.Rd

Stage 1 deliberately labels PCON as **basecaller quality** rather than assuming that every instrument/basecaller/version supplies an identically calibrated Phred implementation. The empirical behavior on the laboratory's own AB1 files will be reviewed before quality is used in consensus decisions.

## alpha.1.1 test hotfix

The original alpha.1 synthetic signal fixture used A=100 and strongest alternative=10 for its first base, which is exactly a 10.0 ratio. The test incorrectly asserted `> 10`. alpha.1.1 corrects the boundary assertion to `>= 10`. Application and evidence-calculation behavior are unchanged.
