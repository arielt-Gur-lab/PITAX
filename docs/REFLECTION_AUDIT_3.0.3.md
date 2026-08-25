# PITAX 3.0.3 — Reflection audit findings

**Date:** 2026-08-25  
**Scope:** Static review + full automated suite + engineering probes (no product code changes)  
**App version:** 3.0.3 · **Project schema:** 6 · **R:** 4.5.2  

---

## Executive verdict

**Go for Alpha 10 multi-assay editor — with conditions.**

Domain logic for Stages 2–4 and taxonomic *decision* rules is in good shape and is the best-tested layer. The AB1 → trim → consensus path works on controlled fixtures.

### Follow-up status (3.0.4)

Implemented after this audit: logo MD5 guard removed; BUG-01 / BUG-03 / BUG-04 / BUG-05 / BUG-06 fixed; suite expanded to 16 groups and is green.

### Top three risks that were blocking Alpha 10 (now addressed in 3.0.4)

1. **BUG-01** — Curating one read cleared *all* consensus while leaving other isolates' BLAST jobs READY.  
2. **BUG-06** — Assignment signature omitted `Assay_ID`.  
3. **BUG-03** — Schema 5→6 mapped unknown/`Other` loci to **ITS** silently.

---

## Baseline (Stage 0)

| Check | Result |
|-------|--------|
| Preflight parse (45 R files) | Pass |
| Groups 1–13 | Pass (~24 s through group 13) |
| Group 14 structure / logo MD5 | **FAIL** |
| Group 15 Alpha 10 assay/schema | Pass (run separately) |
| Doc drift | `docs/stages/V3_STAGE3.md` still says **thirteen** runner groups; runner is **15** |

**ENV-01 (resolved in 3.0.4):** Logo MD5 byte guard removed by product decision; artwork may change without failing the structure contract. `www/logo.png` remains a required file.

---

## Coverage map (Stage 1)

| Area | LOC (approx) | Automated coverage | Gap |
|------|--------------|--------------------|-----|
| `core_sanger.R` | 1431 | Structure markers only | No unit tests for trim/channel/metrics |
| `sequence_tools.R` | 998 | Peak flags + curation slice | Primer map / plots largely untested |
| `taxonomy_tools.R` | 729 | Strong decision-engine smoke | I/O enrichment less covered |
| `export_tools.R` | 638 | Light / name contracts | BLAST XML/CSV parse untested |
| `stage3_consensus.R` | 627 | Strong algorithm smoke | E2E trim→consensus was missing (probed here) |
| `90_blast.R` | 573 | Empty-selector contract | Submit/retrieve/stale behavior untested |
| `100_taxonomy.R` | 695 | Via taxonomy unit only | Server RID binding untested |
| `50_evidence_review.R` | 566 | Domain curation unit | Shiny paths untested |
| `stage4_multilocus.R` | 376 | Strong | — |
| `assay_profiles.R` | 202 | Alpha 10 smoke | Silent Other→ITS edge |
| Stage 2/4 “app contracts” | — | Source-string guards | Not behavioral correctness |

**Rule confirmed:** Stage 2/4 app contracts are refactor nets, not scientific proof.

---

## Correctness findings (Stage 2)

### BUG-01 — Asymmetric invalidation after curation (P0)

**Where:** `R/server/stages/10_project.R` → `commit_curated_result`  
**What:** On sequence change, `rv$consensus_set` is replaced with an empty set for the **entire project**, but `invalidate_downstream_for_sample` runs only for affected consensus IDs + the curated sample.  
**Probe E:** with isolates ISO1 + ISO2, editing `readA_F` clears all consensus and invalidates ISO1; **ISO2 BLAST remains notionally READY** while Stage 3 is empty.  
**Why it matters:** Export/BLAST UI can show leftover READY evidence for isolates whose consensus was wiped; operators may trust stale NCBI results.  
**Recommendation:** Either (a) invalidate BLAST/taxonomy for **all** consensus IDs when clearing the set, or (b) clear only consensus records that include the edited read.

### BUG-02 — Logo contract failure in working tree (P0 env)

See ENV-01 above. Not a runtime science bug; blocks the declared “all tests passed” gate.

### BUG-03 — Silent `Other` / unknown → `ITS` on schema 5→6 (P1)

**Where:** `R/domain/assay/assay_profiles.R` → `pitax_normalize_locus_id(..., fallback)`, `assay_profile_from_legacy_settings`, `assay_migrate_schema5_state`  
**Probe A/B:** `Other` and `SomethingWeird` both become `ITS`; migration force-links every read to that assay.  
**Why it matters:** Legacy free-text loci remount as ITS while BLAST/taxonomy evidence is preserved → locus label can disagree with historical interpretation.  
**Recommendation:** Fail migration (or require explicit operator remap) when legacy target is outside the controlled vocabulary; never default unknown → ITS silently.

### BUG-04 — `latest_blast_rid_for_sample` ignores status (P1)

**Where:** `R/server/stages/100_taxonomy.R` L2–6  
**What:** Returns the **last row** RID for the sample with no READY/STALE filter.  
**Probe D:** job order READY, READY, STALE → taxonomy would pick `RID_STALE`.  
Hit rows for that RID may already have been deleted by invalidation, yielding empty taxonomy rather than the previous READY set — confusing, and dangerous if STALE hits were ever left in place.  
**Recommendation:** Prefer latest job with `status == "READY"` (else explicit empty / message).

### BUG-05 — `retrieve_blast_job` trusts STALE flag only (P1)

**Where:** `R/server/stages/90_blast.R` → `retrieve_blast_job`  
**What:** Blocks retrieve when status is STALE; does **not** re-compare stored `consensus_revision` to the current analysis revision.  
**Why it matters:** If a job is incorrectly left non-STALE (see BUG-01), retrieve can still apply.  
**Recommendation:** Re-check revision (and/or sequence signature) on retrieve, not only status.

### BUG-06 — Assignment signature omits `Assay_ID` (P1 / Alpha 10)

**Where:** `R/server/stages/20_upload.R` → `assignment_state_signature`  
Columns: Source_ID, Isolate, Locus, Direction, Final_Name — **no Assay_ID**.  
**Probe C:** Two TEF1 assays with different primers produce **identical** signatures and Final_Names (`ISO1_TEF1_F`).  
**Why it matters:** Once multi-assay lands, swapping Assay_ID without changing locus/direction will **not** clear consensus/BLAST.  
**Recommendation:** Include `Assay_ID` (and maybe primer identity) in the signature before shipping the multi-assay editor.

### BUG-07 — Stage 4 extract errors look like session drift (P2)

**Where:** `R/domain/multilocus/stage4_multilocus.R` → `stage4_current_project_matches`  
**Probe F:** Stored “Current session” rows + broken project → `FALSE` (same as real staleness).  
**Recommendation:** Surface extract failures as explicit errors, not only “rebuild required”.

### OBS-01 — Closed / recurring classes (still guarded)

| Class | Status |
|-------|--------|
| ASCII / encoding runtime contracts | Covered by group 14 (when logo OK) |
| DataTables scroll + autoWidth | Group 13 pass |
| Plotly empty-trace warnings | Still appear in Stage 3 Shiny startup test (lesson §11.2 not fully extinguished on empty plots) |
| Taxonomy decision conservatism | Strong unit coverage; do not loosen without new tests |

### OBS-02 — `tryCatch` → NULL / empty (P2 pattern)

Present in load architecture, BLAST CSV fallback, taxonomy enrich, AB1 audit isolation. Often intentional (keep pipeline alive). Risk is **silent empty evidence**. Prefer explicit status strings where the user must act.

---

## Efficiency findings (Stage 3)

| ID | Severity | Issue | Recommendation |
|----|----------|-------|----------------|
| PERF-01 | P2 now / Stage 5 | `sync_summary_from_results` rebuilds **all** reads after one curation (`00_state.R`) | Patch only the edited sample |
| PERF-02 | P2 | Chromatogram ships full `scattergl` + per-base text labels (`sequence_tools.R`) | Overview downsample; full detail on zoom window |
| PERF-03 | P2 | Peak flags recomputed for QC UI and Team Summary | Cache on result object |
| PERF-04 | P3 | Consensus traceback uses repeated `c()` (`stage3_consensus.R` ~L308) | Preallocate / rev-at-end; low risk at ~1 kb Sanger |
| PERF-05 | P2 | `www/pitax.js` schedules DT adjust on every `shiny:value` | Prefer explicit `adjustDataTables` after known renders |
| PERF-06 | P2 | Assignment editor `renderUI` builds N×M live inputs (`60_assignment.R`) | DT / fewer inputs for plate-scale uploads |
| PERF-07 | Stage 5 | NCBI throttle is session-local (`90_blast.R`); `Sys.sleep` blocks session | Global scheduler (already on Stage 5 roadmap) |
| PERF-08 | P3 | Duplicate `make_fasta` / `wrap_sequence` in `core_sanger.R` and `export_tools.R` (export wins after source order) | Single helper module |
| PERF-09 | Stage 5 | Bootstrap may `install.packages` at startup | Deploy with frozen package set only |

**Lab-scale judgment:** For typical fungal plates (≲48 reads), scientific cost (trim, overlap DP, Plotly once) is acceptable. PERF-01/02/05/06 become noticeable first as project size grows.

---

## Manual / engineering scenarios (Stage 4)

| Scenario | Result |
|----------|--------|
| Schema Other→ITS | Confirmed silent remap + force-link |
| Assignment Assay_ID swap | Signature unchanged |
| Multi-isolate invalidation logic | Confirmed asymmetry |
| Taxonomy last-RID | Picks STALE when last |
| Stage4 broken extract | Returns FALSE (false “stale”) |
| **Synthetic AB1 → trim → consensus** | **Pass:** lengths 709/672 trimmed; consensus READY, 709 bp, 0 review rows |
| **Conflict AB1 → trim → consensus** | **Pass:** status REVIEW_REQUIRED; IUPAC **M** present; Needs_Review=1. Raw conflict @400 shifts to consensus ~377 after trim_start=24 — expected |

**Note:** Controlled fixtures prove the engineering pipeline. Independently sequenced real paired AB1 truth sets remain deferred (README / roadmap).

---

## Missing tests to add *after* reflection (not now)

1. Behavioral smoke: fixture AB1 → `trim_one_ab1` → `stage3_build_consensus_set` (clean + conflict).  
2. Unit: `commit_curated_result` invalidation breadth with two isolates.  
3. Unit: `latest_blast_rid_for_sample` READY-only selection.  
4. Unit: schema5 target outside vocabulary must fail or require remap.  
5. Unit: assignment signature includes Assay_ID.  
6. Smoke: BLAST XML/CSV parse + accession uniqueness on a stored fixture response.  
7. Refresh stage docs (13 → 15 groups).

---

## Go / no-go for Alpha 10

| Question | Answer |
|----------|--------|
| Is the scientific core stable enough to continue? | **Yes**, with the conditions below |
| Proceed with multi-assay editor? | **Yes after** addressing BUG-01 and BUG-06 (Assay_ID in signature); treat BUG-03 as a migration gate |
| Blocked by suite greenness? | **Yes locally** until ENV-01 (logo) is resolved |
| Defer to Stage 5? | PERF-07 NCBI scheduler, project-size memory, Connect concurrency |

**Suggested fix order (when you approve implementation):**  
1. Restore/approve logo (ENV-01) → green suite  
2. BUG-01 invalidation symmetry  
3. BUG-06 Assay_ID signature (prerequisite for multi-assay)  
4. BUG-03 migration hard-fail for unknown loci  
5. BUG-04/05 taxonomy/BLAST RID hygiene  
6. Optional PERF-01/05 quick wins  

---

## Appendix — probe artifacts

Temporary probe scripts used during this audit (safe to delete):

- `tests/tmp_reflection_probes.R`
- `tests/tmp_e2e_probe.R`
- `tests/tmp_e2e_detail.R`
- `tests/tmp_reflection_out.txt` / `tests/tmp_e2e_out.txt` (if present)
