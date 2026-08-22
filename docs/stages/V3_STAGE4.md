# PITAX 3.0 — Stage 4 multi-locus isolate profile

Version: `3.0.0-alpha.9.1`

## Purpose

Stage 4 combines separately completed single-locus PITAX projects at isolate level while preserving every locus as independent evidence. It does not concatenate genes and it does not treat loci as votes.

## Input contract

- Input files are PITAX `.sangerproject` files using project schema 4 or newer.
- Each source project must contain exactly one Gene/Locus.
- Its Stage 3 analysis sequences must be current and review-complete.
- Isolate identity comes only from the explicit `Isolate` field.
- The same Isolate/Locus combination may occur only once across all sources.
- The current session may be included directly; other loci are added as saved projects.

## Stored profile

Project schema 5 stores `pitax-multilocus-profile-v1` with:

- source-project provenance: filename, application version, save time, project schema, locus and MD5 fingerprint where available;
- one evidence row per Isolate/Locus;
- the exact sequence and consensus revision used;
- Stage 3 status and sequence length;
- locus-level taxonomic recommendation, supported rank and confidence;
- best molecular match, accession, Identity and query coverage;
- reference support, locus discrimination, BLAST Request ID (RID) and taxonomy time;
- one conservative profile row per isolate.

## Profile states

| State | Meaning |
| --- | --- |
| `CONCORDANT_SPECIES` | The same species-level call is supported by at least two loci. |
| `CONCORDANT_GENUS` | At least two loci support the same genus, but species remains unresolved or locus-dependent. |
| `SPECIES_CONFLICT` | Species-level calls differ; the conflict is retained. |
| `GENUS_CONFLICT` | Supported genera differ; no combined identification is reported. |
| `PARTIAL_EVIDENCE` | Some required locus-level taxonomic evidence is missing or insufficient. |
| `NO_TAXONOMY` | Sequences exist but no locus has completed taxonomic interpretation. |
| `SINGLE_LOCUS` | Only one locus is present, so the object is not yet a multi-locus profile. |

The algorithm never selects a taxon by counting loci. A 2:1 split across genera remains `GENUS_CONFLICT`.

## Visual interpretation

The Stage 4 workspace starts with a visual review layer while retaining all audit tables:

- an isolate selector that supports any number of isolate profiles;
- project-wide counts for isolates, distinct loci, concordant profiles and profiles requiring attention;
- one combined status, conclusion and next action for the selected isolate;
- one evidence card per locus with taxonomic call, confidence, best match, accession, reference context and locus discrimination;
- an Identity/query-coverage plot for direct comparison of the selected isolate's analyzed loci.

Selection is based on the explicit Isolate field. Evidence from different isolates is never mixed in one visual profile.

## Staleness

When the current session was included, PITAX compares its sequence, consensus revision, taxonomy status, recommendation and RID with the stored Stage 4 snapshot. Any change marks the profile stale and blocks Stage 4 downloads until it is rebuilt.

Imported source projects are immutable snapshots identified by their stored provenance and fingerprint. Re-import and rebuild after changing an external locus project.

## Exports

- isolate profile CSV;
- complete per-locus evidence CSV;
- multi-locus FASTA with one record per Isolate/Locus;
- Checkpoint F ZIP containing source provenance, profiles, evidence, FASTA and the interpretation contract.

## Current alpha boundary

Alpha.9.1 implements the Stage 4 data architecture, conservative concordance/conflict behavior and visual multi-isolate review. It does not yet automate taxon-specific marker recommendations. That next increment requires a provenance-bearing literature layer and curated Fusarium, Trichoderma, Aspergillus, Penicillium, Pleurotus and general-fungi scenarios.

## Automated tests

Run `run_tests.bat` on Windows with R installed. Stage 4 is covered by groups 11–13:

1. multi-locus logic, provenance, duplicate blocking, missing taxonomy, staleness, migration and FASTA;
2. Shiny integration, visual multi-isolate selection, project persistence and export contracts;
3. global DataTables header/body alignment contracts.
