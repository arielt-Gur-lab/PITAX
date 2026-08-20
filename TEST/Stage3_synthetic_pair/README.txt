PITAX Stage 3 controlled AB1 pair
=================================

Purpose
-------
This fixture tests Upload -> Rename -> Trim/QC -> Stage 3 direction, pairing,
reverse-complement orientation and overlap mechanics when no laboratory pair is
available.

Files
-----
- TESTPAIR001_forward.ab1: a public AB1 trace copied without sequence changes.
- TESTPAIR001_reverse.ab1: a trace-aware reverse-complement derivative.
- assignment_key.csv: ready-to-import Rename key.
- expected_forward_calls.fasta: full PBAS.2 base calls before PITAX trimming.
- BIOPYTHON_LICENSE.rst: license for the public source fixture.

Recommended test
----------------
1. Upload both AB1 files.
2. Choose Paired reads in Project read model.
3. Set Target/Gene to Other in Assay.
4. In Rename, import assignment_key.csv and apply it.
5. Confirm both rows are TESTPAIR001 / TEST_LOCUS with one Forward and one Reverse.
6. Start trimming, review QC, continue to Analysis Sequence and build.
7. Confirm one paired record is produced and the oriented reads overlap.

Important limitation
--------------------
The Reverse file is derived computationally from the Forward trace. Therefore
the pair is deliberately controlled and nearly symmetric. Passing this test
does not validate real laboratory behavior such as different read lengths,
independent noise, allele mixtures, poor overlap or conflicting base calls.
Stage 3 should not be closed without at least one independently sequenced pair.

Provenance
----------
The source trace is Biopython Tests/Abi/3100.ab1, downloaded from the official
Biopython repository at commit-independent master path on 2026-08-20. Original
SHA-256: 78588d824dd967a1e58a9667a8f7ca1687ba38266420bbc5a6562c1893888efd.
The source is distributed under the Biopython License Agreement included here.
