#!/usr/bin/env python3
"""Create a trace-aware reverse-complement AB1 fixture from one AB1 source.

The output is for deterministic pipeline testing. It is not an independently
sequenced biological Forward/Reverse pair.
"""

from __future__ import annotations

import argparse
import shutil
import struct
from dataclasses import dataclass
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


@dataclass(frozen=True)
class Entry:
    name: str
    number: int
    element_type: int
    element_size: int
    count: int
    data_size: int
    data_offset: int
    entry_offset: int


class AbifFile:
    def __init__(self, payload: bytes):
        if payload[:4] != b"ABIF":
            raise ValueError("Input is not an ABIF file")
        self.data = bytearray(payload)
        root = payload[6:34]
        self.directory_entry_size = struct.unpack(">H", root[10:12])[0]
        self.directory_count = struct.unpack(">I", root[12:16])[0]
        self.directory_offset = struct.unpack(">I", root[20:24])[0]
        self.entries: dict[tuple[str, int], Entry] = {}
        for index in range(self.directory_count):
            offset = self.directory_offset + index * self.directory_entry_size
            raw = payload[offset : offset + self.directory_entry_size]
            name = raw[:4].decode("latin1")
            number = struct.unpack(">I", raw[4:8])[0]
            element_type, element_size, count, data_size, data_offset = struct.unpack(">HHIII", raw[8:24])
            self.entries[(name, number)] = Entry(
                name, number, element_type, element_size, count, data_size, data_offset, offset
            )

    def entry(self, name: str, number: int) -> Entry:
        try:
            return self.entries[(name, number)]
        except KeyError as exc:
            raise ValueError(f"Required ABIF tag {name}.{number} is missing") from exc

    def read(self, name: str, number: int) -> bytes:
        entry = self.entry(name, number)
        if entry.data_size <= 4:
            return bytes(self.data[entry.entry_offset + 20 : entry.entry_offset + 20 + entry.data_size])
        return bytes(self.data[entry.data_offset : entry.data_offset + entry.data_size])

    def write(self, name: str, number: int, value: bytes) -> None:
        entry = self.entry(name, number)
        if len(value) != entry.data_size:
            raise ValueError(f"Replacement for {name}.{number} must remain {entry.data_size} bytes")
        offset = entry.entry_offset + 20 if entry.data_size <= 4 else entry.data_offset
        self.data[offset : offset + entry.data_size] = value

    def write_if_present(self, name: str, number: int, value: bytes) -> None:
        if (name, number) in self.entries:
            self.write(name, number, value)


def reverse_elements(raw: bytes, element_size: int) -> bytes:
    if len(raw) % element_size:
        raise ValueError("ABIF element block is not aligned to its declared element size")
    return b"".join(raw[i : i + element_size] for i in range(0, len(raw), element_size))[::-1] if element_size == 1 else b"".join(
        raw[i : i + element_size] for i in range(len(raw) - element_size, -1, -element_size)
    )


def reverse_complement_text(raw: bytes) -> bytes:
    return raw.decode("ascii").translate(COMPLEMENT)[::-1].encode("ascii")


def add_controlled_conflict(abif: AbifFile, forward_position: int) -> None:
    calls = bytearray(abif.read("PBAS", 2))
    reverse_index = len(calls) - forward_position
    if not 0 <= reverse_index < len(calls):
        raise ValueError("Conflict position is outside the called-base sequence")

    raw_to_oriented = str.maketrans("ACGT", "TGCA")
    oriented_to_alternative = {"A": "C", "C": "A", "G": "T", "T": "G"}
    old_raw = chr(calls[reverse_index])
    old_oriented = old_raw.translate(raw_to_oriented)
    new_oriented = oriented_to_alternative[old_oriented]
    new_raw = new_oriented.translate(raw_to_oriented)

    for tag in (("PBAS", 1), ("PBAS", 2), ("P2BA", 1)):
        if tag not in abif.entries:
            continue
        value = bytearray(abif.read(*tag))
        if len(value) == len(calls):
            value[reverse_index] = ord(new_raw)
            abif.write(*tag, bytes(value))

    positions_raw = abif.read("PLOC", 2)
    positions = struct.unpack(f">{len(positions_raw) // 2}h", positions_raw)
    peak = positions[reverse_index]
    trace_order = abif.read("FWO_", 1).decode("ascii")
    old_tag = ("DATA", 9 + trace_order.index(old_raw))
    new_tag = ("DATA", 9 + trace_order.index(new_raw))
    old_trace = list(struct.unpack(f">{abif.entry(*old_tag).count}h", abif.read(*old_tag)))
    new_trace = list(struct.unpack(f">{abif.entry(*new_tag).count}h", abif.read(*new_tag)))
    for sample in range(max(0, peak - 4), min(len(old_trace), peak + 5)):
        old_trace[sample], new_trace[sample] = new_trace[sample], old_trace[sample]
    abif.write(*old_tag, struct.pack(f">{len(old_trace)}h", *old_trace))
    abif.write(*new_tag, struct.pack(f">{len(new_trace)}h", *new_trace))


def build_pair(source: Path, forward: Path, reverse: Path, conflict_forward_position: int | None = None) -> None:
    source_bytes = source.read_bytes()
    forward.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, forward)
    abif = AbifFile(source_bytes)

    for tag in (("PBAS", 1), ("PBAS", 2), ("P2BA", 1)):
        if tag in abif.entries:
            abif.write(*tag, reverse_complement_text(abif.read(*tag)))

    for tag in (("PCON", 1), ("PCON", 2)):
        if tag in abif.entries:
            abif.write(*tag, abif.read(*tag)[::-1])

    trace_order = abif.read("FWO_", 1).decode("ascii")
    trace_entries = {base: ("DATA", 9 + index) for index, base in enumerate(trace_order)}
    trace_length = abif.entry("DATA", 9).count
    original_traces = {base: abif.read(*tag) for base, tag in trace_entries.items()}
    base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for base, tag in trace_entries.items():
        abif.write(*tag, reverse_elements(original_traces[base_complement[base]], 2))

    for tag in (("PLOC", 1), ("PLOC", 2)):
        if tag not in abif.entries:
            continue
        raw = abif.read(*tag)
        positions = list(struct.unpack(f">{len(raw) // 2}h", raw))
        transformed = [trace_length - 1 - position for position in reversed(positions)]
        if transformed != sorted(transformed):
            raise ValueError(f"Transformed {tag[0]}.{tag[1]} positions are not ascending")
        abif.write(*tag, struct.pack(f">{len(transformed)}h", *transformed))

    for tag in (("P1AM", 1), ("P2AM", 1)):
        if tag in abif.entries:
            abif.write(*tag, reverse_elements(abif.read(*tag), 2))

    if conflict_forward_position is not None:
        add_controlled_conflict(abif, conflict_forward_position)

    reverse.write_bytes(abif.data)

    check = AbifFile(reverse.read_bytes())
    expected = reverse_complement_text(AbifFile(source_bytes).read("PBAS", 2))
    mismatch_count = sum(a != b for a, b in zip(check.read("PBAS", 2), expected))
    expected_mismatches = 1 if conflict_forward_position is not None else 0
    if mismatch_count != expected_mismatches:
        raise ValueError("Reverse AB1 PBAS.2 verification failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("forward", type=Path)
    parser.add_argument("reverse", type=Path)
    parser.add_argument("--conflict-forward-position", type=int)
    args = parser.parse_args()
    build_pair(args.source, args.forward, args.reverse, args.conflict_forward_position)


if __name__ == "__main__":
    main()
