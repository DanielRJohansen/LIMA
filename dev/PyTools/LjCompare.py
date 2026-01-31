#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Record:
    pid_pair: Tuple[int, int]   # sorted
    origin: str
    force_scalar: float
    raw_line: str


LINE_RE = re.compile(
    r"""
    \borigin\s+(?P<origin>[A-Za-z0-9_]+)\s+
    pIds:\s+(?P<p0>-?\d+)\s+(?P<p1>-?\d+)
    """,
    re.VERBOSE,
)

FORCE_RE = re.compile(r"\bforceScalar\s+(?P<force>[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\b")


def _parse_float(s: str) -> float:
    v = float(s)
    if math.isnan(v) or math.isinf(v):
        return v
    return v


def parse_file(path: Path) -> List[Record]:
    records: List[Record] = []
    for i, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        if not line.strip():
            continue

        m1 = LINE_RE.search(line)
        m2 = FORCE_RE.search(line)
        if not (m1 and m2):
            continue

        origin = m1.group("origin")
        p0 = int(m1.group("p0"))
        p1 = int(m1.group("p1"))
        pair = tuple(sorted((p0, p1)))
        force = _parse_float(m2.group("force"))

        records.append(Record(pair, origin, force, f"{path.name}:{i}: {line}"))
    return records


def rel_err(a: float, b: float) -> float:
    # symmetric relative difference; stable when magnitudes differ
    denom = abs(a) + abs(b)
    if denom == 0.0:
        return 0.0
    return abs(a - b) / denom


def fmt_float(v: float, width: int = 14, prec: int = 6) -> str:
    if math.isnan(v):
        return f"{'nan':>{width}}"
    if math.isinf(v):
        return f"{'inf' if v > 0 else '-inf':>{width}}"
    return f"{v:{width}.{prec}f}"


def analyze(records: List[Record]) -> Dict[Tuple[int, int], Dict[str, List[Record]]]:
    by_pair: Dict[Tuple[int, int], Dict[str, List[Record]]] = {}
    for r in records:
        by_pair.setdefault(r.pid_pair, {}).setdefault(r.origin, []).append(r)
    return by_pair


def pick_value(recs: List[Record], agg: str) -> float:
    vals = [r.force_scalar for r in recs if not (math.isnan(r.force_scalar) or math.isinf(r.force_scalar))]
    if not vals:
        # if only nan/inf present, return first to report it explicitly
        return recs[0].force_scalar
    if agg == "mean":
        return sum(vals) / len(vals)
    if agg == "median":
        s = sorted(vals)
        n = len(s)
        mid = n // 2
        return s[mid] if (n % 2) else 0.5 * (s[mid - 1] + s[mid])
    if agg == "maxabs":
        return max(vals, key=lambda x: abs(x))
    raise ValueError(f"Unknown agg: {agg}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare forceScalar for same pIds across origins (PP vs ComComIntra).")
    ap.add_argument("file", type=Path, nargs="?", default=Path("LjCompareData.txt"),
                    help="Path to .txt with kernel printf output (default: LjCompareData.txt)")
    ap.add_argument("--origin-a", default="PP", help="First origin name (default: PP)")
    ap.add_argument("--origin-b", default="ComComIntra", help="Second origin name (default: ComComIntra)")
    ap.add_argument("--min-rel", type=float, default=0.001,
                    help="Minimum symmetric relative difference to report (default: 0.001)")
    ap.add_argument("--agg", choices=["mean", "median", "maxabs"], default="mean",
                    help="If multiple lines exist per (pair, origin), aggregate forceScalar this way (default: mean)")
    ap.add_argument("--show-missing", action="store_true",
                    help="Also list pId pairs that only exist in one of the origins")
    ap.add_argument("--show-lines", action="store_true",
                    help="Include the raw source lines for reported pairs")
    args = ap.parse_args()

    records = parse_file(args.file)
    by_pair = analyze(records)

    a = args.origin_a
    b = args.origin_b

    diffs: List[Tuple[float, Tuple[int, int], float, float, int, int]] = []
    missing_a: List[Tuple[Tuple[int, int], int]] = []
    missing_b: List[Tuple[Tuple[int, int], int]] = []

    for pair, origins in by_pair.items():
        has_a = a in origins
        has_b = b in origins

        if not has_a and not has_b:
            continue

        if has_a and has_b:
            va = pick_value(origins[a], args.agg)
            vb = pick_value(origins[b], args.agg)
            r = rel_err(va, vb)
            if r >= args.min_rel:
                diffs.append((r, pair, va, vb, len(origins[a]), len(origins[b])))
        else:
            if args.show_missing:
                if not has_a and has_b:
                    missing_a.append((pair, len(origins[b])))
                if has_a and not has_b:
                    missing_b.append((pair, len(origins[a])))

    diffs.sort(key=lambda t: t[0], reverse=True)

    # Report
    print(f"File: {args.file}")
    print(f"Origins compared: {a} vs {b}")
    print(f"Aggregation: {args.agg}")
    print(f"Report threshold: rel >= {args.min_rel:.3f}  (rel = |a-b|/(|a|+|b|))")
    print()

    header = (
        f"{'p0':>3} {'p1':>3} |"
        f" {('F_'+a):>14} {('F_'+b):>14} |"
        f" {'rel':>8} |"
        f" {'nA':>3} {'nB':>3}"
    )
    print(header)
    print("-" * len(header))

    for r, pair, va, vb, nA, nB in diffs:
        p0, p1 = pair
        print(
            f"{p0:3d} {p1:3d} |"
            f" {fmt_float(va)} {fmt_float(vb)} |"
            f" {r:8.4f} |"
            f" {nA:3d} {nB:3d}"
        )
        if args.show_lines:
            for rec in by_pair[pair].get(a, []):
                print(f"  A: {rec.raw_line}")
            for rec in by_pair[pair].get(b, []):
                print(f"  B: {rec.raw_line}")
            print()

    if not diffs:
        print("(no pairs exceeded threshold)")

    if args.show_missing:
        if missing_a:
            print()
            print(f"Pairs missing in {a} (present in {b}):")
            for (p0, p1), n in sorted(missing_a):
                print(f"  {p0:3d} {p1:3d}  (n={n})")
        if missing_b:
            print()
            print(f"Pairs missing in {b} (present in {a}):")
            for (p0, p1), n in sorted(missing_b):
                print(f"  {p0:3d} {p1:3d}  (n={n})")


if __name__ == "__main__":
    main()
