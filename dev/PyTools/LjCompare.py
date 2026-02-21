#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
from enum import Enum

class Origin(Enum):
    ComCom_ = 0
    PP = 1


Float3 = Tuple[float, float, float]


# -----------------------------
# Parsing
# -----------------------------


@dataclass(frozen=True)
class Record:
    pids: Tuple[int, int]   # sorted
    ljDiff: Float3
    dist: float
    sigma: float
    eps: float
    forceMag: float
    origin: Origin

_LINE_RE = re.compile(
    r"""
    LJ:\s+diff:\s+
        (?P<dx>[-+0-9.eE]+)\s+
        (?P<dy>[-+0-9.eE]+)\s+
        (?P<dz>[-+0-9.eE]+)
    \s+dist\s+(?P<dist>[-+0-9.eE]+)
    \s+sigma:\s+(?P<sigma>[-+0-9.eE]+)
    \s+eps:\s+(?P<eps>[-+0-9.eE]+)
    .*?
    forceMagnitude\s+(?P<force>[-+0-9.eE]+)
    \s+origin\s+(?P<origin>\w+)
    \s+pIds:\s+(?P<p0>-?\d+)\s+(?P<p1>-?\d+)
    """,
    re.VERBOSE,
)

def _ParseOrigin(originToken: str) -> Origin:
    if originToken.startswith("ComCom"):
        return Origin.ComCom_
    if originToken == "PP":
        return Origin.PP
    raise ValueError(f"Unknown origin: {originToken}")


def ParseRow(row: str) -> Record:
    match = _LINE_RE.search(row)
    if not match:
        raise ValueError(f"Failed to parse row: {row}")

    p0 = int(match["p0"])
    p1 = int(match["p1"])

    return Record(
        pids=tuple(sorted((p0, p1))),
        ljDiff=(
            float(match["dx"]),
            float(match["dy"]),
            float(match["dz"]),
        ),
        dist=float(match["dist"]),
        sigma=float(match["sigma"]),
        eps=float(match["eps"]),
        forceMag=float(match["force"]),
        origin=_ParseOrigin(match["origin"]),
    )

def ParseFile(path: str) -> Tuple[Dict[Tuple[int, int], Record],
                                  Dict[Tuple[int, int], Record]]:
    ppRecords: Dict[Tuple[int, int], Record] = {}
    comcomRecords: Dict[Tuple[int, int], Record] = {}

    with open(path, "r", encoding="utf-8") as file:
        for row in file:
            rec = ParseRow(row)
            targetMap = comcomRecords if rec.origin == Origin.ComCom_ else ppRecords
            if rec.pids not in targetMap:
                targetMap[rec.pids] = rec

    return ppRecords, comcomRecords



# -----------------------------
# Analysis
# -----------------------------

@dataclass(frozen=True)
class ForceDiff:
    pids: Tuple[int, int]
    ppForce: float
    comcomForce: float
    absDiff: float


@dataclass(frozen=True)
class MissingEntry:
    origin: str
    pids: Tuple[int, int]
    forceMag: float


@dataclass(frozen=True)
class ReportAnalysis:
    diffs: List[ForceDiff]
    missing: List[MissingEntry]
    maxDiff: float
    avgDiff: float
    maxMissing: float
    avgMissing: float



# Generate and print a report. First report all record sets with the same key, where the difference in forceMagnitude is largest, largest first. Report no more than 20.
# Then report all the cases where an entry is missing in either map, once again prioritise the records with largest forceMagnitudes
def AnalyseRecords(ppRecords: Dict[Tuple[int, int], Record],
                   comcomRecords: Dict[Tuple[int, int], Record]) -> ReportAnalysis:
    diffs: List[ForceDiff] = []

    relThreshold = 0.01 # Dont care bout diff <1%
    maxAbsDiff = 0

    for key in ppRecords.keys() & comcomRecords.keys():
        pp = ppRecords[key]
        cc = comcomRecords[key]
        absDiff = abs(pp.forceMag - cc.forceMag)
        relDiff = absDiff / cc.forceMag
        maxAbsDiff = max(maxAbsDiff, absDiff)
        if (relDiff < relThreshold):
            continue
        diffs.append(
            ForceDiff(
                pids=key,
                ppForce=pp.forceMag,
                comcomForce=cc.forceMag,
                absDiff=diff,
            )
        )

    diffs.sort(key=lambda d: d.absDiff, reverse=True)

    maxDiff = diffs[0].absDiff if diffs else 0.0
    avgDiff = sum(d.absDiff for d in diffs) / len(diffs) if diffs else 0.0

    missing: List[MissingEntry] = []
    absThreshold = 50 # MN or whatever forcesystem we use... Anything such a small force is so far away, the particle was probably just sorted off purposefully in NLists

    for key in ppRecords.keys() - comcomRecords.keys():
        rec = ppRecords[key]
        if abs(rec.forceMag) > absThreshold:
            missing.append(MissingEntry("PP", key, rec.forceMag))

    for key in comcomRecords.keys() - ppRecords.keys():
        rec = comcomRecords[key]
        if abs(rec.forceMag) > absThreshold:
            missing.append(MissingEntry("ComCom", key, rec.forceMag))

    missing.sort(key=lambda m: abs(m.forceMag), reverse=True)

    maxMissing = abs(missing[0].forceMag) if missing else 0.0
    avgMissing = (
        sum(abs(m.forceMag) for m in missing) / len(missing)
        if missing else 0.0
    )

    return ReportAnalysis(
        diffs=diffs,
        missing=missing,
        maxDiff=maxDiff,
        avgDiff=avgDiff,
        maxMissing=maxMissing,
        avgMissing=avgMissing,
    )


def _FmtPids(pids: Tuple[int, int]) -> str:
    return f"{pids[0]:6d},{pids[1]:6d}"


def PrintReport(analysis: ReportAnalysis, maxRows: int = 20):
    def _Hdr(title: str):
        print("\n" + title)
        print("=" * len(title))

    # ---- Diff section
    _Hdr("Force magnitude mismatches (PP vs ComCom)")
    print(f"Total comparable pairs : {len(analysis.diffs)}")
    print(f"Largest |ΔF|           : {analysis.maxDiff:12.6f}")
    print(f"Average |ΔF|           : {analysis.avgDiff:12.6f}")

    print(
        "\n"
        "   pId0 ,   pId1 | "
        "PP force       | "
        "ComCom force   | "
        "|ΔF|\n"
        "----------------+----------------+----------------+----------------"
    )

    for d in analysis.diffs[:maxRows]:
        print(
            f"{_FmtPids(d.pids)} | "
            f"{d.ppForce:14.6f} | "
            f"{d.comcomForce:14.6f} | "
            f"{d.absDiff:14.6f}"
        )

    # ---- Missing section
    _Hdr("Missing entries between PP and ComCom")
    print(f"Total missing pairs    : {len(analysis.missing)}")
    print(f"Largest |F| missing    : {analysis.maxMissing:12.6f}")
    print(f"Average |F| missing    : {analysis.avgMissing:12.6f}")

    print(
        "\n"
        "Origin  |    pId0 ,   pId1 | forceMagnitude\n"
        "--------+-----------------+----------------"
    )

    for m in analysis.missing[:maxRows]:
        print(
            f"{m.origin:6s} | "
            f"{_FmtPids(m.pids)} | "
            f"{m.forceMag:14.6f}"
        )






def main():
    ppRecords, comcomRecords = ParseFile("LjCompareData.txt")
    analysis = AnalyseRecords(ppRecords, comcomRecords)
    PrintReport(analysis)



if __name__ == "__main__":
    main()
