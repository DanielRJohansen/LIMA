#!/usr/bin/env python3
"""Capture solo and batched T4 kernel timings, then write an HTML comparison.

Run from any directory:
    python dev/profile_t4_batch.py
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIR = ROOT / "build" / "x64_Release"
DEFAULT_OUTPUT_DIR = ROOT / "dev" / "t4-batch-profile"
KERNEL_NAMES = (
    "NbNonlocalKernel",
    "BondgroupsKernel",
    "SuperclusterIntegrateKernel",
    "ReserveInteractions",
    "BuildTasks",
    "PclusterBondgroupsGather",
    "BuildNointeractionMatricesKernel",
    "ClusteringKernel",
    "SortReserveInteractionsOutput",
    "ComputeMeanposAndRadiiForEachPclusterInEachSuperclusterKernel",
    "GetPclusterPositions",
    "ApplyBoundaryCondition",
    "ClusteringPretransferKernel",
    "CompressSuperclusters",
    "SortPClusterIndicesInBlocks",
    "SetSuperclusterSimulationId",
    "ComputeKineticEnergyKernel",
)


def Run(command: list[str], cwd: Path | None = None) -> None:
    print("+", subprocess.list2cmdline(command))
    subprocess.run(command, cwd=cwd, check=True)


#def BuildProfileExecutable(build_dir: Path) -> None:
#    build_command = ["cmake", "--build", str(build_dir), "--target", "limaprofile", "--config", "Release", "-j", "8"]
#    if os.name != "nt":
#        Run(build_command, ROOT)
#        return
#    dev_command = Path(os.environ.get("ProgramFiles", r"C:\Program Files")) / (
#        "Microsoft Visual Studio/2022/Community/Common7/Tools/VsDevCmd.bat")
#    if not dev_command.is_file():
#        Run(build_command, ROOT)
#        return
#    command_text = f'call "{dev_command}" -arch=x64 && {subprocess.list2cmdline(build_command)}'
#    Run(["cmd.exe", "/d", "/s", "/c", command_text], ROOT)


def FindNsys(explicit: Path | None) -> Path:
    if explicit:
        return explicit
    if path := shutil.which("nsys"):
        return Path(path)
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    candidates = sorted((program_files / "NVIDIA Corporation").glob(
        "Nsight Systems 2024.5.1/target-windows-x64/nsys.exe"))
    candidates.extend(sorted((program_files / "NVIDIA Corporation").glob(
        "Nsight Systems */target-windows-x64/nsys.exe")))
    if candidates:
        return candidates[0]
    raise FileNotFoundError("Nsight Systems was not found. Pass --nsys <path-to-nsys.exe>.")


def FindNcu(explicit: Path | None) -> Path:
    if explicit:
        return explicit
    if path := shutil.which("ncu"):
        return Path(path)
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    candidates = sorted((program_files / "NVIDIA Corporation").glob(
        "Nsight Compute */target/windows-desktop-win7-x64/ncu.exe"))
    if candidates:
        return candidates[-1]
    raise FileNotFoundError("Nsight Compute was not found. Pass --ncu <path-to-ncu.exe>.")


def KernelLabel(name: str) -> str | None:
    if "DeviceScanInitKernel" in name:
        return "DeviceScanInitKernel"
    if "DeviceScanKernel" in name:
        return "DeviceScanKernel"
    if "DeviceReduceSingleTileKernel" in name:
        return "DeviceReduceSingleTileKernel"
    if "DeviceReduceKernel" in name:
        return "DeviceReduceKernel"
    if "transform_kernel" in name:
        return "transform_kernel"
    for label in sorted(KERNEL_NAMES, key=len, reverse=True):
        if label in name:
            return label
    return None


def ReadKernelSummary(path: Path) -> dict[str, tuple[int, int]]:
    totals: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    with path.open(newline="", encoding="utf-8-sig") as csv_file:
        for row in csv.DictReader(csv_file):
            label = KernelLabel(row["Name"])
            if not label:
                continue
            totals[label][0] += int(row["Instances"])
            totals[label][1] += int(row["Total Time (ns)"])
    return {name: (calls, total_ns) for name, (calls, total_ns) in totals.items()}


def FormatMs(nanoseconds: int) -> str:
    return f"{nanoseconds / 1_000_000:.3f}"


def FormatUs(nanoseconds: int, calls: int) -> str:
    return f"{nanoseconds / calls / 1_000:.2f}" if calls else "—"


def ScalingSignal(per_system_ratio: float | None) -> tuple[str, str, str]:
    if per_system_ratio is None:
        return "unknown", "No comparison", "—"
    if per_system_ratio <= 0.85:
        return "strong", "Strong scaling", "Keep"
    if per_system_ratio <= 1.0:
        return "neutral", "Near linear", "Watch"
    if per_system_ratio <= 1.25:
        return "watch", "Scaling loss", "Watch"
    return "bottleneck", "Scaling regression", "Investigate"


def WriteReport(path: Path, solo: dict[str, tuple[int, int]], batch: dict[str, tuple[int, int]], steps: int, batch_size: int, batch_name: str) -> None:
    rows = []
    for name in sorted(set(solo) | set(batch), key=lambda key: solo.get(key, (0, 0))[1], reverse=True):
        solo_calls, solo_ns = solo.get(name, (0, 0))
        batch_calls, batch_ns = batch.get(name, (0, 0))
        total_ratio = batch_ns / solo_ns if solo_ns else None
        per_system_ratio = total_ratio / batch_size if total_ratio else None
        ratio = f"{total_ratio:.2f}×" if total_ratio else "—"
        per_system = f"{per_system_ratio:.2f}×" if per_system_ratio else "—"
        signal_class, signal_title, signal_text = ScalingSignal(per_system_ratio)
        bar_width = min((per_system_ratio or 0) / 5 * 100, 100)
        rows.append(
            f"<tr class=\"{signal_class}\">"
            f"<td>{html.escape(name)}</td><td>{solo_calls:,}</td><td>{FormatMs(solo_ns)}</td><td>{FormatUs(solo_ns, solo_calls)}</td>"
            f"<td>{batch_calls:,}</td><td>{FormatMs(batch_ns)}</td><td>{FormatUs(batch_ns, batch_calls)}</td>"
            f"<td>{ratio}</td><td><div class=\"ratio\"><span>{per_system}</span><div class=\"ratio-track\"><div class=\"ratio-bar {signal_class}\" style=\"width:{bar_width:.1f}%\"></div></div></div></td>"
            f"<td><span class=\"signal {signal_class}\" title=\"{signal_title}\">{signal_text}</span></td></tr>"
        )

    integration_solo = solo.get("SuperclusterIntegrateKernel", (0, 0))[1]
    integration_batch = batch.get("SuperclusterIntegrateKernel", (0, 0))[1]
    integration_ratio = integration_batch / integration_solo / batch_size if integration_solo else 0
    document = f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>T4 Batch Kernel Comparison</title>
<style>
body{{font:15px/1.5 system-ui,sans-serif;max-width:1220px;margin:36px auto;padding:0 20px;color:#1f2933}}h1{{margin-bottom:4px}}p{{color:#52606d}}.note{{padding:14px 18px;background:#fff5e5;border-left:4px solid #d97706;margin:22px 0}}.legend{{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0;color:#52606d;font-size:13px}}.signal{{display:inline-block;padding:2px 7px;border-radius:999px;font-size:12px;font-weight:600}}.signal.strong{{background:#dcfce7;color:#166534}}.signal.neutral{{background:#e0f2fe;color:#075985}}.signal.watch{{background:#fef3c7;color:#92400e}}.signal.bottleneck{{background:#fee2e2;color:#991b1b}}.table-wrap{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:9px 10px;border-bottom:1px solid #d9e2ec;text-align:right;white-space:nowrap}}th{{background:#f0f4f8;position:sticky;top:0}}th:first-child,td:first-child{{text-align:left}}tbody tr.bottleneck{{background:#fff7f7}}tbody tr.watch{{background:#fffdf6}}.ratio{{min-width:112px}}.ratio>span{{display:block;font-variant-numeric:tabular-nums;font-weight:600}}.ratio-track{{height:5px;background:#e5e7eb;border-radius:4px;margin-top:4px;overflow:hidden}}.ratio-bar{{height:100%;border-radius:4px}}.ratio-bar.strong{{background:#22c55e}}.ratio-bar.neutral{{background:#38bdf8}}.ratio-bar.watch{{background:#f59e0b}}.ratio-bar.bottleneck{{background:#ef4444}}a{{color:#0969da;margin-right:16px}}
</style><body>
<h1>T4 GPU kernel comparison</h1>
<p>Nsight Systems capture · RTX 4090 · {steps:,} timesteps · solo versus batch of {batch_size} · logging and neighbor-list intervals set by the command line.</p>
<div class="note">Integration per-system ratio: <b>{integration_ratio:.2f}×</b>. Values below 1 improve aggregate throughput. Kernel totals may overlap across CUDA streams, so they are not wall-clock time.</div>
<div class="legend"><span class="signal strong">Keep ≤0.85×</span><span class="signal neutral">Watch 0.85–1.00×</span><span class="signal watch">Watch 1.00–1.25×</span><span class="signal bottleneck">Investigate &gt;1.25×</span></div>
<p>Rows are sorted by solo GPU time. The per-system bar uses 5× as its full scale.</p>
<div class="table-wrap"><table><thead><tr><th>Kernel</th><th>Solo calls</th><th>Solo ms</th><th>Solo μs/call</th><th>Batch calls</th><th>Batch ms</th><th>Batch μs/call</th><th>Total ratio</th><th>Per-system ratio</th><th>Signal</th></tr></thead><tbody>
{''.join(rows)}</tbody></table></div>
<p><a href="solo.nsys-rep">Solo capture</a><a href="{batch_name}.nsys-rep">Batch capture</a><a href="{batch_name}.ncu-rep">Batch Nsight Compute</a><a href="solo-stats_cuda_gpu_kern_sum.csv">Solo CSV</a><a href="{batch_name}-stats_cuda_gpu_kern_sum.csv">Batch CSV</a></p>
</body></html>"""
    path.write_text(document, encoding="utf-8")


def Capture(nsys: Path, executable: Path, output_dir: Path, name: str, batch_size: int, steps: int, logging_interval: int, nlist_interval: int) -> Path:
    prefix = output_dir / name
    Run([
        str(nsys), "profile", "--trace=cuda,nvtx", "--capture-range=cudaProfilerApi",
        "--capture-range-end=stop-shutdown", "--force-overwrite=true", "--output", str(prefix),
        str(executable), "--batch-t4", str(batch_size), str(steps), str(logging_interval), str(nlist_interval),
    ])
    report = prefix.with_suffix(".nsys-rep")
    Run([str(nsys), "stats", "--force-export=true", "--force-overwrite=true", "--report", "cuda_gpu_kern_sum", "--format", "csv", "--output", str(output_dir / f"{name}-stats"), str(report)])
    return output_dir / f"{name}-stats_cuda_gpu_kern_sum.csv"


def CaptureCompute(ncu: Path, executable: Path, output_dir: Path, batch_name: str, batch_size: int, steps: int, logging_interval: int, nlist_interval: int) -> None:
    Run([
        str(ncu), "--set", "full", "--profile-from-start", "off", "--target-processes", "all",
        "--force-overwrite", "--export", str(output_dir / batch_name), str(executable),
        "--batch-t4", str(batch_size), str(steps), str(logging_interval), str(nlist_interval),
    ])


def CaptureBaselineCompute(ncu: Path, executable: Path, output_dir: Path) -> None:
    Run([
        str(ncu), "--set", "full", "--profile-from-start", "on", "--target-processes", "all",
        "--force-overwrite", "--export", str(output_dir / "baseline"), str(executable), "--baseline",
    ])


def CaptureProfileBaselineCompute(ncu: Path, executable: Path, output_dir: Path, logging_interval: int, nlist_interval: int) -> None:
    Run([
        str(ncu), "--set", "full", "--profile-from-start", "off", "--target-processes", "all",
        "--force-overwrite", "--export", str(output_dir / "profile-baseline"), str(executable),
        "--batch-t4", "1", "4", str(logging_interval), str(nlist_interval),
    ])


def Main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--nsys", type=Path, help="Path to nsys.exe")
    parser.add_argument("--ncu", type=Path, help="Path to ncu.exe")
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=1, help="Number of systems in the batched capture")
    parser.add_argument("--capture-baseline", action="store_true", help="Capture only the isolated four-step NCU baseline")
    parser.add_argument("--capture-profile-baseline", action="store_true", help="Capture only the four-step ProfileT4 NCU baseline")
    parser.add_argument("--logging-interval", type=int, default=40)
    parser.add_argument("--nlist-interval", type=int, default=20)
    parser.add_argument("--skip-build", action="store_true", help="Use the existing limaprofile executable")
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be at least 1")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    #if not args.skip_build:
    #    BuildProfileExecutable(args.build_dir)
    build_dirs = (args.build_dir.resolve(), ROOT / "build" / "x64-Release")
    executable = next((build_dir / "code" / "LIMA_TESTS" / "limaprofile.exe"
        for build_dir in build_dirs if (build_dir / "code" / "LIMA_TESTS" / "limaprofile.exe").is_file()),
        build_dirs[0] / "code" / "LIMA_TESTS" / "limaprofile.exe")
    if not executable.is_file():
        raise FileNotFoundError(f"Profiling executable was not built: {executable}")
    if args.capture_baseline:
        CaptureBaselineCompute(FindNcu(args.ncu), executable, output_dir)
        print(f"Wrote {output_dir / 'baseline.ncu-rep'}")
        return
    if args.capture_profile_baseline:
        CaptureProfileBaselineCompute(FindNcu(args.ncu), executable, output_dir, args.logging_interval, args.nlist_interval)
        print(f"Wrote {output_dir / 'profile-baseline.ncu-rep'}")
        return
    nsys = FindNsys(args.nsys)
    ncu = FindNcu(args.ncu)
    solo_csv = Capture(nsys, executable, output_dir, "solo", 1, args.steps, args.logging_interval, args.nlist_interval)
    batch_name = f"batch{args.batch_size}"
    batch_csv = Capture(nsys, executable, output_dir, batch_name, args.batch_size, args.steps, args.logging_interval, args.nlist_interval)
    report = output_dir / f"T4_KERNEL_COMPARISON_{batch_name}.html"
    WriteReport(report, ReadKernelSummary(solo_csv), ReadKernelSummary(batch_csv), args.steps, args.batch_size, batch_name)

    computeSteps = 4 # compute is muuuch slower
    CaptureCompute(ncu, executable, output_dir, batch_name, args.batch_size, computeSteps, args.logging_interval, args.nlist_interval)
    print(f"Wrote {report}")


if __name__ == "__main__":
    Main()
