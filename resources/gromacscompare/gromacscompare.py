"""lima gromacscompare: LIMA vs GROMACS on the GROMACS heterogeneous-parallelization benchmark systems (Pall et al.,
J. Chem. Phys. 153, 134110 (2020)), for force parity, structure and performance. Writes report.html (report_page.py, the
same page as LIMA's published benchmark report) and results.json; nothing is pass/fail.

Run as `lima gromacscompare [--systems SET] [--structure] [-o DIR]`, which runs this script. Needs Python 3.8+
(standard library only) and GROMACS with CUDA: natively, or in WSL on Windows.

Per system, in <output-dir>/<system>/:
- forces/   step-0 forces of both engines from the same conf.gro/topol.top (single point, no constraints), compared as
            relative RMS error sqrt(sum |F_lima - F_gmx|^2 / sum |F_gmx|^2), overall and for the worst residue name
- perf/     ns/day of both engines, alternating, on this machine
- md/       only with --structure: a 0.2 ns LIMA run and two 0.2 ns GROMACS runs differing only in the thermostat seed.
            For each trajectory gmx makes the protein whole and unwrapped, then the C-alpha RMSD to the start and the
            per-atom C-alpha RMSF over the second half. LIMA's RMSF is compared to GROMACS' by Pearson correlation, next to
            the correlation between the two GROMACS runs, ie. how much runs differ by chance. Without --structure the report
            shows these numbers from the 1 ns reference runs (reference_parity.json), marked as such

The inputs come from https://lima-dynamics.com/benchdata/v1, made by
LIMA_data/salad_gromacs_bench/bundle/make_bundle.py, whose README describes it."""
import argparse
import hashlib
import html
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path, PureWindowsPath

import report_page

HERE = Path(__file__).parent
BASE_URL = os.environ.get("LIMA_BENCHDATA_URL", "https://lima-dynamics.com/benchdata/v1")
DATASET = "v1"
SYSTEMS = {  # name: (tier, atoms, bytes, sha256), from make_bundle.py
    "villin": ("small", 5006, 212764, "585214930d7b408316c42cf4b4b6c43f8fab182a75188752f0c2ba2bf27d5b64"),
    "rnase_cubic": ("small", 24040, 918469, "911047bc119d34d391e5b4c367c665a216ddf916d9d975d816329bafd49d8dbd"),
    "adh_dodec": ("all", 191122, 7433786, "875492bb2a1186190dfd30c680c250fc0912d77b98c18c049229030ba1baaf84"),
    "aqp_ensemble": ("all", 219984, 9744534, "894c6e817ab92d72d283dd1cce7d468ca7766b5e2a1ac939821c827faa19a9e8"),
    "ion_channel": ("all", 283354, 10049822, "63bf9a45f3880072147b37c71b250c313a006d3d9adaab7c56fcdc0c70588dad"),
    "stmv": ("all", 1066628, 54429044, "e5f51dd3da6758931f905ab4f51c72c7545d7d6ecc29be76f69d368e275cbfb8"),
}
# Steps per performance run, ~30-60 s each on an RTX 4090 (as LIMA_data/salad_gromacs_bench/bench.sh)
PERF_STEPS = {"villin": 150000, "rnase_cubic": 120000, "adh_dodec": 30000, "aqp_ensemble": 25000, "ion_channel": 20000, "stmv": 6000}
MD_STEPS = 100000           # 0.2 ns at 2 fs
MD_FRAME_INTERVAL = 1000    # 2 ps
PERF_ROUNDS = 2
DT_PS = 0.002


class Failure(Exception):
    pass


def log(message=""):
    print(message, flush=True)


def lima_mdp_with(mdp, **values):
    """An .mdp's text with keys replaced or added, eg. nsteps."""
    lines, seen = [], set()
    for line in Path(mdp).read_text().splitlines():
        key = line.split("=")[0].strip().lower().replace("_", "-")
        if key in values:
            seen.add(key)
            line = f"{key} = {values[key]}"
        lines.append(line)
    lines += [f"{key} = {value}" for key, value in values.items() if key not in seen]
    return "\n".join(lines) + "\n"


def write_mdp(path, source, **values):
    path.write_text(lima_mdp_with(source, **{k.replace("_", "-"): v for k, v in values.items()}))
    return path


# ------------------------------------------------ Data ------------------------------------------------ #

def default_data_dir():
    if os.environ.get("LIMA_BENCHDATA_DIR"):
        return Path(os.environ["LIMA_BENCHDATA_DIR"])
    if os.name == "nt" and os.environ.get("LOCALAPPDATA"):
        return Path(os.environ["LOCALAPPDATA"]) / "LIMA" / "benchdata"
    cache = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(cache) / "lima" / "benchdata"


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(system, data_dir):
    """The system's directory, downloaded, verified and extracted on first use."""
    _, _, size, sha256 = SYSTEMS[system]
    root = data_dir / DATASET
    target = root / system
    marker = target / ".sha256"
    if marker.exists() and marker.read_text().strip() == sha256:
        return target

    root.mkdir(parents=True, exist_ok=True)
    archive = root / f"{system}.tar.gz.part"
    url = f"{BASE_URL.rstrip('/')}/{system}.tar.gz"
    log(f"Downloading {url} ({size / 1e6:.1f} MB)")
    try:
        with urllib.request.urlopen(url, timeout=60) as response, open(archive, "wb") as out:
            shutil.copyfileobj(response, out, 1 << 20)
    except OSError as error:
        archive.unlink(missing_ok=True)
        raise Failure(f"Downloading {url} failed: {error}")
    if archive.stat().st_size != size or sha256_of(archive) != sha256:
        archive.unlink()
        raise Failure(f"{url} does not match the expected size and SHA-256: the download is damaged, or the dataset "
                      f"changed and this LIMA is outdated")

    staging = root / f".extract-{system}"
    shutil.rmtree(staging, ignore_errors=True)
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            name = PureWindowsPath(member.name)
            if not (member.isfile() or member.isdir()) or name.is_absolute() or ".." in name.parts or name.parts[0] != system:
                raise Failure(f"{archive.name} contains an unexpected entry {member.name}")
        tar.extractall(staging, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
    shutil.rmtree(target, ignore_errors=True)
    (staging / system).rename(target)
    staging.rmdir()
    archive.unlink()
    marker.write_text(sha256 + "\n")
    return target


# ------------------------------------------------ Engines ------------------------------------------------ #

def run(command, cwd, logfile, stdin=None):
    """Runs a command with its output in logfile; raises Failure with the log's tail if it fails."""
    with open(logfile, "w", errors="replace") as out:
        result = subprocess.run(command, cwd=cwd, stdout=out, stderr=subprocess.STDOUT, input=stdin,
                                text=True, errors="replace")
    if result.returncode != 0:
        tail = "\n".join(Path(logfile).read_text(errors="replace").replace("\0", "").splitlines()[-15:])
        raise Failure(f"'{' '.join(map(str, command))}' failed with exit code {result.returncode}; see {logfile}:\n{tail}")
    return Path(logfile).read_text(errors="replace").replace("\0", "")


class Gromacs:
    """gmx, natively or in WSL. Every file argument is relative to cwd, so the paths work in both."""

    def __init__(self, path, wsl_distro):
        self.path = path
        self.wsl = wsl_distro

    def command(self, args, cwd):
        if self.wsl is None:
            return [self.path] + args
        distro = ["-d", self.wsl] if self.wsl else []
        return ["wsl.exe"] + distro + ["--cd", str(cwd), "-e", self.path] + args

    def run(self, args, cwd, logfile, stdin=None):
        return run(self.command(args, cwd), cwd, logfile, stdin)

    def describe(self):
        return f"{self.path} (in WSL{' ' + self.wsl if self.wsl else ''})" if self.wsl is not None else self.path


def find_gromacs(option, scratch):
    """The GROMACS to use and its version line. Fails unless it is a CUDA build."""
    candidates = []
    if option:
        if option.startswith("wsl:"):
            candidates.append(Gromacs(option[4:], ""))
        else:
            candidates.append(Gromacs(option, None))
    else:
        for name in ("gmx", "gmx_mpi"):
            if shutil.which(name):
                candidates.append(Gromacs(shutil.which(name), None))
        if os.name == "nt" and shutil.which("wsl.exe"):
            for name in ("gmx", "gmx_mpi"):
                found = subprocess.run(["wsl.exe", "-e", "bash", "-lc", f"command -v {name}"],
                                       capture_output=True, text=True, errors="replace")
                if found.returncode == 0 and found.stdout.strip():
                    candidates.append(Gromacs(found.stdout.strip().splitlines()[-1], ""))
    if not candidates:
        raise Failure("GROMACS was not found. Install a GROMACS built with CUDA, or point to it with --gmx PATH "
                      "(--gmx wsl:/path/to/gmx for one in WSL, whose PATH is often not set up for non-interactive shells)")

    reasons = []
    for gmx in candidates:
        try:
            text = gmx.run(["--version"], scratch, scratch / "gmx_version.log")
        except (Failure, OSError) as error:
            reasons.append(f"{gmx.describe()}: {error}")
            continue
        version = re.search(r"GROMACS version:\s*(\S+)", text)
        gpu = re.search(r"GPU support:\s*(\S+)", text)
        precision = re.search(r"Precision:\s*(\S+)", text)
        if gpu and gpu.group(1).upper() == "CUDA":
            return gmx, {"version": version.group(1) if version else "?", "precision": precision.group(1) if precision else "?"}
        reasons.append(f"{gmx.describe()}: GPU support is {gpu.group(1) if gpu else 'unknown'}, not CUDA")
    raise Failure("No GROMACS with CUDA was found. The comparison runs GROMACS on the GPU, as it is normally run.\n  "
                  + "\n  ".join(reasons))


def find_lima(option):
    if option:
        return option
    found = shutil.which("lima")
    if not found:
        raise Failure("lima was not found on the PATH; give it with --lima PATH")
    return found


def lima_mdrun(lima, cwd, inputs, mdp, logfile, *extra):
    text = run([lima, "--no-daemon", "mdrun", "-c", str(inputs / "conf.gro"), "-t", str(inputs / "topol.top"), "-s", str(mdp),
                "--conf-out", "out.gro"] + list(extra), cwd, logfile)
    seconds = re.search(r"^\s*Time:\s*([\d.]+)", text, re.MULTILINE)
    if not seconds:
        raise Failure(f"No engine time in {logfile}")
    return float(seconds.group(1))


def gmx_grompp(gmx, cwd, mdp, tpr, logfile):
    gmx.run(["grompp", "-f", mdp, "-c", "../inputs/conf.gro", "-p", "../inputs/topol.top", "-o", tpr, "-maxwarn", "5",
             "-quiet"], cwd, logfile)


GMX_MDRUN_GPU = ["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu", "-ntmpi", "1", "-notunepme"]


def gmx_mdrun_flags(gmx):
    # Thread pinning is slower under WSL
    return GMX_MDRUN_GPU + (["-pin", "off"] if gmx.wsl is not None else [])


# ------------------------------------------------ Measurements ------------------------------------------------ #

def read_forces_csv(path):
    with open(path) as f:
        return [tuple(map(float, line.split(","))) for line in f if line.strip()]


def read_xvg(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line and line[0] not in "#@":
            rows.append([float(v) for v in line.split()])
    return rows


def residue_names(gro):
    lines = Path(gro).read_text().splitlines()
    return [line[5:10].strip() for line in lines[2:2 + int(lines[1])]]


def force_parity(system, data, work, lima, gmx):
    work.mkdir(parents=True, exist_ok=True)
    # The reference single point: no constraints, the system's cutoff and PME spacing, nstfout 1
    shutil.copy(data / "reference" / "forces" / "forces.mdp", work / "gmx_forces.mdp")
    gmx_grompp(gmx, work, "gmx_forces.mdp", "gmx_forces.tpr", work / "gmx_grompp.log")
    gmx.run(["mdrun", "-s", "gmx_forces.tpr", "-deffnm", "gmx_forces", "-quiet", "-nb", "gpu", "-pme", "gpu", "-ntmpi", "1"],
            work, work / "gmx_mdrun.log")
    gmx.run(["traj", "-s", "gmx_forces.tpr", "-f", "gmx_forces.trr", "-of", "gmx_forces.xvg", "-quiet"], work,
            work / "gmx_traj.log", stdin="System\n")
    values = read_xvg(work / "gmx_forces.xvg")[0][1:]
    reference = [tuple(values[i:i + 3]) for i in range(0, len(values), 3)]

    mdp = write_mdp(work / "lima_forces.mdp", data / "inputs" / "lima.mdp", nsteps=1, constraints="none")
    lima_mdrun(lima, work, work.parent / "inputs", mdp, work / "lima_mdrun.log", "--forces-out", "lima_forces.csv")
    forces = read_forces_csv(work / "lima_forces.csv")

    residues = residue_names(data / "inputs" / "conf.gro")
    if not (len(forces) == len(reference) == len(residues)):
        raise Failure(f"{system}: {len(forces)} LIMA forces, {len(reference)} GROMACS forces, {len(residues)} atoms")
    error = norm = 0.0
    per_residue = {}
    for f, g, residue in zip(forces, reference, residues):
        e = sum((a - b) ** 2 for a, b in zip(f, g))
        r = sum(b * b for b in g)
        error += e
        norm += r
        acc = per_residue.setdefault(residue, [0.0, 0.0])
        acc[0] += e
        acc[1] += r
    worst = max(per_residue.items(), key=lambda kv: kv[1][0] / kv[1][1] if kv[1][1] > 0 else 0)
    return {"relativeRmsError": math.sqrt(error / norm), "worstResidue": worst[0],
            "worstResidueRelativeRmsError": math.sqrt(worst[1][0] / worst[1][1])}


def performance(system, data, work, lima, gmx):
    """ns/day of both engines, alternating per round so clock drift hits both alike. LIMA from the engine-time difference
    of a 1/5 and a full run, so setup and the first steps are excluded; GROMACS with -resethway."""
    work.mkdir(parents=True, exist_ok=True)
    steps = PERF_STEPS[system]
    short_mdp = write_mdp(work / "lima_short.mdp", data / "inputs" / "lima.mdp", nsteps=steps // 5)
    full_mdp = write_mdp(work / "lima_full.mdp", data / "inputs" / "lima.mdp", nsteps=steps)
    write_mdp(work / "gmx.mdp", data / "inputs" / "gmx.mdp", nsteps=steps)
    gmx_grompp(gmx, work, "gmx.mdp", "gmx.tpr", work / "gmx_grompp.log")
    inputs = work.parent / "inputs"
    results = {"steps": steps, "lima": [], "gromacs": []}
    for round_ in range(1, PERF_ROUNDS + 1):
        short = lima_mdrun(lima, work, inputs, short_mdp, work / f"lima_short_{round_}.log")
        full = lima_mdrun(lima, work, inputs, full_mdp, work / f"lima_full_{round_}.log")
        results["lima"].append((steps - steps // 5) * DT_PS * 1e-3 / ((full - short) / 86400))
        gmx.run(["mdrun", "-s", "gmx.tpr", "-deffnm", f"gmx_{round_}", "-quiet", "-resethway", "-noconfout"] + gmx_mdrun_flags(gmx),
                work, work / f"gmx_mdrun_{round_}.log")
        text = (work / f"gmx_{round_}.log").read_text(errors="replace")
        match = re.search(r"^Performance:\s*([\d.]+)", text, re.MULTILINE)
        if not match:
            raise Failure(f"No performance line in {work / f'gmx_{round_}.log'}")
        results["gromacs"].append(float(match.group(1)))
        threads = re.search(r"^Using (\d+) OpenMP threads", text, re.MULTILINE)
        if threads:
            results["gromacsOpenmpThreads"] = int(threads.group(1))
        log(f"  round {round_}: LIMA {results['lima'][-1]:.1f}, GROMACS {results['gromacs'][-1]:.1f} ns/day")
    for engine in ("lima", "gromacs"):
        values = sorted(results[engine])
        results[engine + "Median"] = values[len(values) // 2] if len(values) % 2 else (values[len(values) // 2 - 1] + values[len(values) // 2]) / 2
    return results


def analyze_trajectory(gmx, work, tpr, trajectory, name):
    """C-alpha RMSD to the start (fitted) and per-atom RMSF over the second half, as LIMA_data's parity/analyze.sh."""
    gmx.run(["trjconv", "-s", tpr, "-f", trajectory, "-pbc", "whole", "-o", f"{name}_whole.xtc", "-quiet"], work,
            work / f"{name}_trjconv_whole.log", stdin="System\n")
    gmx.run(["trjconv", "-s", tpr, "-f", f"{name}_whole.xtc", "-pbc", "nojump", "-o", f"{name}_nojump.xtc", "-quiet"], work,
            work / f"{name}_trjconv_nojump.log", stdin="System\n")
    gmx.run(["rms", "-s", tpr, "-f", f"{name}_nojump.xtc", "-o", f"{name}_rmsd.xvg", "-quiet"], work,
            work / f"{name}_rms.log", stdin="C-alpha\nC-alpha\n")
    end = read_xvg(work / f"{name}_rmsd.xvg")[-1][0]
    gmx.run(["rmsf", "-s", tpr, "-f", f"{name}_nojump.xtc", "-o", f"{name}_rmsf.xvg", "-b", str(end / 2), "-quiet"], work,
            work / f"{name}_rmsf.log", stdin="C-alpha\n")
    for suffix in ("whole", "nojump"):
        (work / f"{name}_{suffix}.xtc").unlink(missing_ok=True)
    return work / f"{name}_rmsd.xvg", work / f"{name}_rmsf.xvg"


def pearson(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))


def structure(system, data, work, lima, gmx):
    work.mkdir(parents=True, exist_ok=True)
    inputs = work.parent / "inputs"
    # The .mdp of the dataset's two GROMACS runs (they differ only in ld-seed), shortened. The tpr of the first also gives
    # the topology for analyzing LIMA's trajectory
    xvgs = {}
    for seed in (1, 2):
        write_mdp(work / f"gmx_seed{seed}.mdp", data / "reference" / f"md_seed{seed}" / "md.mdp", nsteps=MD_STEPS,
                  nstxout_compressed=MD_FRAME_INTERVAL)
        gmx_grompp(gmx, work, f"gmx_seed{seed}.mdp", f"gmx_seed{seed}.tpr", work / f"gmx_seed{seed}_grompp.log")
        log(f"  GROMACS {MD_STEPS * DT_PS * 1e-3:g} ns, seed {seed}")
        gmx.run(["mdrun", "-s", f"gmx_seed{seed}.tpr", "-deffnm", f"gmx_seed{seed}", "-quiet"] + gmx_mdrun_flags(gmx),
                work, work / f"gmx_seed{seed}_mdrun.log")
        xvgs[f"gmx{seed}"] = analyze_trajectory(gmx, work, f"gmx_seed{seed}.tpr", f"gmx_seed{seed}.xtc", f"gmx_seed{seed}")

    log(f"  LIMA {MD_STEPS * DT_PS * 1e-3:g} ns")
    mdp = write_mdp(work / "lima_md.mdp", data / "inputs" / "lima.mdp", nsteps=MD_STEPS)
    seconds = lima_mdrun(lima, work, inputs, mdp, work / "lima_md_mdrun.log", "--trajectory", "lima_md.xtc",
                         "--set", f"nstxout_compressed={MD_FRAME_INTERVAL}")
    xvgs["lima"] = analyze_trajectory(gmx, work, "gmx_seed1.tpr", "lima_md.xtc", "lima_md")

    result = {"simulatedNs": MD_STEPS * DT_PS * 1e-3, "limaNsPerDay": MD_STEPS * DT_PS * 1e-3 / (seconds / 86400)}
    rmsf = {}
    for engine, (rmsd_path, rmsf_path) in xvgs.items():
        rmsd = read_xvg(rmsd_path)
        end = rmsd[-1][0]
        late = [r for t, r in rmsd if t >= end / 2]
        result[f"rmsdNm_{engine}"] = sum(late) / len(late)
        rmsf[engine] = [r for _, r in read_xvg(rmsf_path)]
    if not (len(rmsf["lima"]) == len(rmsf["gmx1"]) == len(rmsf["gmx2"])):
        raise Failure(f"{system}: RMSF over {len(rmsf['lima'])}, {len(rmsf['gmx1'])} and {len(rmsf['gmx2'])} C-alpha atoms")
    result["calphaAtoms"] = len(rmsf["lima"])
    result["rmsfCorrelationLimaGromacs"] = pearson(rmsf["gmx1"], rmsf["lima"])
    result["rmsfCorrelationGromacsGromacs"] = pearson(rmsf["gmx1"], rmsf["gmx2"])
    return result


def reference_structure(system):
    """structure()'s numbers from the reference runs, for a report without --structure"""
    reference = json.loads((HERE / "reference_parity.json").read_text())
    r = reference["systems"][system]
    return {"reference": f"LIMA {reference['lima']}, GROMACS {reference['gromacs']}, {reference['gpu']}",
            "simulatedNs": reference["simulatedNs"], "calphaAtoms": r["calphaAtoms"], "rmsdNm_lima": r["rmsdNm_lima"],
            "rmsdNm_gmx1": r["rmsdNm_gmx"], "rmsdNm_gmx2": r["rmsdNm_gmx2"], "rmsfCorrelationLimaGromacs": r["rmsfCorrelation"],
            "rmsfCorrelationGromacsGromacs": r["rmsfCorrelationGmxReference"]}


# ------------------------------------------------ Report ------------------------------------------------ #

REPORT_CSS = ".bars .error { grid-column: 1 / -1; margin: 0; color: var(--soft); font-size: 12px; }"
STRUCTURE_MARK = "†"


def report_rows(results):
    rows = []
    for system, r in results.items():
        row = {"system": system, "atoms": SYSTEMS[system][1]}
        if "performance" in r:
            row.update(lima=r["performance"]["limaMedian"], gromacs=r["performance"]["gromacsMedian"])
        if "forces" in r:
            row["forces"] = r["forces"]["relativeRmsError"]
        if "structure" in r:
            s = r["structure"]
            row.update(rmsdLima=s["rmsdNm_lima"], rmsdGmx1=s["rmsdNm_gmx1"], rmsdGmx2=s["rmsdNm_gmx2"],
                       rmsfR=s["rmsfCorrelationLimaGromacs"], rmsfRGmx=s["rmsfCorrelationGromacsGromacs"])
            if "reference" in s:
                row["structureMark"] = STRUCTURE_MARK
        if "error" in r:
            row["error"] = r["error"]
        rows.append(row)
    return rows


def method_html(results, setup):
    """The footer: this machine's setup, the physics, and what the parity numbers mean"""
    threads = next((r["performance"]["gromacsOpenmpThreads"] for r in results.values()
                    if "gromacsOpenmpThreads" in r.get("performance", {})), None)
    gmx = setup["gromacs"]
    structures = [r["structure"] for r in results.values() if "structure" in r]
    measured = [s for s in structures if "reference" not in s]
    md_ns = (measured or structures or [{"simulatedNs": MD_STEPS * DT_PS * 1e-3}])[0]["simulatedNs"]
    reference = next((s["reference"] for s in structures if "reference" in s), None)
    esc = html.escape
    parts = [f"""
    <div>
      <h3>Setup</h3>
      <p>{esc(setup['gpu'])}, {esc(setup['cpu'])}, {esc(setup['os'])}. {esc(setup['limaVersion'])}. GROMACS {esc(gmx['version'])} with CUDA{', under WSL' if setup['gromacsWsl'] else ''},
        everything on the GPU (-nb -pme -bonded -update gpu){f', {threads} OpenMP threads' if threads else ''}.</p>
      <p>Median of {PERF_ROUNDS} runs per engine, run alternately. Throughput from the steady state: GROMACS -resethway, LIMA a long minus a short run.</p>
    </div>
    <div>
      <h3>Same physics</h3>
      <p>2 fs steps, h-bond constraints, PME with each system's original cutoff (0.9-1.2 nm) and grid spacing, v-rescale thermostat.
        LIMA integrates with velocity Verlet, GROMACS with leap-frog, its fastest path. GROMACS keeps the cutoff fixed (-notunepme).</p>
      <p>Systems: Páll et al., J. Chem. Phys. 153, 134110 (2020), force fields as published. ADH, aquaporin and GluCl run as
        2 copies, the exact rectangular cell of their triclinic box.</p>
    </div>
    <div>
      <h3>Parity</h3>
      <p><b>Forces</b>: the relative RMS difference of every atom's force, LIMA vs GROMACS, for the same structure and settings.
        1×10<sup>−4</sup> means a typical force differs by 0.01%, ten times less than PME's own approximation of the exact Ewald sum (~10<sup>−3</sup>).</p>
      <p><b>RMSD</b>: how far the protein's Cα atoms drift from the start, the mean over the second half of a {md_ns:g} ns run, next to two GROMACS
        runs. <b>RMSF r</b>: how well LIMA and GROMACS agree on which residues move most (Pearson correlation of per-Cα fluctuations,
        1 is identical), next to the agreement between the two GROMACS runs, which differ only in their thermostat seed.</p>"""]
    if reference:
        parts.append(f"""
      <p><b>{STRUCTURE_MARK}</b> RMSD and RMSF r from the reference runs ({esc(reference)}), not measured on this machine.
        <b>lima gromacscompare --structure</b> measures them.</p>""")
    parts.append("""
    </div>
  """)
    return "".join(parts)


def write_report(out_dir, results, setup):
    (out_dir / "results.json").write_text(json.dumps({"setup": setup, "systems": results}, indent=2) + "\n")
    rows = report_rows(results)
    rows_html = [report_page.row(r, report_page.data_uri(HERE / "assets" / f"{r['system']}.webp", "image/webp")) for r in rows]
    page = report_page.page("dark", "page", REPORT_CSS, (HERE / "assets" / "logo.svg").read_text(encoding="utf-8"),
                            f"Benchmark report · {setup['date']}", report_page.headline(rows),
                            "A one-to-one comparison on GROMACS' own benchmark systems, run on this machine.",
                            report_page.stats(rows), report_page.table(rows_html), "", method_html(results, setup))
    (out_dir / "report.html").write_text(page, encoding="utf-8")


def summary(results):
    """The report's numbers as text, for the console (ASCII: a Windows console may not have the HTML's mark)"""
    lines = [f"{'System':13s} {'LIMA':>9s} {'GROMACS':>9s} {'Speedup':>8s} {'Forces':>9s} {'RMSF r':>14s}",
             f"{'':13s} {'[ns/day]':>9s} {'[ns/day]':>9s} {'':8s} {'':9s} {'LIMA (GMX)':>14s}"]
    for row in report_rows(results):
        if "error" in row:
            lines.append(f"{row['system']:13s} not completed: {row['error'].splitlines()[0]}")
            continue
        lines.append(f"{row['system']:13s} {row['lima']:9.1f} {row['gromacs']:9.1f} {row['lima'] / row['gromacs']:7.2f}x "
                     f"{row['forces']:9.2e} {row['rmsfR']:6.2f} ({row['rmsfRGmx']:.2f}){'*' if 'structureMark' in row else ''}")
    if any("structureMark" in row for row in report_rows(results)):
        lines.append(f"* from the 1 ns reference runs, not this machine; --structure measures it")
    return "\n".join(lines)


def machine_info():
    gpu = "unknown"
    try:
        smi = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=30)
        if smi.returncode == 0 and smi.stdout.strip():
            name, driver = (smi.stdout.strip().splitlines()[0].split(",") + ["?"])[:2]
            gpu = f"{name.strip()} (driver {driver.strip()})"
    except (OSError, subprocess.SubprocessError):
        pass
    cpu = platform.processor() or platform.machine()
    if os.name == "nt":
        try:
            ps = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                                capture_output=True, text=True, timeout=30)
            cpu = ps.stdout.strip() or cpu
        except (OSError, subprocess.SubprocessError):
            pass
    elif Path("/proc/cpuinfo").exists():
        match = re.search(r"^model name\s*:\s*(.+)$", Path("/proc/cpuinfo").read_text(), re.MULTILINE)
        cpu = match.group(1) if match else cpu
    return gpu, cpu


# ------------------------------------------------ Main ------------------------------------------------ #

def estimate_minutes(systems, with_structure):
    """Rough wall time on an RTX 4090, from the measured ns/day of each system."""
    lima_ns_per_day = {"villin": 2400, "rnase_cubic": 1600, "adh_dodec": 360, "aqp_ensemble": 245, "ion_channel": 230, "stmv": 51}
    gmx_ns_per_day = {"villin": 1100, "rnase_cubic": 920, "adh_dodec": 180, "aqp_ensemble": 96, "ion_channel": 130, "stmv": 33}
    days = lambda ns, nsday: ns / nsday * 24 * 60
    minutes = 0.0
    for s in systems:
        perf_ns = PERF_STEPS[s] * DT_PS * 1e-3
        minutes += PERF_ROUNDS * (days(1.2 * perf_ns, lima_ns_per_day[s]) + days(perf_ns, gmx_ns_per_day[s])) + 0.5  # + forces, setup
        if with_structure:
            md_ns = MD_STEPS * DT_PS * 1e-3
            minutes += days(md_ns, lima_ns_per_day[s]) + 2 * days(md_ns, gmx_ns_per_day[s]) + 1  # + analysis
    return minutes


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="lima gromacscompare", description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--systems", default="small",
                        help="small (villin, rnase_cubic), all (also adh_dodec, aqp_ensemble, ion_channel, stmv), or a "
                             "comma-separated list (default: small)")
    parser.add_argument("--structure", action="store_true",
                        help="Also run both engines for 0.2 ns and compare RMSD and RMSF (default: show those from the "
                             "1 ns reference runs)")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("gromacscompare"),
                        help="Where the runs and report.html go (default: ./gromacscompare)")
    parser.add_argument("--gmx", help=argparse.SUPPRESS)  # Only if gmx is not found on the PATH or in WSL
    parser.add_argument("--lima", help=argparse.SUPPRESS)  # Passed by lima
    return parser.parse_args(argv)


def select_systems(selection):
    if selection in ("small", "all"):
        return [s for s, (tier, *_) in SYSTEMS.items() if selection == "all" or tier == "small"]
    names = [s for s in selection.split(",") if s]
    unknown = [s for s in names if s not in SYSTEMS]
    if unknown or not names:
        raise Failure(f"Unknown system(s) {', '.join(unknown) or selection!r}; choose small, all, or from {', '.join(SYSTEMS)}")
    return [s for s in SYSTEMS if s in names]


def main(argv):
    args = parse_args(argv)
    systems = select_systems(args.systems)
    out_dir = Path(os.path.abspath(args.output_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = default_data_dir()

    gmx, gmx_version = find_gromacs(args.gmx, out_dir)
    lima = find_lima(args.lima)
    log(f"GROMACS {gmx_version['version']} with CUDA: {gmx.describe()}")
    log(f"LIMA: {lima}")
    minutes = estimate_minutes(systems, args.structure)
    log(f"Systems: {', '.join(systems)}{', with structure runs' if args.structure else ''}. Roughly "
        f"{f'{minutes / 60:.1f} h' if minutes > 90 else f'{minutes:.0f} min'} on an RTX 4090; slower GPUs take longer")

    lima_version = subprocess.run([lima, "--version"], capture_output=True, text=True).stdout.strip()
    gpu, cpu = machine_info()
    setup = {"date": time.strftime("%Y-%m-%d %H:%M"), "lima": f"{lima_version} ({lima})", "limaVersion": lima_version or "LIMA",
             "gromacs": gmx_version, "gromacsPath": gmx.describe(), "gromacsWsl": gmx.wsl is not None, "gpu": gpu, "cpu": cpu,
             "os": f"{platform.system()} {platform.release()}"}
    results = {}
    for system in systems:
        log(f"\n== {system} ({SYSTEMS[system][1]:,} atoms)")
        results[system] = {}
        try:
            data = fetch(system, data_dir)
            work = out_dir / system
            shutil.rmtree(work / "inputs", ignore_errors=True)
            # A copy, so both engines read the inputs by paths relative to their run directory, also from WSL
            shutil.copytree(data / "inputs", work / "inputs", ignore=shutil.ignore_patterns("*.bin"))
            log("  forces")
            r = results[system]["forces"] = force_parity(system, data, work / "forces", lima, gmx)
            log(f"  force relative RMS error {r['relativeRmsError']:.2e}")
            log("  performance")
            results[system]["performance"] = performance(system, data, work / "perf", lima, gmx)
            if args.structure:
                r = results[system]["structure"] = structure(system, data, work / "md", lima, gmx)
                log(f"  RMSF correlation LIMA-GROMACS {r['rmsfCorrelationLimaGromacs']:.3f}, "
                    f"GROMACS-GROMACS {r['rmsfCorrelationGromacsGromacs']:.3f}")
            else:
                results[system]["structure"] = reference_structure(system)
        except Failure as error:
            results[system]["error"] = str(error)
            log(f"  {system} not completed: {error}")
        write_report(out_dir, results, setup)  # After every system, so an interrupted run keeps what it has

    log("\n" + summary(results))
    log(f"\nReport: {out_dir / 'report.html'}")
    return 1 if any("error" in r for r in results.values()) else 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Failure as error:
        print(f"lima gromacscompare: {error}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
