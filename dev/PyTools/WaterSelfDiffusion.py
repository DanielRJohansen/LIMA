import argparse
import csv
from pathlib import Path
import sys

import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def ExperimentalDiffusion(temperature):
    # Holz, Heil & Sacco, PCCP 2 (2000) 4740: D = D0 (T/Ts - 1)^gamma, fitted to measurements at 0-100 C [m^2/s]
    return 1.635e-8 * (temperature / 215.05 - 1.0) ** 2.063


def PlotWaterSelfDiffusion(csv_path: Path, show: bool) -> None:
    with csv_path.open(newline="") as file:
        rows = list(csv.DictReader(file))
    temperatures = np.array([float(row["temperature_K"]) for row in rows])
    diffusion = np.array([float(row["D_m2s"]) for row in rows])

    figure, axis = plt.subplots(figsize=(8, 4.5))
    curve = np.linspace(273.15, 373.15, 200)
    axis.plot(curve, ExperimentalDiffusion(curve) * 1e9, label="Experiment (Holz et al. 2000)", color="#64748b", linewidth=2)
    axis.scatter(temperatures, diffusion * 1e9, label=f"LIMA, SPC/E: {len(rows)} simulations", color="#3b82f6", zorder=3)
    axis.set_xlabel("Temperature (K)")
    axis.set_ylabel("Self-diffusion coefficient (10$^{-9}$ m$^2$/s)")
    axis.set_title("Water self-diffusion")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(csv_path.with_suffix(".png"), dpi=180)
    if show:
        plt.show()
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    PlotWaterSelfDiffusion(args.csv_path, args.show)
