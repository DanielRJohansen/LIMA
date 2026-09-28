import argparse
import csv
from pathlib import Path
import sys

import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def PlotDensityProfile(csv_path: Path, show: bool) -> None:
    profile = np.genfromtxt(csv_path, delimiter=",", names=True)
    figure, axis = plt.subplots(figsize=(8, 4.5))
    axis.plot(profile["z_nm"], profile["water_density"], label="Water", color="#3b82f6", linewidth=2)
    axis.plot(profile["z_nm"], profile["head_density"], label="Lipid heads", color="#f59e0b", linewidth=2)
    axis.plot(profile["z_nm"], profile["tail_density"], label="Lipid tails", color="#8b5cf6", linewidth=2)
    axis.set_xlabel("Position through membrane (nm)")
    axis.set_ylabel("Number density (nm$^{-3}$)")
    axis.set_title(csv_path.parent.name.replace("_", " "))
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(csv_path.with_suffix(".png"), dpi=180)
    if show:
        plt.show()
    plt.close(figure)


def PlotDensityComparison(csv_path: Path, show: bool) -> None:
    groups = {}
    with csv_path.open(newline="") as file:
        for row in csv.DictReader(file):
            key = (row["composition"], float(row["temperature"]))
            groups.setdefault(key, []).append(row)

    compositions = list(dict.fromkeys(composition for composition, _ in groups))
    temperatures = sorted({temperature for _, temperature in groups})
    figure, axes = plt.subplots(len(compositions), len(temperatures), figsize=(6 * len(temperatures), 4 * len(compositions)), sharex=True, sharey=True, squeeze=False)
    for row, composition in enumerate(compositions):
        for column, temperature in enumerate(temperatures):
            axis = axes[row][column]
            profile = groups[(composition, temperature)]
            z = [float(point["z_nm"]) for point in profile]
            axis.plot(z, [float(point["water_density"]) for point in profile], label="Water", color="#3b82f6", linewidth=2)
            axis.plot(z, [float(point["head_density"]) for point in profile], label="Lipid heads", color="#f59e0b", linewidth=2)
            axis.plot(z, [float(point["tail_density"]) for point in profile], label="Lipid tails", color="#8b5cf6", linewidth=2)
            axis.set_title(f"{composition}: {temperature:.0f} K")
            axis.grid(alpha=0.25)
    for axis in axes[-1]:
        axis.set_xlabel("Position through membrane (nm)")
    for axis in axes[:, 0]:
        axis.set_ylabel("Number density (nm$^{-3}$)")
    axes[0][0].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(csv_path.with_suffix(".png"), dpi=180)
    if show:
        plt.show()
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path, nargs="?")
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    if args.comparison:
        PlotDensityComparison(args.comparison, args.show)
    elif args.csv_path:
        PlotDensityProfile(args.csv_path, args.show)
    else:
        parser.error("provide a density profile CSV or --comparison")
