import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_binary_float_matrices(num_slices, gridpoints_per_dim, center_slice, spacing, filename="C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin"):
    with open(filename, "rb") as file:
        data = np.frombuffer(file.read(), dtype=np.float32)


    data = data.reshape(num_slices, gridpoints_per_dim, gridpoints_per_dim)
    slices = [data[i] for i in range(data.shape[0])]

    indices = [center_slice + i * spacing for i in range(-num_slices // 2, num_slices // 2 + 1)]

    vmin = max(data.min(), 1e-10)  # Avoid log(0) issues
    vmax = data.max()

    fig, axes = plt.subplots(num_slices, 1, figsize=(15, num_slices * 5))
#    axes = axes.flatten()

    for ax, slice_data, slice_idx in zip(axes, slices, indices):
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        #im = ax.imshow(slice_data, cmap='viridis', aspect='equal', norm=norm)
        im = ax.imshow(slice_data, cmap='viridis', aspect='equal', vmin=vmin, vmax=vmax)
        ax.set_title(f"Slice {slice_idx}")
        ax.axis('off')
        fig.colorbar(im, ax=ax, shrink=0.6)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python plot_slices.py <num_slices> <gridpoints_per_dim> <center_slice> <spacing>")
        sys.exit(1)

    num_slices = int(sys.argv[1])
    gridpoints_per_dim = int(sys.argv[2])
    center_slice = int(sys.argv[3])
    spacing = int(sys.argv[4])

    plot_binary_float_matrices(num_slices, gridpoints_per_dim, center_slice, spacing)
