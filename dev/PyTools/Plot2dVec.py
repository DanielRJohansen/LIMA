import numpy as np
import matplotlib.pyplot as plt
import sys


def read_binary_float_matrices(filename, num_slices, gridpoints_per_dim):
    with open(filename, "rb") as file:
        data = np.frombuffer(file.read(), dtype=np.float32)
    return data.reshape(num_slices, gridpoints_per_dim, gridpoints_per_dim)


def plot_colormaps(axes, slices, indices, data_min, data_max):
    for ax, slice_data, slice_idx in zip(axes, slices, indices):
        im = ax.imshow(slice_data, cmap='viridis', aspect='equal', vmin=data_min, vmax=data_max)
        ax.set_title(f"Slice {slice_idx}")
        ax.axis('off')
        plt.colorbar(im, ax=ax, shrink=0.6)


def extract_line_plots(data):
    num_slices = data.shape[0]
    x_lines = []
    y_lines = []

    for slice_idx in range(num_slices):
        slice_data = data[slice_idx]
        x_idx = np.argmax(slice_data.max(axis=1))  # Column index of max
        y_idx = np.argmax(slice_data.max(axis=0))  # Row index of max

        x_line = slice_data[x_idx, :]  # Horizontal line (row at max row index)
        y_line = slice_data[:, y_idx]  # Vertical line (column at max column index)

        x_lines.append(x_line)
        y_lines.append(y_line)

    return x_lines, y_lines


def plot_line_plots(axes, x_lines, y_lines, indices, data_min, data_max):
    for ax, x_line, y_line, slice_idx in zip(axes, x_lines, y_lines, indices):
        ax.plot(x_line, label="X-Dim Line")
        ax.plot(y_line, label="Y-Dim Line")
        ax.set_title(f"Line Plot for Slice {slice_idx}")
        ax.legend()
        ax.set_ylim(data_min, data_max)  # Set the same y-limits for all plots



def plot_binary_float_matrices(num_slices, gridpoints_per_dim, center_slice, spacing,
                               filename="C:/Users/Daniel/git_repo/LIMA_data/Pool/PmePot_AllSlices.bin"):
    data = read_binary_float_matrices(filename, num_slices, gridpoints_per_dim)
    slices = [data[i] for i in range(data.shape[0])]

    indices = [center_slice + (i - num_slices // 2) * spacing for i in range(num_slices)]
    print(indices)
    x_lines, y_lines = extract_line_plots(data)

    # Create subplots with an extra column for line plots
    fig, axes = plt.subplots(num_slices, 2, figsize=(15, num_slices * 5))
    data_min, data_max = data.min(), data.max()

    # Ensure axes is always iterable
    axes = axes if num_slices > 1 else [axes]

    # Plot colormaps in the first column
    plot_colormaps([ax[0] for ax in axes], slices, indices, data_min, data_max)

    # Plot line plots in the second column
    plot_line_plots([ax[1] for ax in axes], x_lines, y_lines, indices, data_min, data_max)

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

    print("Num slices", num_slices, " Center slice ", center_slice)

    plot_binary_float_matrices(num_slices, gridpoints_per_dim, center_slice, spacing)
