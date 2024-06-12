from MakeLipids import MakeLipids
import struct
import matplotlib.pyplot as plt
import numpy as np

def read_histogram_data(filename):
    with open(filename, 'rb') as f:
        # Read the number of steps
        numSteps = struct.unpack('i', f.read(4))[0]

        all_bins = []
        all_counts = []

        # Read bins size
        f.seek(0, 2)  # Move to the end of the file to get file size
        file_size = f.tell()
        f.seek(4)  # Move back to the first histogram data after the steps size
        remaining_size = file_size - 4  # Subtract the size of the steps_size integer

        for _ in range(numSteps):
            # Calculate number of bins based on remaining size
            num_bins = remaining_size // (8 + 4) // numSteps

            # Read bins
            bins = struct.unpack('q' * num_bins, f.read(8 * num_bins))
            all_bins.append(bins)

            # Read counts
            counts = struct.unpack('i' * num_bins, f.read(4 * num_bins))
            all_counts.append(counts)

    return all_bins, all_counts


def plot_histogram(all_bins, all_counts):
    nPlots = len(all_counts)
    fig, axes = plt.subplots(nPlots, 1)

    if nPlots == 1:
        axes = [axes]

    for ax, bins, counts in zip(axes, all_bins, all_counts):
        print(counts)
        x = np.arange(len(bins))
        ax.bar(x, counts, width=0.8, align='center', alpha=0.7)
        ax.set_yscale('log')  # Set the y-axis to a logarithmic scale
        ax.set_xticks(x)
        ax.set_xticklabels(bins, rotation=45)
        ax.set_xlabel('Bins')
        ax.set_ylabel('Frequency')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    #MakeLipids

    bins, counts = read_histogram_data('C:/Users/Daniel/git_repo/LIMA_data/psome2/histogram_data.bin')
    plot_histogram(bins, counts)