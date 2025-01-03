import matplotlib.pyplot as plt
import numpy as np
import struct
import os
from pathlib import Path

def read_binary_data(filepath):
    """
    Reads binary data written in the format:
    - n (number of vectors)
    - size of each vector
    - vector labels
    - optional independent variable (if included)
    - vectors (as floats).
    """
    with open(filepath, "rb") as file:
        # Read number of vectors
        n = struct.unpack('i', file.read(4))[0]

        # Read size of vectors
        size = struct.unpack('i', file.read(4))[0]

        # Read labels
        labels = []
        for _ in range(n):
            label_length = struct.unpack('i', file.read(4))[0]
            label = file.read(label_length).decode('utf-8')
            labels.append(label)

        # Read optional independent variable flag
        has_x_axis = struct.unpack('i', file.read(4))[0]
        print("has axis", has_x_axis)

        # Read the optional x-axis vector
        x_axis = None
        if has_x_axis:
            x_axis = np.frombuffer(file.read(size * 4), dtype=np.float32)

        # Read the vectors
        data = []
        for _ in range(n):
            data.append(np.frombuffer(file.read(size * 4), dtype=np.float32))

        return labels, np.array(data), x_axis


def plot_data(labels, data, x_axis=None):
    """
    Plots multiple vectors with x-axis labels reflecting the provided x_axis, while using indices for plotting.
    """
    plt.figure(figsize=(10, 6))

    # Use indices as the actual x-values for plotting
    indices = np.arange(data.shape[1])

    # Print provided x_axis for debugging
    if x_axis is not None:
        print(f"x_axis (labels): {x_axis}")

    # Plot each data vector using indices as x-values
    for i in range(data.shape[0]):
        plt.plot(indices, data[i], label=labels[i])

    # Set x-axis ticks to match indices, with x_axis values as labels
    if x_axis is not None:
        plt.xticks(indices, labels=[f"{x:.2f}" for x in x_axis], rotation=45)  # Format tick labels

    # Label the axes and add other plot details
    plt.ylabel('Value')
    plt.title('Data Plot')
    plt.legend()
    plt.grid(True)
    plt.show()


def PlotData():
    filepath = Path(__file__).parent.parent / "tmp.bin"

    if os.path.exists(filepath):
        labels, data, x_axis = read_binary_data(filepath)
        plot_data(labels, data, x_axis)
    else:
        print(f"File not found: {filepath}")


if __name__ == "__main__":
    PlotData()
