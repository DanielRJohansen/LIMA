import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import os
import sys

def plot_probability_density_function(filepath, num_bins=100):
    # Load data from binary file
    data = np.fromfile(filepath, dtype=np.float32)

    data_min = min(data)
    data_max = max(data)
    print(data_min, data_max)
    plt.hist(data, bins=100, range=(data_min, data_max), density=True)  # setting the range
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.title("PDF of Data (Unnormalized)")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filepath = sys.argv[1]

        if os.path.exists(filepath):
            plot_probability_density_function(filepath)

        else:
            print(f"File not found: {filepath}")