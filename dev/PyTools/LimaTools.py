from MakeLipids import MakeLipids
import struct
import matplotlib.pyplot as plt
import numpy as np
import CompareParametersWithGromacs
import Potentials
import os
import Plot2dVec
import Energies
import math
from scipy.special import erfc
import Plotdata
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




def count_lines_in_directory():
    directory = R"C:\Users\Daniel\git_repo\LIMA\code"

    def count_lines_in_file(file_path):
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
            return sum(1 for _ in file)

    total_lines = 0
    extensions = ['.cpp', '.h', '.cuh', '.cu']
    excluded_dirs = {'LIMA_DEPENDS', 'dependencies'}

    for root, dirs, files in os.walk(directory):
        # Remove excluded directories from the walk
        dirs[:] = [d for d in dirs if d not in excluded_dirs]
        
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                file_path = os.path.join(root, file)
                total_lines += count_lines_in_file(file_path)

    print(f"Total lines of code: {total_lines}")


def PlotErfcScalar():
    kappa = 2.5
    x = np.linspace(0.001, 3, 500)  # Avoid division by zero

    # Compute the scalar term
    erfc_term = erfc(x * kappa)
    # scalar = erfc_term / x ** 2
    scalar = erfc_term + (2.0 * kappa / np.sqrt(3.14)) * x * np.exp(-kappa * kappa * x ** 2)
    # scalar = erfc_term + (2.0 * kappa / np.sqrt(3.14)) * (np.exp(-kappa*kappa * x**2))
    # scalar = erfc(x*kappa) + 2*kappa/(np.sqrt(3.14)*x)*np.exp(-kappa * kappa * x**2) * x
    # scalar = erfc(x*kappa)

    # Plot the function
    plt.figure(figsize=(8, 5))
    plt.plot(x, scalar, label=r'$\mathrm{Scalar} = \frac{\mathrm{erfc}(|x| \cdot \kappa)}{x^2}$')
    plt.title('Scalar Derived from erfc(x)')
    plt.xlabel('x')
    plt.ylabel('Scalar')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":

    #print(math.erfc(3))

    Plotdata.PlotData()

    #count_lines_in_directory()
    exit(0)

    # Get the current folder path
    folder_path = R"C:\Users\Daniel\git_repo\LIMA\resources\Slipids"

    # Loop over all .itp files in the current folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".itp"):
            file_path = os.path.join(folder_path, filename)

            # Read the file
            with open(file_path, 'r') as file:
                lines = file.readlines()

            # Replace "lipid_section" with ";lipid_section"
            updated_lines = [line.replace('lipid_section', ';lipid_section') for line in lines]

            # Write the changes back to the file
            with open(file_path, 'w') as file:
                file.writelines(updated_lines)

    print("Replacement complete.")
    #Potentials.ShowLJ()
    #Energies.ShowEnergies()
    #MakeLipids

    #bins, counts = read_histogram_data('C:/Users/Daniel/git_repo/LIMA_data/psome2/histogram_data.bin')
    #plot_histogram(bins, counts)

    #Potentials.ShowLJ()
    #CompareParametersWithGromacs.CompareParameters("C:/Users/Daniel/git_repo/LIMA_data/Forcefieldtest")

    #pot_energy = [38976.5, 33643.6, 26578.7, 18757.6, 11260.7, 5123.36, 1193.4, 13.6142, 1746.96, 6154.03, 12626.1, 20269.2, 28027.6, 34829.7, 39736, 42068.7, 41505.7, 38124.7, 32392.8, 25101.5, 17258.2, 9946.05, 4175.13, 742.553, 122.428, 2400.4, 7261.83, 14035.2, 21785, 29440.7, 35945, 40399.3, 42188.5, 41065.4, 37185, 31083.6, 23603.7, 15778.5, 8688.96, 3314.26, 396.775, 339.488, 3150.32, 8441.02, 15480.8, 23297.3, 30810.9, 36983.7, 40963.1, 42199.5, 40522, 36162.5, 29722.9, 22092.9, 14326.4, 7495.93, 2545.15, 157.815, 663.688, 3992.9, 9685.6, 16955.4, 24798.3, 32131, 37940.5, 41424.5, 42101.7, 39878.5, 35062.2, 28317.9, 20577.1, 12909.2, 6373.15, 1871.82, 26.9351, 1093.33, 4923.71, 10989, 18451.5, 26280.3, 33394.2, 38810.4, 41781, 41895.5, 39138.2, 33890, 26875.7, 19064.1, 11534.3, 5326.36, 1297.69, 4.79429, 1626.23, 5938.05, 12344.7, 19961.2, 27735.5, 34593.9]
    #plt.plot(pot_energy)
    #plt.show()


   # plotStuff()
