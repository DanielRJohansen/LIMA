from UpgradeableFileFormat import *
import numpy as np
import matplotlib.pyplot as plt

def LoadData():
    data_path = r"C:\Users\Daniel\git_repo\LIMA_data\CoordCompression\h2o\trajectory.uff"
    file = UpgradeableFileParser(data_path)

    numAtoms = file.get_section("numAtoms", 'int32')[0]
    numFrames = file.get_section("numFrames", 'int32')[0]
    print("Num frames ", numFrames, " num atoms ", numAtoms)

    data = np.array(file.get_section("trajectory", 'float32'), dtype=np.float32)

    # Sanity check
    expected_size = numFrames * numAtoms * 3
    if data.size != expected_size:
        raise ValueError(f"Expected {expected_size} floats, got {data.size}")

    # Reshape: (numFrames *  numMolecules, 9)
    data_np = data.reshape(numFrames * (numAtoms // 3), 9)
    return data_np


def PlotHydrogenDensity(bins=128):
    data = LoadData()

    # data shape: [N, 9]
    # For each molecule: [ox, oy, oz, h1x, h1y, h1z, h2x, h2y, h2z]

    # Reshape
    oxygen = data[:, 0:3]
    h1 = data[:, 3:6]
    h2 = data[:, 6:9]

    # Compute relative positions
    rel_h1 = h1 - oxygen
    rel_h2 = h2 - oxygen
    rel_h = np.vstack((rel_h1, rel_h2))  # Shape: [2N, 3]

    # Choose projection plane (e.g. XY)
    x = rel_h[:, 0]
    y = rel_h[:, 1]

    # 2D histogram
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = np.log1p(heatmap)  # log(1 + count) to avoid log(0)

    # Plot
    plt.imshow(heatmap.T, origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               aspect='equal')
    plt.title("Log Density of Hydrogen Atoms Relative to Oxygen (XY plane)")
    plt.xlabel("X [relative]")
    plt.ylabel("Y [relative]")
    plt.colorbar(label='log(count + 1)')
    plt.show()



if __name__ == '__main__':
    PlotHydrogenDensity()