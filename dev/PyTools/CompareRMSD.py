import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import matplotlib.pyplot as plt
import sys

def CompareRMSD(trr1Path, trr2Path, tprPath):
    # Load two trajectories
    ref = mda.Universe("reference.pdb")  # Reference structure
    traj1 = mda.Universe(tprPath, trr1Path)
    traj2 = mda.Universe(tprPath, trr2Path)

    # Align trajectories to reference
    rmsd1 = RMSD(traj1, ref, select="backbone")
    rmsd2 = RMSD(traj2, ref, select="backbone")
    rmsd1.run()
    rmsd2.run()

    # Plot RMSD comparison
    plt.plot(rmsd1.rmsd[:, 1], rmsd1.rmsd[:, 2], label="Run 1")
    plt.plot(rmsd2.rmsd[:, 1], rmsd2.rmsd[:, 2], label="Run 2")
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (nm)")
    plt.legend()
    plt.title("RMSD Comparison")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python CompareRMSD.py <traj1.trr> <traj2.trr>")
        sys.exit(1)

    CompareRMSD(sys.argv[1], sys.argv[2])