import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import matplotlib.pyplot as plt
import sys

def CompareRMSD(trr1Path, trr2Path, groPath):
    # Load two trajectories
    traj1 = mda.Universe(groPath, trr1Path)
    print("traj1")
    traj2 = mda.Universe(groPath, trr2Path)
    print("traj2")

    # Load reference structure from gro file
    ref_universe = mda.Universe(groPath)
    ref = ref_universe.select_atoms("backbone")  # Reference selection

    # Align trajectories to reference
    rmsd1 = RMSD(traj1, ref, select="backbone")
    rmsd1.run()
    print("rmsd1")
    rmsd2 = RMSD(traj2, ref, select="backbone")
    rmsd2.run()
    print("rmsd2")

    # Plot RMSD comparison
    plt.plot(rmsd1.results.rmsd[:, 1], rmsd1.results.rmsd[:, 2], label=trr1Path)
    plt.plot(rmsd2.results.rmsd[:, 1], rmsd2.results.rmsd[:, 2], label=trr2Path)
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (nm)")
    plt.legend()
    plt.title("RMSD Comparison")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python CompareRMSD.py <traj1.trr> <traj2.trr> <conf.gro>")
        sys.exit(1)

    CompareRMSD(sys.argv[1], sys.argv[2], sys.argv[3])