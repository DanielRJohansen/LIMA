import os
import subprocess
import re
from typing import List



class Bond:
    def __init__(self, b0: float, kB: float):
        self.b0 = b0
        self.kB = kB
class Angle:
    def __init__(self, theta0: float, kTheta: float):
        self.theta0 = theta0
        self.kTheta = kTheta
class Dihedral:
    def __init__(self, phi0: float, kPhi: float, multiplicity: int):
        self.phi0 = phi0
        self.kPhi = kPhi
        self.multiplicity = multiplicity
class ImproperDihedral:
    def __init__(self, psi0: float, kPsi: float):
        self.psi0 = psi0
        self.kPsi = kPsi






def compare_bonds(bonds1: List[Bond], bonds2: List[Bond], tol=0.01) -> bool:
    if len(bonds1) != len(bonds2):
        return False
    for b1, b2 in zip(bonds1, bonds2):
        if abs(b1.b0 - b2.b0) > tol or abs(b1.kB - b2.kB) > tol:
            return False
    return True

def compare_dihedrals(dihedrals1: List[Dihedral], dihedrals2: List[Dihedral], tol=0.01) -> bool:
    if len(dihedrals1) != len(dihedrals2):
        return False
    for d1, d2 in zip(dihedrals1, dihedrals2):
        if abs(d1.phi0 - d2.phi0) > tol or abs(d1.kPhi - d2.kPhi) > tol or d1.multiplicity != d2.multiplicity:
            return False
    return True

def compare_angles(angles1: List[Angle], angles2: List[Angle], tol=0.01) -> bool:
    if len(angles1) != len(angles2):
        return False
    for a1, a2 in zip(angles1, angles2):
        if abs(a1.theta0 - a2.theta0) > tol or abs(a1.kTheta - a2.kTheta) > tol:
            return False
    return True

def compare_improper_dihedrals(improper_dihedrals1: List[ImproperDihedral], improper_dihedrals2: List[ImproperDihedral], tol=0.01) -> bool:
    if len(improper_dihedrals1) != len(improper_dihedrals2):
        return False
    for id1, id2 in zip(improper_dihedrals1, improper_dihedrals2):
        if abs(id1.psi0 - id2.psi0) > tol or abs(id1.kPsi - id2.kPsi) > tol:
            return False
    return True





def parseForcefieldFromTpr(file_path: str):
    print("There")
    bonds = []
    angles = []
    dihedrals = []
    improper_dihedrals = []

    with open(file_path, 'r') as file:
        for line in file:
            if 'functype' in line and 'BONDS' in line:
                bond_match = re.search(r'b0A\s*=\s*([\d\.\+e\-]+),\s*cbA\s*=\s*([\d\.\+e\-]+)', line)
                if bond_match:
                    b0, kB = map(float, bond_match.groups())
                    bonds.append(Bond(b0, kB))
            elif 'functype' in line and 'ANGLES' in line:
                angle_match = re.search(r't0A\s*=\s*([\d\.\+e\-]+),\s*ctA\s*=\s*([\d\.\+e\-]+)', line)
                if angle_match:
                    theta0, kTheta = map(float, angle_match.groups())
                    angles.append(Angle(theta0, kTheta))
            elif 'functype' in line and 'PDIHS' in line:
                dihedral_match = re.search(r'phiA\s*=\s*([\d\.\+e\-]+),\s*cpA\s*=\s*([\d\.\+e\-]+),\s*phiB\s*=\s*([\d\.\+e\-]+),\s*cpB\s*=\s*([\d\.\+e\-]+),\s*mult\s*=\s*(\d+)', line)
                if dihedral_match:
                    phiA, cpA, phiB, cpB, multiplicity = dihedral_match.groups()
                    dihedrals.append(Dihedral(float(phiA), float(cpA), int(multiplicity)))
            elif 'functype' in line and 'IDIHS' in line:
                improper_dihedral_match = re.search(r'psi0\s*=\s*([\d\.\+e\-]+),\s*kPsi\s*=\s*([\d\.\+e\-]+)', line)
                if improper_dihedral_match:
                    psi0, kPsi = map(float, improper_dihedral_match.groups())
                    improper_dihedrals.append(ImproperDihedral(psi0, kPsi))

    return bonds, dihedrals, angles, improper_dihedrals

def ParseForcefieldFromItp(file_path: str):
    bonds = []
    angles = []
    dihedrals = []
    improper_dihedrals = []


    return bonds, angles, dihedrals, improper_dihedrals


def CompareAllSubdirs():
    # Get all subdirectories in the current directory
    subdirs = [d for d in os.listdir() if os.path.isdir(d)]

    # Loop through each subdir and run the commands
    for subdir in subdirs:
        os.chdir(subdir)

        # First get bondparameters chosen by GROMACS
        subprocess.run(['gmx', 'grompp', '-f', 'sim.mdp', '-c', 'molecule/conf.gro', '-p', 'molecule/topol.top', '-o', 'output.tpr'], check=True)
        subprocess.run(['gmx', 'dump', '-s', 'output.tpr', '>', 'tpr_content.txt'], shell=True, check=True)
        bondsG, anglesG, dihedralsG, impropersG = parseForcefieldFromTpr("tpr_content.txt")

        # Now get bondparameters from LIMA
        subprocess.run(['lima', 'getforcefieldparams'])
        bondsL, anglesL, dihedralsL, impropersL = ParseForcefieldFromItp("appliedForcefield.itp")

        compare_bonds(bondsL, bondsG)
        compare_angles(anglesL, anglesG)
        compare_dihedrals(dihedralsL, dihedralsG)
        compare_improper_dihedrals(impropersL, impropersG)

        os.chdir('..')







if __name__ == '__main__':


    exit()






