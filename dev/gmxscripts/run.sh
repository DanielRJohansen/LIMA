#!/bin/bash
set -e  # Exit on any error

# Prepare energy minimization
gmx grompp -f em.mdp -c membrane.gro -p membrane.top -o em.tpr
gmx mdrun -deffnm em

# Prepare simulation
gmx grompp -f run.mdp -c em.gro -p membrane.top -o run.tpr
gmx mdrun -deffnm run

echo "EM and simulation completed successfully."
