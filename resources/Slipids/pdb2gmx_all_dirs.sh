#!/bin/bash

# Loop through each .pdb file in the current directory
for pdb_file in *.pdb; do
  # Extract the lipid name without the extension
  lipidname=$(basename "$pdb_file" .pdb)

  # Convert the pdb file to gro format using editconf
  gmx editconf -f "${lipidname}.pdb" -o "${lipidname}.gro"
done
