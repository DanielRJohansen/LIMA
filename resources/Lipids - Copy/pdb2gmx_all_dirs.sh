#!/bin/bash

# Loop through each directory in the current directory
for dir in */ ; do
    if [ -d "$dir" ]; then
        # Remove the trailing slash to get the directory name
        dirname=${dir%/}

        # Change into the directory
        cd "$dirname"

        # Execute the gmx command with the directory name
        gmx pdb2gmx -f "${dirname}.pdb" -o "${dirname}.gro" -ff Slipids_2020 -water spc
        echo gmx pdb2gmx -f "${dirname}.pdb" -o "${dirname}.gro" -ff Slipids_2020 -water spc
        # Change back to the parent directory
        cd ..

    fi
    break
done
