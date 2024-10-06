#!/bin/bash

# This scripts installs all the dependencies LIMA needs: gcc, cuda, cmake, make,
# Then it installs itself in /opt/LIMA/
# Finally it executes 2 tests so ensure everything is working correctly

#if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

set -e

echo "\nWelcome to the LIMA Dynamics installer\n"


## -- INSTALL DEPENDENCIES  -- ##

# Determine the distribution
if [ -f /etc/arch-release ]; then
    DISTRO="Arch"
elif [ -f /etc/lsb-release ]; then
    DISTRO="Ubuntu"
else
    echo "Unsupported distribution"
    exit 1
fi

# Check if we should install external dependencies
    # Check if the user provided exactly one argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <-none|-all> (install external dependencies)"
    exit 1
fi
if [ "$1" = "-all" ]; then
    echo "Installing dependencies"

    case $DISTRO in
    "Arch")
        sudo pacman -S cmake --noconfirm
        sudo pacman -S make --noconfirm
        sudo pacman -S cuda --noconfirm
        sudo pacman -S cuda-tools --noconfirm
        sudo pacman -S base-devel --noconfirm
        ;;
    "Ubuntu")
        sudo apt-get install -y make
        sudo apt-get install -y nvidia-cuda-toolkit
        sudo apt-get install -y build-essential
        sudo apt-get install -y gcc-13 g++-13
        sudo apt-get install -y cmake
        ;;
    esac
elif [ "$1" = "-none" ]; then
    echo "No dependencies will be installed."
else
    echo "Usage: $0 <-none|-all>"
    exit 1
fi
## -- INSTALL DEPENDENCIES done  -- ##






## -- INSTALL LIMA  -- ##

LIMA_DIR="/usr/share/LIMA"    # Used by Cmake to save all resource files
sudo rm -rf $LIMA_DIR

BUILD_DIR="$(pwd)/build"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake ../ 
if [ $? -ne 0 ]; then
    echo "CMake failed"
    exit 1
fi
make -j
sudo make install -j
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi
#rm -rf "$BUILD_DIR"/*


echo -e "\n\t LIMA has been installed successfully \n\n\n"
## -- INSTALL LIMA done  -- ##



# Move files to LIMAMD, prepared for distribution
rm ~/Downloads/LIMAMD/lima
cp /usr/bin/lima ~/Downloads/LIMAMD/
rm -rf ~/Downloads/LIMAMD/resources
cp -r /usr/share/LIMA/resources ~/Downloads/LIMAMD

lima selftest


# Run small sim
if [ "$1" != "-notest" ]; then
    SIMS_DIR=/$HOME/LIMA/simulations
    echo "Running self test in dir $SIMS_DIR"
    mkdir -p "$SIMS_DIR"

    # Clone the data repository only if the directory doesn't exist
    cd /$HOME/Downloads
    [ ! -d "LIMA_data" ] && git clone --quiet https://github.com/DanielRJohansen/LIMA_data 2>/dev/null

    cp -r ./LIMA_data/* $SIMS_DIR/ #exclude .gitignore


    cd "$SIMS_DIR/T4Lysozyme"
    lima makesimparams
    sed -i 's/n_steps=[0-9]*/n_steps=50000/' sim_params.txt
    lima mdrun
    #cd "$sims_dir"/manyt4

    #lima mdrun # doesnt work, because this scrip has sudo, and the program must run as normal user
    #$SUDO_USER -u lima mdrun  # Doesnt work, because 2nd arg must be mdrun, otherwise the program doesnt know what to do
fi
