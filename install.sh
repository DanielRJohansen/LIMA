#!/bin/bash

# This scripts installs all the dependencies LIMA needs: gcc, cuda, cmake, make,
# Then it installs itself in /opt/LIMA/
# Finally it executes 2 tests so ensure everything is working correctly



if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

echo "\nWelcome to the LIMA Dynamics installer\n"

install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"

echo "Using $program_dir as install directory"
rm -rf "$program_dir"/



## -- INSTALL DEPENDENCIES  -- ##
echo "Installing dependencies"

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
    echo "Welcome to the LIMA Dynamics installer"
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
        sudo apt-get update
        sudo apt-get install -y make
        sudo apt-get install -y nvidia-cuda-toolkit
        sudo apt-get install -y build-essential
        sudo apt-get install -y gcc-13
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

# Prepare the source code
mkdir "$program_dir"
mkdir "$program_dir/build"
mkdir "$program_dir/source"
#mkdir -p "$source_dir"/build
cp -r "$install_dir"/code/* "$program_dir/source"
cp -r "$install_dir"/resources "$program_dir/."

# Build the public "lima" executable
cd "$program_dir"/build
cmake "$program_dir/source/LIMA_APP/"
make install
echo -e "\n\tLIMA client have been installed\n\n"


# Build LIMA once in /opt/, to ensure everything works
cd "$program_dir/build"
rm -rf ./*
cmake ../source/ 
if [ $? -ne 0 ]; then
    echo "CMake failed"
    exit 1
fi
make install 
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi

echo -e "\n\tAll LIMA applications have been installed\n\n\n"

## -- INSTALL LIMA done  -- ##










# Run Self Test
# check cuda works
$program_dir"/build/LIMA_ENGINE/engine_self_test"

# Run small sim
cd "$install_dir"
if [ "$1" != "-notest" ]; then
    su -c "./selftest.sh" $SUDO_USER
fi
