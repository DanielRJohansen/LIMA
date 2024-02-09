#!/bin/bash

# This scripts installs all the dependencies LIMA needs: gcc, cuda, cmake, make,
# Then it installs itself in /opt/LIMA/
# Finally it executes 2 tests so ensure everything is working correctly

#if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

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

# Prepare the source code
install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"

echo "Using $program_dir as install directory"
sudo rm -rf "$program_dir"
sudo mkdir "$program_dir"

# copy everything from installdir to program_dir
sudo cp -r "$install_dir"/* "$program_dir"/
sudo chmod -R a+rwx $program_dir

# Build the public "lima" executable
cd "$program_dir"/build
cmake "$program_dir/code/LIMA_APP/"
if [ $? -ne 0 ]; then
    echo "CMake failed"
    exit 1
fi
sudo make install
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi
# Make this readable for all users
#chmod 777 /opt/LIMA -R
echo -e "\n\tLIMA client have been installed\n\n"


# Build LIMA once in ~/LIMA, to ensure everything works
userprogram_dir="$HOME/LIMA"
# Make this dir readable and writable for the user and the group, but only readable for others
mkdir -p "$userprogram_dir"
cp -r "$install_dir"/* "$userprogram_dir"/

cd "$userprogram_dir/build"
rm -rf ./*
cmake ../ 
if [ $? -ne 0 ]; then
    echo "CMake failed"
    exit 1
fi
make install -j
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi
#chmod 777 $userprogram_dir -R

echo -e "\n\tAll LIMA applications have been installed\n\n\n"

## -- INSTALL LIMA done  -- ##










# Run Self Test
# check cuda works
$userprogram_dir"/build/code/LIMA_ENGINE/engine_self_test"
if [ $? -ne 0 ]; then
    echo "engine_self_test failed"
    exit 1
fi

# Run small sim
cd "$userprogram_dir"
if [ "$1" != "-notest" ]; then
    #su -c "./selftest.sh" $SUDO_USER
    sims_dir=/$HOME/LIMA/simulations
    echo "Running self test in dir $sims_dir"

    mkdir -p "$sims_dir"

    cd /$HOME/LIMA
    git clone --quiet https://github.com/DanielRJohansen/LIMA_data 2>/dev/null

    cp -r ./LIMA_data/* $sims_dir/ #exclude .gitignore
    #rsync -q -av --exclude '.*' ./LIMA_data/ "$sims_dir/"  # Exclude hidden files/directories

#    chmod 777 /home/$SUDO_USER/LIMA -R


    cd "$sims_dir"/T4Lysozyme
    #cd "$sims_dir"/manyt4

    #lima mdrun # doesnt work, because this scrip has sudo, and the program must run as normal user
    #$SUDO_USER -u lima mdrun  # Doesnt work, because 2nd arg must be mdrun, otherwise the program doesnt know what to do
fi
