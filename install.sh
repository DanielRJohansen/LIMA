#!/bin/bash



if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

echo "Welcome to the LIMA Dynamics installer"

install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"
apps_dir="$program_dir"/Applications


echo "Using $program_dir as install directory"
rm -rf "$program_dir"/




echo "Installing dependencies"
mkdir -p "$apps_dir"/dependencies
cp -r ./dependencies/* "$apps_dir"/dependencies/



# Check if we should install external dependencies
    # Check if the user provided exactly one argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <-none|-all> (install external dependencies)"
    exit 1
fi
if [ "$1" = "-all" ]; then
    echo "Welcome to the LIMA Dynamics installer"
    echo "Installing dependencies"
    pacman -S cmake --noconfirm
    pacman -S make --noconfirm
    pacman -S cuda --noconfirm
    pacman -S cuda-tools --noconfirm
    #pacman -S glfw-x11 --noconfirm
elif [ "$1" = "-none" ]; then
    echo "No dependencies will be installed."
else
    echo "Usage: $0 <-none|-all>"
    exit 1
fi











#if(UNIX)
    #set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
#fi()








# Prepare the source code
mkdir -p "$apps_dir"
mkdir "$apps_dir"/build
cp -r ./code/* "$apps_dir"/



# Build the public "lima" executable
cd "$apps_dir"/build
cmake "$apps_dir"/LIMA_APP/
make install

# Add the following line to add lima to the PATH
#echo 'export PATH="$program_dir/bin:$PATH"' >> ~/.bashrc
#source ~/.bashrc  # Reload the Bash configuration

echo "Current PATH: $PATH"
echo 'export PATH="/opt/LIMA/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
echo "Updated PATH: $PATH"



echo "All LIMA applications have been installed"




lima
exit 0

# Run Self Test

if [ "$1" != "-notest" ]; then
    printf "Running self test"

    sims_dir="$program_dir"/Simulations
    mkdir -p "$sims_dir"

    git clone https://github.com/DanielRJohansen/LIMA_data

    #cp -r --exclude '.*' ./LIMA_data/* $sims_dir/ #exclude .gitignore
    rsync -q -av --exclude '.*' ./LIMA_data/ "$sims_dir/"  # Exclude hidden files/directories

    cd "$sims_dir"
    lima
fi

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1
