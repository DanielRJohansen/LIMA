#!/bin/bash

echo "\nWelcome to the LIMA Dynamics installer\n"

echo "Installing dependencies"
pacman -S cmake

pacman -S cuda --noconfirm
pacman -S cuda-tools --noconfirm












install_dir=$PWD
program_dir=~/Desktop/LIMA/
apps_dir="$program_dir"/Applications/

echo "Using $program_dir as install directory"
rm -rf "$program_dir"/
mkdir -p "$apps_dir"/dependencies



# Install glfw
pacman -S glfw-x11 --noconfirm
cp -r ./dependencies/* "$program_dir"/Dependencies/





mkdir -p "$program_dir/Applications"
mkdir -p "$program_dir/Simulation"



######### Now make Quantom

mkdir "$apps_dir"/src
mkdir "$apps_dir"/build

cp ./code/* "$apps_dir"/src/
mv "$apps_dir"/src/CMakeLists.txt "$apps_dir/"


cd "$apps_dir"/build
cmake ../
make
mv mdrun ../

printf "All LIMA applications have been installed\n"



## Run DEMO

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1

