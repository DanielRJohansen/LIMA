#!/bin/bash

if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi


echo "Welcome to the LIMA Dynamics installer"

echo "Installing dependencies"
#pacman -S cmake --noconfirm
#pacman -S make --noconfirm
#pacman -S cuda --noconfirm
#pacman -S cuda-tools --noconfirm












install_dir=$PWD
program_dir=/home/"$SUDO_USER"/Desktop/LIMA
apps_dir="$program_dir"/Applications

echo "Using $program_dir as install directory"
rm -rf "$program_dir"/


mkdir -p "$apps_dir"
mkdir -p "$program_dir/Simulation"
mkdir -p "$apps_dir"/dependencies



# Install glfw
#pacman -S glfw-x11 --noconfirm
cp -r ./dependencies/* "$apps_dir"/dependencies/







######### Now make Quantom

mkdir "$apps_dir"/src
mkdir "$apps_dir"/build

cp -r ./code/* "$apps_dir"/src/
mv "$apps_dir"/src/CMakeLists.txt "$apps_dir/"


cd "$apps_dir"/build
#cmake -DCMAKE_CUDA_FLAGS=”-arch=sm_89” ../
cmake ../
make
mv mdrun ../

printf "All LIMA applications have been installed"



## Run DEMO

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1

