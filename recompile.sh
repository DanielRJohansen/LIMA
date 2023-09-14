#!/bin/bash


install_dir="$PWD"
#program_dir=/home/"$SUDO_USER"/Desktop/LIMA
program_dir="$PWD"/../../LIMA
apps_dir="$program_dir"/Applications
sims_dir="$program_dir"/Simulations""


#echo "Using $program_dir as install directory"
#rm -rf "$program_dir"/


mkdir -p "$apps_dir"
mkdir -p "$sims_dir"
#mkdir -p "$apps_dir"/dependencies



# Install glfw
#pacman -S glfw-x11 --noconfirm
#cp -r ./dependencies/* "$apps_dir"/dependencies/







######### Now make Quantom

mkdir "$apps_dir"/
mkdir "$apps_dir"/build

cp -r ./code/* "$apps_dir"/
#mv "$apps_dir"/src/CMakeLists.txt "$apps_dir/"


cd "$apps_dir"/build
#cmake -DCMAKE_CUDA_FLAGS=”-arch=sm_89” ../
#export CC=/opt/cuda/bin/gcc
#export CXX=/opt/cuda/bin/g++
cmake ../ -Wno-dev



make
mv LIMA_TESTS/limatests ../
cp -r "$install_dir"/../LIMA_data/* "$sims_dir"/.

cd $apps_dir
