#!/bin/bash

# DO NOT USE THIS SCRIPT
# It is for developers to find compilebugs only, please use install.sh



if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

echo "Welcome to the LIMA Dynamics installer"

install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"
apps_dir="$program_dir"/Applications
sims_dir="$program_dir"/Simulations


echo "Using $program_dir as install directory"
rm -rf "$program_dir"/




echo "Installing dependencies"
mkdir -p "$apps_dir"/dependencies
cp -r ./dependencies/* "$apps_dir"/dependencies/



# Check if we should install external dependencies
    # Check if the user provided exactly one argument
if [ "$#" -ne 1 ]; then
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
mkdir -p "$sims_dir"
mkdir -p "$apps_dir"
mkdir "$apps_dir"/build
cp -r ./code/* "$apps_dir"/



# Build the public "lima" executable
cd "$apps_dir"/build
cmake "$apps_dir"/LIMA_APP/
make install
#mv lima "$program_dir"



#mv "$apps_dir"/src/CMakeLists.txt "$apps_dir/"










#cmake -DCMAKE_CUDA_FLAGS=”-arch=sm_89” ../
#export CC=/opt/cuda/bin/gcc
#export CXX=/opt/cuda/bin/g++
cmake ../ 

printf "Make the self-test files \n"
cd LIMA_ENGINE
make
./engine_self_test
cd ..

cd LIMA_FORCEFIELDMAKER
make
./ffm_self_test
cd ..

cd LIMA_MD
make
./md_self_test
cd ..



#make
#mv LIMA_TESTS/limatests ../
#cp -r "$install_dir"/../LIMA_data/* "$sims_dir"/.



cd $apps_dir
#./limatests

#mv mdrun ../

printf "All LIMA applications have been installed"



## Run DEMO

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1
