#!/bin/bash

echo "\nWelcome to the LIMA Dynamics installer\n"

echo "Installing dependencies"
pacman -S cmake

pacman -S cuda --noconfirm
pacman -S cuda-tools --noconfirm

pacman -S glfw-x11 --noconfirm














program_dir=~/Desktop/LIMA/
echo "Using $program_dir as install directory"
rm -rf "$program_dir"/


mkdir -p "$program_dir/Applications"
mkdir -p "$program_dir/Simulation"



######### Now make Quantom
apps_dir="$program_dir"/Applications/
mkdir -p "$apps_dir"/include
mkdir "$apps_dir"/src
mkdir "$apps_dir"/build

cp ./Quantom/Quantom/*.*h "$apps_dir"/include/
cp ./Quantom/Quantom/*.cu "$apps_dir"/src/
cp ./Quantom/Quantom/*.cpp "$apps_dir"/src/

cp ./LIMA/Applications/Quantom/* "$apps_dir"/



# Handle GLFW for Quantom. This is a quite bad way...
GLFW_CMAKER=~/Downloads/glfw-3.3.6/CMakeLists.txt
if [ ! -f "$GLFW_CMAKER" ]; then
	echo "glfw-3.3.6 not present in Downloads folder. Terminating install."
	exit 1
fi
mkdir -p "$program_dir"/Dependencies/GLFW/
cp -r ~/Downloads/glfw-3.3.6/* "$program_dir"/Dependencies/GLFW/


#exit 0

######### Now make FFM
FFM_dir="$program_dir"/Applications/Forcefieldmaker
mkdir -p "$FFM_dir"/include
mkdir "$FFM_dir"/src
mkdir "$FFM_dir"/build

cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.h "$FFM_dir"/include/
cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.cpp "$FFM_dir"/src/
cp ./LIMA/Applications/Forcefieldmaker/* "$FFM_dir"/				# CMakeLists.txt and build.sh


## Now make SimPostprocessor
SPP_dir="$program_dir"/Applications/SimPreprocessor
mkdir -p "$SPP_dir"/include
mkdir "$SPP_dir"/src
mkdir "$SPP_dir"/build

cp ./LIMA_services/LIMA_services/*.h "$SPP_dir"/include/
cp ./LIMA_services/LIMA_services/*.cpp "$SPP_dir"/src/


## Make test sim
s_dir="$program_dir"/Simulation
mkdir -p "$s_dir/"Molecule
mkdir "$s_dir"/Forcefield
#cp ~/Downloads/QnD/* "$s_dir"/Forcefield/	# nb and b ff
#cp ~/Downloads/QnD/*/* "$s_dir"/Molecule/	# conf and topol
cp ./Demo/6lzm/* "$s_dir"/Molecule/		# conf and topol"
cp ./Demo/CHARMM_FF/* "$s_dir"/Forcefield/	# nb and b ff



## Compile all applications
cd "$FFM_dir"
chmod +x build.sh
./build.sh
#./ffmrun




cd "$apps_dir"
chmod +x build.sh
chmod +x mdrun.sh
./build.sh

printf "\n\n\n\n\n\n\n\n\n\n"
printf "All LIMA applications have been installed\n"



## Run DEMO

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1

cd "$FFM_dir"
./ffmrun


cd "$apps_dir"
./mdrun.sh
