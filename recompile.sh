#!/bin/bash


install_dir="$PWD"
program_dir="$PWD"/../../LIMA
apps_dir="$program_dir"/Applications
sims_dir="$program_dir"/Simulations""


mkdir -p "$apps_dir"
mkdir -p "$sims_dir"

mkdir "$apps_dir"/
mkdir "$apps_dir"/build

cp -r ./code/* "$apps_dir"/


cd "$apps_dir"/build
cmake ../ -Wno-dev



make
mv LIMA_TESTS/limatests ../

cd $apps_dir
