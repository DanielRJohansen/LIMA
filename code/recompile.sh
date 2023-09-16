#!/bin/bash


#This script can be called from anywhere so
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# We are currently in opt which we dont have privileges to alter. So we move everything to local

mkdir -p ~/LIMA/source
mkdir -p ~/LIMA/applications
rm -rf ~/LIMA/source/*
cp -rf "$SCRIPT_DIR"/* ~/LIMA/source

cd ~/LIMA/source/build
rm -rf ./*

log_file=./build.log
cmake ../ -Wno-dev #> "$log_file" 2>&1

make install #> "$log_file" 2>&1
mv LIMA_TESTS/limatests ../
