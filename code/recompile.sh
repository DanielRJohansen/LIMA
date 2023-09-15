#!/bin/bash


#This script can be called from anywhere so
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"/build
rm -rf ./*

log_file="$SCRIPT_DIR/build.log"
cmake ../ -Wno-dev > "$log_file" 2>&1

make > "$log_file" 2>&1
mv LIMA_TESTS/limatests ../
