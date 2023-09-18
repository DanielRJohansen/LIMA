#!/bin/bash

# This script is called after the file user constants file is made, so start with copying that
cp UserConstants.h ~/LIMA/source/LIMA_BASE/include/

mkdir -p ~/LIMA/source/build
cd ~/LIMA/source/build
rm -rf ./*


log_file=./limabuild.log
cmake ../ -Wno-dev > "$log_file" 2>&1
if [ $? -ne 0 ]; then
    echo "CMake failed"
    exit 1
fi

make install > "$log_file" 2>&1
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi

mv LIMA_TESTS/limatests ../
exit 0
