#!/bin/bash

# This script is called after the file user constants file is made, so start with copying that
cp UserConstants.h ~/LIMA/source/LIMA_BASE/include/


program_dir="~/LIMA/"

mkdir -p ~/LIMA/build
cd ~/LIMA/build
rm -rf ./*


log_file=./limabuild.log
cmake "$program_dir/source" -Wno-dev > "$log_file" 2>&1
if [ $? -ne 0 ]; then
    echo "CMake failed"
    cat "$log_file"
    exit 1
fi

#cores=$(lscpu | grep "Core(s) per socket" | awk '{print $4}')
#make -j$cores

make install > "$log_file" 2>&1
if [ $? -ne 0 ]; then
    echo "Make failed"
    cat "$log_file"
    exit 1
fi

mv LIMA_TESTS/limatests ../
exit 0
