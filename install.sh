#!/bin/bash



if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

echo "Welcome to the LIMA Dynamics installer"

install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"
source_dir="$program_dir"/source


echo "Using $program_dir as install directory"
rm -rf "$program_dir"/




echo "Installing dependencies"
#mkdir -p "$source_dir"/dependencies
#cp -r ./dependencies/* "$source_dir"/dependencies/

# Determine the distribution
if [ -f /etc/arch-release ]; then
    DISTRO="Arch"
elif [ -f /etc/lsb-release ]; then
    DISTRO="Ubuntu"
else
    echo "Unsupported distribution"
    exit 1
fi

# Check if we should install external dependencies
    # Check if the user provided exactly one argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <-none|-all> (install external dependencies)"
    exit 1
fi
if [ "$1" = "-all" ]; then
    echo "Welcome to the LIMA Dynamics installer"
    echo "Installing dependencies"
    #pacman -S cmake --noconfirm
    #pacman -S make --noconfirm
    #pacman -S cuda --noconfirm
    #pacman -S cuda-tools --noconfirm
    ##pacman -S glfw-x11 --noconfirm

    case $DISTRO in
    "Arch")
        sudo pacman -S cmake --noconfirm
        sudo pacman -S make --noconfirm
        sudo pacman -S cuda --noconfirm
        sudo pacman -S cuda-tools --noconfirm
        sudo pacman -S base-devel --noconfirm
        ;;
    "Ubuntu")
        sudo apt-get update
        
        sudo apt-get install -y make
        sudo apt-get install -y nvidia-cuda-toolkit
        sudo apt-get install -y build-essential

        # apt's version of cmake is behind, so overwrite it
        sudo apt-get install -y cmake
        
        #sudo apt remove cmake -y
        #sudo apt install python3-pip
        #pip install cmake --upgrade
        #echo 'export PATH="$PATH:$HOME/.local/bin"' >> $HOME/.bashrc
        #source $HOME/.bashrc
        ;;
    esac
elif [ "$1" = "-none" ]; then
    echo "No dependencies will be installed."
else
    echo "Usage: $0 <-none|-all>"
    exit 1
fi



# Prepare the source code
mkdir "$program_dir"
mkdir -p "$source_dir"/build
cp -r "$install_dir"/code/* "$source_dir"

cp -r "$install_dir"/resources "$program_dir/."

# Build the public "lima" executable
cd "$source_dir"/build
cmake "$source_dir"/LIMA_APP/
make install
echo -e "\n\tAll LIMA applications have been installed\n\n\n"





# Run Self Test
# check cuda works
$source_dir"/build/LIMA_ENGINE/engine_self_test"
# Run small sim
cd "$install_dir"
if [ "$1" != "-notest" ]; then
    su -c "./selftest.sh" $SUDO_USER
fi
