#!/bin/bash


#This script can be called from anywhere so
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# We are currently in opt which we dont have privileges to alter. So we move everything to local

mkdir -p ~/LIMA/source
mkdir -p ~/LIMA/applications
rm -rf ~/LIMA/source/*
#cp -rf "$SCRIPT_DIR"/* ~/LIMA/source
cp -rf /opt/LIMA/source/* ~/LIMA/source

exit 0
