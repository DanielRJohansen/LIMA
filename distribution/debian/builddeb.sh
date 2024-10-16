#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

# Variables
DOCKER_IMAGE="ubuntu:latest"
CONTAINER_NAME="ubuntu_deb_build"
DEBIAN_DIR="./debian"  # Update this to the location of your debian files

# Pull the latest Ubuntu Docker image if not already present
if ! docker image inspect "$DOCKER_IMAGE" > /dev/null 2>&1; then
    echo "Ubuntu Docker image not found. Downloading..."
    docker pull "$DOCKER_IMAGE"
fi

# Run Ubuntu container
echo "Starting Ubuntu container..."
docker run --name "$CONTAINER_NAME" -itd "$DOCKER_IMAGE" /bin/bash

# Install necessary dependencies in the container
echo "Installing dependencies..."
docker exec "$CONTAINER_NAME" bash -c "apt-get update &&
apt-get install -y devscripts build-essential curl debhelper-compat libtbb12 libgl1 libglu1-mesa libglx0"

# Copy Debian directory into the container
echo "Copying Debian directory into the container..."
docker cp "$DEBIAN_DIR" "$CONTAINER_NAME":/root/debian

# Download source files
echo "Downloading and preparing source files..."
docker exec "$CONTAINER_NAME" bash -c "\
    cd /root && \
    curl -L https://github.com/DanielRJohansen/LIMAMD/archive/main.tar.gz -o source.tar.gz && \
    tar xzf source.tar.gz && \
    rm source.tar.gz"

# Move debian folder into the root of the extracted source
echo "Moving Debian directory to the source root..."
docker exec "$CONTAINER_NAME" bash -c "mv /root/debian /root/LIMAMD-main/debian"

# Rename the directory to match the package name (lima)
echo "Renaming directory to match package name..."
docker exec "$CONTAINER_NAME" bash -c "mv /root/LIMAMD-main /root/lima"

# Build the deb package
echo "Building the .deb package..."
docker exec "$CONTAINER_NAME" bash -c "cd /root/lima && yes y | debuild -us -uc"

# Move the built package to a specific directory within the container
echo "Moving the built .deb package..."
docker exec "$CONTAINER_NAME" bash -c "mkdir -p /root/lima/build && mv ../lima* /root/lima/build/"

# Copy the built .deb package back to the host
echo "Copying the built .deb package to the host..."
mkdir -p build
docker cp "$CONTAINER_NAME":/root/lima/build/lima_*.deb ./build/

# Clean up unnecessary files
echo "Cleaning up unnecessary files..."
docker exec "$CONTAINER_NAME" bash -c "\
    rm -rf /root/lima/resources /root/lima/lima /root/lima/license.md /root/lima/README.md && \
    rm /root/lima/debian/lima.substvars /root/lima/debian/debhelper-build-stamp && \
    rm -rf /root/lima/debian/lima"

# Install the .deb package
echo "Installing the .deb package..."
docker exec "$CONTAINER_NAME" bash -c "dpkg -i /root/lima/build/lima_*.deb || apt-get install -f -y"

echo "Deb package installed successfully."

# Stop and remove the container after completion
docker rm -f "$CONTAINER_NAME"
