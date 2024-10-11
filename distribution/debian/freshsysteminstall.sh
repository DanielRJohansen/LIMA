#!/bin/bash

# Variables
DOCKER_IMAGE="ubuntu:23.04"
CONTAINER_NAME="fresh_ubuntu_install_container"
DEB_PACKAGE_PATH="./build/lima_*.deb"  # Update this path to the location of your .deb package
PKG_NAME="lima"

# Check if Docker image exists
if ! docker image inspect "$DOCKER_IMAGE" > /dev/null 2>&1; then
    echo "Docker image $DOCKER_IMAGE not found. Downloading..."
    docker pull "$DOCKER_IMAGE"
    if [ $? -ne 0 ]; then
        echo "Failed to pull the Docker image."
        exit 1
    fi
fi

# Run Docker container
echo "Starting Docker container..."
docker run --gpus all --name "$CONTAINER_NAME" -itd "$DOCKER_IMAGE" /bin/bash

if [ $? -ne 0 ]; then
    echo "Failed to start the Docker container."
    exit 1
fi

# Update package lists
echo "Updating package lists..."
docker exec "$CONTAINER_NAME" bash -c "apt-get update"
if [ $? -ne 0 ]; then
    echo "Failed to update package lists."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Install necessary packages
echo "Installing necessary packages..."
docker exec "$CONTAINER_NAME" bash -c "apt-get install -y wget gnupg software-properties-common"
if [ $? -ne 0 ]; then
    echo "Failed to install necessary packages."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Add NVIDIA CUDA repository
echo "Adding NVIDIA CUDA repository..."
docker exec "$CONTAINER_NAME" bash -c "wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin"
docker exec "$CONTAINER_NAME" bash -c "mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600"
docker exec "$CONTAINER_NAME" bash -c "wget https://developer.download.nvidia.com/compute/cuda/12.2.1/local_installers/cuda-repo-ubuntu2204-12-2-local_12.2.1-535.86.10-1_amd64.deb"
docker exec "$CONTAINER_NAME" bash -c "dpkg -i cuda-repo-ubuntu2204-12-2-local_12.2.1-535.86.10-1_amd64.deb"
docker exec "$CONTAINER_NAME" bash -c "cp /var/cuda-repo-ubuntu2204-12-2-local/cuda-*-keyring.gpg /usr/share/keyrings/"

# Update package lists again
echo "Updating package lists after adding CUDA repository..."
docker exec "$CONTAINER_NAME" bash -c "apt-get update"
if [ $? -ne 0 ]; then
    echo "Failed to update package lists after adding CUDA repository."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Install CUDA Toolkit
echo "Installing CUDA Toolkit..."
docker exec "$CONTAINER_NAME" bash -c "apt-get install -y cuda"
if [ $? -ne 0 ]; then
    echo "Failed to install CUDA Toolkit."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Copy .deb package into the container
echo "Copying .deb package into the container..."
docker cp $DEB_PACKAGE_PATH "$CONTAINER_NAME":/tmp/
if [ $? -ne 0 ]; then
    echo "Failed to copy the .deb package into the container."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Install the .deb package
echo "Installing the .deb package..."
docker exec "$CONTAINER_NAME" bash -c "apt-get install -y ./tmp/lima_*.deb"
if [ $? -ne 0 ]; then
    echo "Failed to install the .deb package."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Verify package installation
echo "Verifying installation..."
docker exec "$CONTAINER_NAME" bash -c "dpkg -l | grep $PKG_NAME"
if [ $? -eq 0 ]; then
    echo "Package installed successfully."
else
    echo "Package installation failed."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Execute lima -help and lima selftest
echo "Exec: lima -help"
docker exec "$CONTAINER_NAME" bash -c "lima -help"
if [ $? -ne 0 ]; then
    echo "Failed to execute 'lima -help'."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

echo "Exec: lima selftest"
docker exec "$CONTAINER_NAME" bash -c "lima selftest"
if [ $? -ne 0 ]; then
    echo "Failed to execute 'lima selftest'."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Clean up
docker rm -f "$CONTAINER_NAME"
echo "Test completed successfully."
exit 0
