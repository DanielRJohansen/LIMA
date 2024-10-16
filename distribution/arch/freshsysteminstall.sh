#!/bin/bash

# Variables
DOCKER_IMAGE="manjarolinux/base"
CONTAINER_NAME="manjaro_test_container"
PKGBUILD_PATH="./PKGBUILD"  # Update this path to the location of your PKGBUILD file
PKG_NAME="lima"            # Replace with the actual package name
USER_NAME="testuser"

# Check if Docker image exists
if ! docker image inspect "$DOCKER_IMAGE" > /dev/null 2>&1; then
    echo "Manjaro Docker image not found. Downloading..."
    docker pull "$DOCKER_IMAGE"
    if [ $? -ne 0 ]; then
        echo "Failed to pull the Manjaro Docker image."
        exit 1
    fi
fi

# Run Manjaro container
echo "Starting Manjaro container..."
docker run --gpus all --name "$CONTAINER_NAME" -itd "$DOCKER_IMAGE" /bin/bash

if [ $? -ne 0 ]; then
    echo "Failed to start the Manjaro container."
    exit 1
fi

# Update system
echo "Updating system and install base-devel..."
docker exec "$CONTAINER_NAME" bash -c "pacman -Syu --noconfirm && pacman -S --noconfirm base-devel"
if [ $? -ne 0 ]; then
    echo "Failed to install base-devel in the container."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Create a non-root user
echo "Creating non-root user..."
docker exec "$CONTAINER_NAME" bash -c "useradd -m $USER_NAME && echo '$USER_NAME ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers"
if [ $? -ne 0 ]; then
    echo "Failed to create a non-root user."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Copy PKGBUILD into the container
echo "Copying PKGBUILD into the container..."
docker cp "$PKGBUILD_PATH" "$CONTAINER_NAME":/home/$USER_NAME/PKGBUILD
if [ $? -ne 0 ]; then
    echo "Failed to copy the PKGBUILD file into the container."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Build and install the package as the non-root user
echo "Building and installing the package..."
docker exec -u "$USER_NAME" "$CONTAINER_NAME" bash -c "cd /home/$USER_NAME && makepkg -si --noconfirm"
if [ $? -ne 0 ]; then
    echo "Failed to build or install the package."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

# Verify package installation
echo "Verifying installation..."
docker exec "$CONTAINER_NAME" bash -c "pacman -Qs $PKG_NAME"
if [ $? -eq 0 ]; then
    echo "Package installed successfully."
else
    echo "Package installation failed."
    docker rm -f "$CONTAINER_NAME"
    exit 1
fi

echo "Exec: lima -help"
docker exec "$CONTAINER_NAME" bash -c "export XDG_RUNTIME_DIR=/tmp && lima -help"
echo "Exec: lima selftest"
docker exec "$CONTAINER_NAME" bash -c "export XDG_RUNTIME_DIR=/tmp && lima selftest"
docker rm -f "$CONTAINER_NAME"
exit 0
