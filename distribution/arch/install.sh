#!/bin/bash

set -e

script_dir=$(dirname "$(realpath "$0")")
cd "$script_dir"

PushBuildToGit() {
    # Navigate to the Downloads/LIMAMD directory
    cd ~/Downloads/LIMAMD || { echo "Directory Downloads/LIMAMD not found."; exit 1; }

    # Check for changes in LIMAMD repository
    git_status=$(git status --porcelain)
    if [ -z "$git_status" ]; then
        echo "No changes detected in LIMAMD. Do you want to continue? (Yes will skip the git part)"
        while true; do
            read -p "Continue without pushing changes? (y/n): " yn
            case $yn in
                [Yy]* ) return;;
                [Nn]* ) echo "Aborting."; exit 1;;
                * ) echo "Please answer yes or no.";;
            esac
        done
    fi

    # Start by pushing the new build to LIMAMD if there are changes
    git status

    # Ask for confirmation before continuing
    while true; do
        read -p "Do you wish to push the update to LIMAMD? (y/n): " yn
        case $yn in
            [Yy]* ) break;;
            [Nn]* ) echo "Aborting."; exit 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done

    git add -A
    read -p "Enter commit message: " commit_msg
    git commit -m "$commit_msg"
    git push origin main
}






# Start by pushing the updated build
PushBuildToGit

cd "$script_dir"


# Update PKGBUILD
# Automatically increment the minor version number in the PKGBUILD
pkgver=$(grep "^pkgver=" PKGBUILD | cut -d= -f2)
major=$(echo $pkgver | cut -d. -f1)
minor=$(echo $pkgver | cut -d. -f2)
new_minor=$((minor + 1))
new_pkgver="$major.$new_minor"

# Update the pkgver in PKGBUILD
sed -i "s/pkgver=$pkgver/pkgver=$new_pkgver/" PKGBUILD

# Fetch and calculate new checksum
wget https://github.com/DanielRJohansen/LIMAMD/archive/main.tar.gz -O main.tar.gz
new_checksum=$(sha256sum main.tar.gz | awk '{ print $1 }')

# Update the checksum in PKGBUILD
sed -i "s/sha256sums=('.*')/sha256sums=('$new_checksum')/" PKGBUILD

# Cleanup
rm main.tar.gz

# Notify user
echo "PKGBUILD updated with new version $new_pkgver and checksum."






# Now test that we can install from the PKGBUILD

## Preinstall cleanup
sudo pacman -R lima
makepkg -Cf
sudo rm -rf /usr/share/LIMA
sudo rm /usr/bin/lima
rm -rf src
rm -rf pkg
rm lima*
rm main*
ls

makepkg -si

## Postinstall cleanup
rm -rf pkg
rm -rf src
rm main.tar.gz
rm lima-*
