# First cleanup from previous builds
rm -rf resources
rm lima license.md README.md
rm -rf build
sudo rm /usr/bin/lima
sudo rm -rf /usr/share/LIMA

rm debian/lima.substvars debian/debhelper-build-stamp
rm -rf debian/lima

set -e

# Now download files and move to debian dir
curl -L https://github.com/DanielRJohansen/LIMAMD/archive/main.tar.gz -o source.tar.gz
tar xzf source.tar.gz
mv ./LIMAMD-main/* .
rm -r LIMAMD-main
rm source.tar.gz

# Build the deb package
echo -e "\n\n### Building deb package"
debuild -us -uc
mkdir build
mv ../lima* ./build/
rm -rf resources
rm lima license.md README.md

rm debian/lima.substvars debian/debhelper-build-stamp
rm -rf debian/lima

# Try to install using the deb package
echo -e "\n\n### Installing from deb package"
#sudo apt install ./build/lima_*.deb
sudo dpkg -i ./build/lima_*.deb
sudo apt-get install -f


lima selftest
