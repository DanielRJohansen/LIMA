How to update the distribution

# Compile LIMA
First compile LIMA using install.sh in the toplevel dir ../../install.sh.
This updates the executable and resources dir in the Downloads/LIMAMD repo
commit and push that update

# Update PKGBUILD
Update the SHA checksum in PKGBUILD. Use the cmd below to get the new checksum:
wget https://github.com/DanielRJohansen/LIMAMD/archive/main.tar.gz && sha256sum main.tar.gz
Update version numbers

# Verify that we can correctly install using the PKGBUILD
./install.sh


# Publish PKGBUILD file to AUR
