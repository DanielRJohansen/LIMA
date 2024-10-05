# Preinstall cleanup
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


# Postinstall cleanup
rm -rf pkg
rm -rf src
rm main.tar.gz
rm lima-*
