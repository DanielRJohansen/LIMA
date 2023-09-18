sims_dir=~/LIMA/simulations
echo "Running self test in dir $sims_dir"

mkdir -p "$sims_dir"

cd ~/LIMA
git clone --quiet https://github.com/DanielRJohansen/LIMA_data 2>/dev/null

cp -r ./LIMA_data/* $sims_dir/ #exclude .gitignore
#rsync -q -av --exclude '.*' ./LIMA_data/ "$sims_dir/"  # Exclude hidden files/directories

#cd "$sims_dir"/T4Lysozyme
cd "$sims_dir"/manyt4

lima mdrun
