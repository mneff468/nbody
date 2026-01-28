#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --partition=Centaurus
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

echo "Node: $(hostname)"
echo "Dir:  $(pwd)"
ls -la

make clean
make

# benchmark cases
echo "Running solar..."
./nbody solar.tsv 200 5000000 5000 > results_solar.tsv 2> time_solar.txt

echo "Running N=100..."
./nbody 100 1 10000 100 > results_100.tsv 2> time_100.txt

echo "Running N=1000..."
./nbody 1000 1 10000 100 > results_1000.tsv 2> time_1000.txt

# Plot one
python3 plot.py results_100.tsv plot_100.pdf 1

echo "DONE"
ls -l results_*.tsv time_*.txt plot_100.pdf slurm-*.out slurm-*.err
