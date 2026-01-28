#!/usr/bin/env bash
set -euo pipefail

make

# test to make sure it runs 
./nbody 3 0.01 50 5 > test.tsv 2> test.time
python3 plot.py test.tsv test.pdf 1

# Solar short test
./nbody solar.tsv 200 2000 50 > solar_short.tsv 2> solar_short.time
python3 plot.py solar_short.tsv solar_short.pdf 10000

echo "Done. Created: test.tsv/test.pdf and solar_short.tsv/solar_short.pdf"
