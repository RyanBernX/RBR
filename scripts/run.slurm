#!/bin/bash
#SBATCH -J RBR
#SBATCH -o rbr.log
#SBATCH -t 2-1:00:00
#SBATCH -c 24


export OMP_NUM_THREADS=12
#./test_gset.rb
./test_dimacs10.rb
