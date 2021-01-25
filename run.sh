#!/bin/bash
#SBATCH -n 1
#SBATCH -t 00:01:00
#SBATCH -J hsmc
#SBATCH -A snic2019-7-45

export OMP_NUM_THREADS=1
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=2
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=4
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=8
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=16
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=32
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=48
./lmp_corr -i "test/*.dat.gz"

export OMP_NUM_THREADS=64
./lmp_corr -i "test/*.dat.gz"
