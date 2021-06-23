#!/bin/bash
#SBATCH --job-name=PFT
#SBATCH --time=0:30:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end
#SBATCH --mail-user=wwang138@jh.edu
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=48
#SBATCH --ntasks=36
#SBATCH --cpus-per-task=1
##SBATCH -o %A_%a.out

### exit on error
set -e 

mpiexec ./test.exe
wait
mpiexec python py-plot/plot_mpi.py
wait
date
