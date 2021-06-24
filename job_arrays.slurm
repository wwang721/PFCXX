#!/bin/bash
#SBATCH --job-name=PFT
#SBATCH --time=0:30:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end
#SBATCH --mail-user=wwang138@jh.edu
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=48
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH -o %A_%a.out
#SBATCH --array=0-47

### exit on error
set -e 

./test.exe $SLURM_ARRAY_TASK_ID
wait
##mpiexec python py-plot/plot_mpi.py 
##wait
date
