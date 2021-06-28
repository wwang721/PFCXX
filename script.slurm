#!/bin/bash
#SBATCH --job-name=PFT
#SBATCH --time=1:00:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end
#SBATCH --mail-user=wwang138@jh.edu
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=48
#SBATCH --ntasks=11
#SBATCH --cpus-per-task=1
##SBATCH -o %A_%a.out

### for other partition, use `sacctmgr show qos format=name,priority` 
### to see all possible values
###SBATCH -p bigmem
###SBATCH --qos=qos_bigmem


### exit on error
set -e 

# mpiexec ./test.exe
# wait
mpiexec --mca mpi_warn_on_fork 0 python py-plot/plot_mpi.py 
wait
date
