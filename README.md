# Phase Field Model


![tag](https://img.shields.io/badge/Version-v0.1-green.svg) ![gcc](https://img.shields.io/badge/GCC-9.3.0-lightgrey.svg) 
![openmpi](https://img.shields.io/badge/OpenMPI-3.1.6-red.svg) 
![python](https://img.shields.io/badge/Python-3.8.6-brightgreen.svg) ![hdf5](https://img.shields.io/badge/HDF5-1.10.7-ff69b4.svg) 
![cmake](https://img.shields.io/badge/CMake-3.18.4-eacd76.svg) ![license](https://img.shields.io/badge/License-MIT-bddd22.svg)

The resource management system for the HPC clusters is
![slurm](https://img.shields.io/badge/SLURM-19.05.7-44cef6.svg).

Python is too slow, so I choose C++!

## How to compile
1. `cd build`
2. `cmake ..`
3. `make`

If everything goes smoothly, after `cd ..`, you will see an executable file *test.exe* in the main directory.

## How to run
1. create a *data* directory `mkdir data`
2. submit your job to the cluster `sbatch script.slurm` 

## How to plot
1. creat a *figures* directory `mkdir figures`
2. `python py-plot/plot.py` or use the MPI version `plot_mpi.py`



<br />

>	Created by W. Wang on 2021-6-23.
>
>	Copyright © Wei Wang.