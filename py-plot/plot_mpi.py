import numpy as np 
import matplotlib.pyplot as plt
from many_fields import system
import h5py
from tqdm import tqdm
from tools import save_ani
from mpi4py import MPI
import os

def animation(tt, k_index, omega_index, file_name):
    image_list = []
    for t in tt:
        image_list.append("figures/k%d_o%d_t%d.png" % (k_index, omega_index, t))

    duration = 0.2
    fps = 5
    save_ani(image_list, file_name, duration, fps)
    for image in image_list:
        os.remove(image)



comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

nomegas = 6
nks = 6

ntasks = nomegas * nks
ntasks_per_cpu = int(np.ceil(float(ntasks) / size))


for task in range(ntasks_per_cpu):

    index = rank * ntasks_per_cpu + task
    k_index = int(index / nomegas)
    omega_index = index % nomegas
    
    file_name = "data/therm_k_%d_omega_%d.h5" % (k_index, omega_index)
    file = h5py.File(file_name, "r")
    sys = system()
    sys.conf = file["confinement"][:]

    tt = np.arange(0, 1001, 20)

    pbar = tqdm(total=102, desc="Processing")

    plt.figure(figsize=(6,8))
    for t in tt:
        index = 0
        plt.clf()
        for f in sys.sys:
            group = file["field_%d" % index]
            dataset = group["t_%g" % t]
            f.phi = dataset[:]
            
            f.subdomain = dataset.attrs["subdomain"]
            f.center = dataset.attrs["center"]

            f.v0 = dataset.attrs["velocity"]

            f.fulllattice()
            f.vv = np.linalg.norm(f.v0)
            f.polarity = np.arctan2(f.v0[1], f.v0[0])
            index += 1
        
        sys.showplot(t, "k%d_o%d_t%d" % (k_index, omega_index, t), therm=True)
        pbar.update(1)
    file.close()

    file_name = "data/simu_k_%d_omega_%d.h5" % (k_index, omega_index)
    file = h5py.File(file_name, "r")
    sys.confiment = file["confinement"][:]

    for t in tt:
        index = 0
        plt.clf()
        for f in sys.sys:
            group = file["field_%d" % index]
            dataset = group["t_%g" % t]
            f.phi = dataset[:]
            
            f.subdomain = dataset.attrs["subdomain"]
            f.center = dataset.attrs["center"]

            f.v0 = dataset.attrs["velocity"]

            f.fulllattice()
            f.vv = np.linalg.norm(f.v0)
            f.polarity = np.arctan2(f.v0[1], f.v0[0])
            index += 1
        
        sys.showplot(t+100, "k%d_o%d_t%d" % (k_index, omega_index, t+1000))
        pbar.update(1)
        
    file.close()
    
    tt = np.arange(0, 2001, 20)

    animation(tt, k_index, omega_index, "figures/k_%d_omega_%d" % (k_index, omega_index))



