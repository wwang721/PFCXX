#include <iostream>
#include <ctime>
#include <mpi.h>
#include "field.hpp"


int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);

    time_t start,end;
    start=time(NULL);   

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nomegas = 6;
	float * omegas;
	omegas = new float [nomegas];
    for(int i=0; i<nomegas; i++)
        omegas[i] = 0.01 * i;

	int nks = 6;
    float * ks;
	ks = new float [nks];
    for(int i=0; i<nks; i++)
        ks[i] = 0.002 * i;

    int ntasks = nomegas * nks;

    int ntasks_per_cpu = ceil(double(ntasks) / size);

    for(int task=0; task<ntasks_per_cpu; task++)
    {   
        int index = ntasks_per_cpu * rank + task;

        if(index >= ntasks)
            break;

        int k_index = index / nomegas;
        int omega_index = index % nomegas;
		
        System cells;
        
		char fileName[50];
        sprintf(fileName, "data/therm_k_%d_omega_%d", k_index, omega_index);
        
        cells.simulation(fileName, 1000, 0.5, true, ks[k_index], omegas[omega_index]);  // thermalization
		
        sprintf(fileName, "data/simu_k_%d_omega_%d", k_index, omega_index);
        cells.simulation(fileName, 1000, 0.5, false, ks[k_index], omegas[omega_index]);
    }

    // can't delete these when using MPI
    //delete[] omegas;
    //delete[] ks;

    // time
    end = time(NULL);

    if(rank == 0)
    { 
        int T,second,minute,hour;
        T = end-start;
        second = T%60;
        minute = T/60;
        hour = minute/60;
        minute = minute%60; 
        std::cout<<"\nRun Time: "<<hour<<" h "<<minute<<" min "<<second<<" s.\n"; 
    }
	MPI_Finalize();

    return 0;
}
