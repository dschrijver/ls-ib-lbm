#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/initialize.h"

void initialize_MPI(ParamBag *params)
{
    MPI_Comm_size(MPI_COMM_WORLD, &params->number_of_processes);
    int dims[3] = {params->number_of_processes, 1, 1};
    int xperiodic = 0;
#ifdef XPERIODIC
    xperiodic = 1;
#endif
    int periods[3] = {xperiodic, 0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &params->comm_xslices);
    MPI_Comm_rank(params->comm_xslices, &params->process_rank);
    MPI_Cart_get(params->comm_xslices, 3, dims, periods, params->process_coords);
    MPI_Cart_shift(params->comm_xslices, 0, 1, &params->process_neighbors[0], &params->process_neighbors[1]);
    params->i_start = (float)(params->process_coords[0]) / (float)params->number_of_processes * params->NX;
    params->i_end = (float)(params->process_coords[0] + 1) / (float)params->number_of_processes * params->NX;
    params->NX_proc = params->i_end - params->i_start;

    if (params->NX < 2 * params->number_of_processes)
    {
        if (params->process_rank == 0)
        {
            printf("NX (%d) must be at least twice as big as number of processes (%d)!\n", params->NX, params->number_of_processes);
            fflush(stdout);
        }

        MPI_Finalize();
        exit(1);
    }
}