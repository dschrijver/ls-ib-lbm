#include <mpi.h>

#include "include/datatypes.h"
#include "include/memory.h"
#include "params.h"
#include "include/initialize.h"
#include "include/stencil.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    SimulationBag *sim;
    DistributionBag *dists;
    FieldBag *fields;
    ParamBag *params;
    Stencil *stencil;

    allocate_bags(&sim, &dists, &fields, &params, &stencil);

    set_params(params);

    initialize_MPI(params);

    initialize_stencil(sim);

    allocate_distributions(sim);

    allocate_fields(sim);

    params->t = 0;
    params->t_output = 0;
    params->t_log = 0;
    params->n_output = 0;

    MPI_Finalize();
    return 0;
}