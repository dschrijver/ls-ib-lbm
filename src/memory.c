#include <stdlib.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../include/memory.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil, LSMBag **lsm)
{
    *sim = (SimulationBag *)malloc(sizeof(SimulationBag));
    *dists = (DistributionBag *)malloc(sizeof(DistributionBag));
    *fields = (FieldBag *)malloc(sizeof(FieldBag));
    *params = (ParamBag *)malloc(sizeof(ParamBag));
    *stencil = (Stencil *)malloc(sizeof(Stencil));
    *lsm = (LSMBag *)malloc(sizeof(LSMBag));

    (*sim)->params = (*params);
    (*sim)->dists = (*dists);
    (*sim)->fields = (*fields);
    (*sim)->stencil = (*stencil);
    (*sim)->lsm = (*lsm);
}

void allocate_stencil(SimulationBag *sim)
{
    Stencil *stencil = sim->stencil;

    int NP = stencil->NP;

    stencil->cx = (int *)malloc(NP * sizeof(int));
    stencil->cy = (int *)malloc(NP * sizeof(int));
    stencil->cz = (int *)malloc(NP * sizeof(int));
    stencil->wp = (double *)malloc(NP * sizeof(double));
    stencil->p_bounceback = (int *)malloc(NP * sizeof(int));
}

void allocate_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;
    int malloc_size = (NX_proc + 2) * NY * NZ * NP * sizeof(double);

    dists->f1 = (double *)malloc(malloc_size);
    dists->f2 = (double *)malloc(malloc_size);
}

void allocate_fields(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;
    int malloc_size = (NX_proc + 4) * NY * NZ * sizeof(double);

    fields->rho = (double *)malloc(malloc_size);
    fields->pressure = (double *)malloc(malloc_size);
    
    fields->u = (double *)malloc(malloc_size);
    fields->v = (double *)malloc(malloc_size);
    fields->w = (double *)malloc(malloc_size);

    fields->Fx = (double *)malloc(malloc_size);
    fields->Fy = (double *)malloc(malloc_size);
    fields->Fz = (double *)malloc(malloc_size);

    fields->Fx_grav = (double *)malloc(malloc_size);
    fields->Fy_grav = (double *)malloc(malloc_size);
    fields->Fz_grav = (double *)malloc(malloc_size);
    fields->Fx_IBM = (double *)malloc(malloc_size);
    fields->Fy_IBM = (double *)malloc(malloc_size);
    fields->Fz_IBM = (double *)malloc(malloc_size);
    fields->Fx_rigid = (double *)malloc(malloc_size);
    fields->Fy_rigid = (double *)malloc(malloc_size);
    fields->Fz_rigid = (double *)malloc(malloc_size);

    memset(fields->Fx_IBM, 0, malloc_size);
    memset(fields->Fy_IBM, 0, malloc_size);
    memset(fields->Fz_IBM, 0, malloc_size);

    memset(fields->Fx_rigid, 0, malloc_size);
    memset(fields->Fy_rigid, 0, malloc_size);
    memset(fields->Fz_rigid, 0, malloc_size);
}