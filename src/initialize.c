#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/utils.h"
#include "../include/lsm.h"
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

void initialize_fields(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;
    Stencil *stencil = sim->stencil;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double rho_0 = params->rho_0;

    double cs2 = stencil->cs2;

    double *rho = fields->rho;
    double *pressure = fields->pressure;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

#ifdef INI_POISEUILLE
    FOR_DOMAIN
    {
        u[INDEX(i, j, k)] = 0.0;
        v[INDEX(i, j, k)] = 0.0;
        w[INDEX(i, j, k)] = 0.0;
        rho[INDEX(i, j, k)] = rho_0;
    }
#endif

    FOR_DOMAIN
    {
        pressure[INDEX(i, j, k)] = rho[INDEX(i, j, k)] * cs2;
    }
}

void initialize_distribution(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    FieldBag *fields = sim->fields;
    Stencil *stencil = sim->stencil;

    double rho_i, uhat, vhat, what, uc, u2;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *rho = fields->rho;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *f1 = dists->f1;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX(i, j, k)];
        uhat = u[INDEX(i, j, k)] - 1.0 / (2.0 * rho_i) * Fx[INDEX(i, j, k)];
        vhat = v[INDEX(i, j, k)] - 1.0 / (2.0 * rho_i) * Fy[INDEX(i, j, k)];
        what = w[INDEX(i, j, k)] - 1.0 / (2.0 * rho_i) * Fz[INDEX(i, j, k)];
        u2 = uhat * uhat + vhat * vhat + what * what;

        for (int p = 0; p < NP; p++)
        {
            uc = uhat * (double)cx[p] + vhat * (double)cy[p] + what * (double)cz[p];
            f1[INDEX_F(i, j, k, p)] = rho_i * wp[p] * (uc / cs2 + uc * uc / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
        }
    }
}

void initialize_lsm(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    LSMBag *lsm = sim->lsm;

    double x, y, z;

    int NY = params->NY;
    int NZ = params->NZ;

    lsm->n_particles = 0;
    lsm->n_springs = 0;
    lsm->particles = NULL;
    lsm->springs = NULL;

#ifdef INI_POISEUILLE
    int x_left = params->NX / 2 - WIDTH / 2;
    for (int i = x_left; i < x_left + WIDTH; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                x = (double)i;
#ifndef LEFT_NEBB_VELOCITY
                x += 0.5;
#endif
                y = (double)j;
#ifndef BOTTOM_NEBB_VELOCITY
                y += 0.5;
#endif
                z = (double)k;
#ifndef BACK_NEBB_VELOCITY
                z += 0.5;
#endif

                create_particle(x, y, z, 0.0, 0.0, 0.0, lsm);
            }
        }
    }
#endif

    add_springs(sim);

    int n_particles = lsm->n_particles;
    int n_springs = lsm->n_springs;

    if (n_particles < params->number_of_processes)
    {
        if (params->process_rank == 0)
        {
            printf("The number of processes (%d) must be smaller than the number of LSM particles (%d)!\n", params->number_of_processes, n_particles);
            fflush(stdout);
        }

        MPI_Finalize();
        exit(1);
    }

    compute_chi(sim);

    // MPI for LSM
    lsm->p_start = (float)(params->process_rank) / (float)params->number_of_processes * n_particles;
    lsm->p_end = (float)(params->process_rank + 1) / (float)params->number_of_processes * n_particles;
    lsm->p_proc = lsm->p_end - lsm->p_start;

    lsm->s_start = (float)(params->process_rank) / (float)params->number_of_processes * n_springs;
    lsm->s_end = (float)(params->process_rank + 1) / (float)params->number_of_processes * n_springs;
    lsm->s_proc = lsm->s_end - lsm->s_start;

    lsm->p_recvcounts = (int *)malloc(params->number_of_processes * sizeof(int));
    lsm->p_displs = (int *)malloc(params->number_of_processes * sizeof(int));

    lsm->s_recvcounts = (int *)malloc(params->number_of_processes * sizeof(int));
    lsm->s_displs = (int *)malloc(params->number_of_processes * sizeof(int));

    MPI_Allgather(&lsm->p_proc, 1, MPI_INT, lsm->p_recvcounts, 1, MPI_INT, params->comm_xslices);
    MPI_Allgather(&lsm->s_proc, 1, MPI_INT, lsm->s_recvcounts, 1, MPI_INT, params->comm_xslices);

    lsm->p_displs[0] = 0;
    lsm->s_displs[0] = 0;
    for (int i = 1; i < params->number_of_processes; i++)
    {
        lsm->p_displs[i] = lsm->p_displs[i - 1] + lsm->p_recvcounts[i - 1];
        lsm->s_displs[i] = lsm->s_displs[i - 1] + lsm->s_recvcounts[i - 1];
    }
}