#ifndef DATATYPES_H
#define DATATYPES_H

#include <mpi.h>

typedef struct ParamBag
{
    // Simulation variables
    int t;       
    int t_output; 
    int t_log;   
    int n_output; 

    // General parameters
    int NTIME;  
    int NSTORE; 
    int NLOG;  
    int NX;   
    int NY;   
    int NZ;    

    // Relaxation times
    double tau;
    double rho_0;

    double gx;
    double gy;
    double gz;

    // MPI
    int number_of_processes;  ///< Stores number of processes
    int process_rank;         ///< Process rank
    int process_coords[3];    ///< Coords of processor in virtual MPI topology
    int process_neighbors[2]; ///< Left and right processor neighbors in virtual MPI topology.
    int i_start;              ///< Starting index of slab owned by current processor.
    int i_end;                ///< Ending index of slab owned by current processor.
    int NX_proc;              ///< Number of nodes owned by current processor, NX_proc = i_end - i_start.
    MPI_Comm comm_xslices;    ///< Communicator of slab decomposition.
} ParamBag;

typedef struct DistributionBag
{
    double *f1;
    double *f2;
} DistributionBag;

typedef struct FieldBag
{
    // Densities and pressure
    double *rho;
    double *pressure;

    // Velocities
    double *u;
    double *v;
    double *w;

    // Forces
    double *Fx;
    double *Fy;
    double *Fz;

    double *Fx_grav;
    double *Fy_grav;
    double *Fz_grav;
} FieldBag;

typedef struct Stencil
{
    double cs2; // 1/3

    // Stencil
    int NP;
    int *cx;
    int *cy;
    int *cz;
    double *wp;
    int *p_bounceback;

    double C_norm;
    double C_par;
} Stencil;

typedef struct SimulationBag
{
    struct ParamBag *params;
    struct DistributionBag *dists;
    struct FieldBag *fields;
    struct Stencil *stencil;
} SimulationBag;

#endif