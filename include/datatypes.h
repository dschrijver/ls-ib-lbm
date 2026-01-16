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

    // Relaxation time
    double tau;

    // Starting density
    double rho_0;

    // Gravity
    double gx;
    double gy;
    double gz;

    // LSM
    int N_connections_bulk;

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

    // NEBB parameters
    double C_norm;
    double C_par;
} Stencil;

typedef struct LSMParticle
{
    int id;

    double x; 
    double y; 
    double z; 

    double u; 
    double v; 
    double w; 

    double u_predict; 
    double v_predict;
    double w_predict; 

    double Fx_ext;
    double Fy_ext; 
    double Fz_ext; 

    double Fx_spring; 
    double Fy_spring; 
    double Fz_spring; 

    double Fx_damp; 
    double Fy_damp;
    double Fz_damp;

    double Fx_FSI;
    double Fy_FSI;
    double Fz_FSI; 

    double Fx; 
    double Fy;
    double Fz;

    int N_connections; 
    double chi;     

    struct LSMParticle *next;     
} LSMParticle;

typedef struct LSMSpring
{
    int id; 

    struct LSMParticle *particle_1; 
    struct LSMParticle *particle_2; 

    int particle_1_id;
    int particle_2_id;

    double l_eq;                    
    double S;                 

    int active; 

    struct LSMSpring *next; 
} LSMSpring;

typedef struct LSMBag
{
    int n_particles;
    int n_springs;        

    struct LSMParticle *particle_first;
    struct LSMSpring *spring_first;         
} LSMBag;

typedef struct SimulationBag
{
    struct ParamBag *params;
    struct DistributionBag *dists;
    struct FieldBag *fields;
    struct Stencil *stencil;
    struct LSMBag *lsm;
} SimulationBag;

#endif