#include <mpi.h>
#include <stdio.h>

#include "include/datatypes.h"
#include "include/memory.h"
#include "params.h"
#include "include/initialize.h"
#include "include/stencil.h"
#include "include/forcing.h"
#include "include/utils.h"
#include "include/output.h"
#include "include/collide.h"
#include "include/communicate.h"
#include "include/stream.h"
#include "include/fields.h"
#include "definitions.h"
#include "include/wetnode.h"
#include "include/lsm.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    SimulationBag *sim;
    DistributionBag *dists;
    FieldBag *fields;
    ParamBag *params;
    Stencil *stencil;
    LSMBag *lsm;

    allocate_bags(&sim, &dists, &fields, &params, &stencil, &lsm);

    set_params(params);

    initialize_MPI(params);

    initialize_stencil(sim);

    allocate_distributions(sim);

    allocate_fields(sim);

    params->t = 0;
    params->t_output = 0;
    params->t_log = 0;
    params->n_output = 0;

    initialize_fields(sim);

#ifdef LSIBM
    initialize_lsm(sim);
#endif

    evaluate_gravity_force(sim);

    evaluate_total_force(sim);

    initialize_distribution(sim);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    double start_timestep, duration_timestep;
    double start_substep, duration_substep;
    char output_info[128];

    while (params->t < params->NTIME)
    {

        MPI_Barrier(MPI_COMM_WORLD);
        start_timestep = MPI_Wtime();

        // LOGGING TIME
        if ((params->process_rank == 0) && (params->t_log == params->t)) 
        {
            printf("================================================================================\n");
            printf("Time: %d\n", params->t);
            printf("--------------------------------------------------------------------------------\n");
        }

        // OUTPUT
        if (params->t_output == params->t)
        {
            TIME_OUTPUT(
                output_data(sim);
                params->t_output += params->NSTORE;
            )
        }

        TIME("> Collision...",
            collide_distributions(sim);
        )

        TIME("> Communicate distributions...",
            communicate_dists(sim);
        )

        TIME("> Streaming...",
            stream_distributions(sim);
        )

#ifdef WETNODE
        TIME("> Wetnode boundary conditions...",
            wetnode_boundary_conditions(sim);
        )
#endif

        TIME("> Extracting moments...",
            extract_moments(sim);
        )

        TIME("> Evaluating fluid forces...",
            evaluate_gravity_force(sim);
#ifdef LSIBM
            update_preliminary_velocity(sim);
#endif
        )

#ifdef LSIBM
        TIME("> Updating particle positions...", 
            update_particle_positions(sim);
            update_particle_predicted_velocities(sim);
        )

        TIME("> Evaluating particle forces...",
            evaluate_particle_external_forces(sim);
            evaluate_spring_forces(sim);
            evaluate_damping_forces(sim);
            update_particle_preliminary_velocities(sim); 
        )

        TIME("> Evaluating FSI forces...",
            evaluate_particle_weights(sim);
            evaluate_FSI_forces(sim);
        )
#endif

        TIME("> Computing final velocity", 
            evaluate_total_force(sim);
#ifndef LSIBM
            update_final_velocity(sim);
#endif
        )

        duration_timestep = MPI_Wtime() - start_timestep;

        // LOGGING INFORMATION
        if ((params->process_rank == 0) && (params->t_log == params->t))   
        {
            printf("--------------------------------------------------------------------------------\n");
            printf("Step completed!\n");
            printf("    Duration of time step: %.4fs\n", duration_timestep);
            printf("    Total simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0);
            printf("    Expected remaining simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0 / (double)(params->t + 1) * (double)(params->NTIME - params->t - 1));
            params->t_log += params->NLOG;
        }

        params->t++;
    }

    output_data(sim);

    if (params->process_rank == 0)
    {
        printf("\nSimulation done!\n");
    }

    MPI_Finalize();
    return 0;
}