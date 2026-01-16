#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/lsm.h"

void create_particle(double x, double y, double z, double u, double v, double w, LSMBag *lsm)
{
    LSMParticle *particle_old = lsm->particle_first;
    LSMParticle *particle_new = (LSMParticle *)malloc(sizeof(LSMParticle));
    lsm->particle_first = particle_new;

    particle_new->next = particle_old;

    particle_new->x = x;
    particle_new->y = y;
    particle_new->z = z;
    particle_new->u = u;
    particle_new->v = v;
    particle_new->w = w;

    particle_new->Fx_ext = 0.0;
    particle_new->Fy_ext = 0.0;
    particle_new->Fz_ext = 0.0;
    particle_new->Fx_spring = 0.0;
    particle_new->Fy_spring = 0.0;
    particle_new->Fz_spring = 0.0;
    particle_new->Fx_damp = 0.0;
    particle_new->Fy_damp = 0.0;
    particle_new->Fz_damp = 0.0;
    particle_new->Fx_FSI = 0.0;
    particle_new->Fy_FSI = 0.0;
    particle_new->Fz_FSI = 0.0;
    particle_new->Fx = 0.0;
    particle_new->Fy = 0.0;
    particle_new->Fz = 0.0;

    particle_new->N_connections = 0;
    particle_new->chi = 1.0;

    lsm->n_particles++;
}

void add_springs(SimulationBag *sim)
{
    double x, y, z, r, r2, l_eq;

    LSMParticle *particle_i, *particle_j;
    LSMSpring *spring;
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int N = lsm->n_particles;

    if (N == 0)
        return;

    MPI_Comm comm_xslices = params->comm_xslices;

    double r_max = pow((double)(params->N_connections_bulk + 1) / (4.0 / 3.0 * DS_PI), 1.0 / 3.0);
    double r_max_2 = r_max * r_max;

    double start_part, duration_part;

    MPI_Barrier(MPI_COMM_WORLD);
    start_part = MPI_Wtime();
    if (params->process_rank == 0)
    {
        printf("%-78s", "\033[2;37m\n  > Add springs...");
    }

    char print_string[128];

    if (N < params->number_of_processes)
    {
        particle_i = lsm->particle_first;
        while (particle_i != NULL)
        {

            particle_j = particle_i->next;
            while (particle_j != NULL)
            {
                x = particle_i->x - particle_j->x;
                y = particle_i->y - particle_j->y;
                z = particle_i->z - particle_j->z;
                r = sqrt(x * x + y * y + z * z);

                if (r < r_max)
                {

                    spring = spring_exists(lsm, particle_i, particle_j);

                    if (spring == NULL)
                    {
                        create_spring(lsm, particle_i, particle_j, r);
                    }
                    else if (!spring->active)
                    {
                        l_eq = spring->l_eq;
                        if ((fabs(r - l_eq) / l_eq < params->S_reactivate))
                        {
                            spring->active = 1;
                            particle_i->N_connections++;
                            particle_j->N_connections++;
                        }
                    }
                }

                particle_j = particle_j->next;
            }

            particle_i = particle_i->next;
        }
        return;
    }

    // Create local particle list and identify global particles
    int particle_start = (float)(params->process_coords[0]) / (float)params->number_of_processes * N;
    int particle_end = (float)(params->process_coords[0] + 1) / (float)params->number_of_processes * N;

    lsm->particle_first_local = NULL;
    lsm->n_particles_local = 0;
    LSMParticle *particle = lsm->particle_first;
    LSMParticle *particle_local;
    int id = 0;
    while (particle != NULL)
    {
        particle->id = id;
        if ((particle->id >= particle_start) && (particle->id < particle_end))
        {
            particle_local = lsm->particle_first_local;
            lsm->particle_first_local = particle;
            particle->next_local = particle_local;
            lsm->n_particles_local++;
        }
        particle->id_local = -1;
        id++;
        particle = particle->next;
    }

    // Identify local particles
    particle_local = lsm->particle_first_local;
    int id_local = 0;
    while (particle_local != NULL)
    {
        particle_local->id_local = id_local;
        id_local++;

        particle_local = particle_local->next_local;
    }

    LSMSpring **A = (LSMSpring **)malloc(lsm->n_particles_local * N * sizeof(LSMSpring *));
    memset(A, 0, lsm->n_particles_local * N * sizeof(LSMSpring *));

    spring = lsm->spring_first;
    int spring_id = 0;
    int id_1, id_2, id_local_1, id_local_2;
    while (spring != NULL)
    {
        spring->id = spring_id;
        spring_id++;

        id_1 = spring->particle_1->id;
        id_2 = spring->particle_2->id;

        id_local_1 = spring->particle_1->id_local;
        id_local_2 = spring->particle_2->id_local;

        if (id_local_1 > -1)
        {
            A[N * id_local_1 + id_2] = spring;
        }

        if (id_local_2 > -1)
        {
            A[N * id_local_2 + id_1] = spring;
        }

        spring = spring->next;
    }

    PreliminarySpring *local_springs_to_add = (PreliminarySpring *)malloc(sizeof(PreliminarySpring));
    int *local_springs_to_activate = (int *)malloc(sizeof(int));
    int n_added_local = 0;
    int n_activated_local = 0;

    particle_i = lsm->particle_first_local;
    while (particle_i != NULL)
    {

        particle_j = particle_i->next;
        while (particle_j != NULL)
        {
            x = particle_i->x - particle_j->x;
            y = particle_i->y - particle_j->y;
            z = particle_i->z - particle_j->z;
            r2 = x * x + y * y + z * z;

            if (r2 < r_max_2)
            {
                r = sqrt(r2);

                spring = A[N * particle_i->id_local + particle_j->id];

                if (spring == NULL)
                {
                    n_added_local++;
                    local_springs_to_add = (PreliminarySpring *)realloc(local_springs_to_add, n_added_local * sizeof(PreliminarySpring));
                    local_springs_to_add[n_added_local - 1].l_eq = r;
                    local_springs_to_add[n_added_local - 1].particle_1_id = particle_i->id;
                    local_springs_to_add[n_added_local - 1].particle_2_id = particle_j->id;
                }
                else if (!spring->active)
                {
                    l_eq = spring->l_eq;
                    if ((fabs(r - l_eq) / l_eq < params->S_reactivate))
                    {
                        n_activated_local++;
                        local_springs_to_activate = (int *)realloc(local_springs_to_activate, n_activated_local * sizeof(int));
                        local_springs_to_activate[n_activated_local - 1] = spring->id;
                    }
                }
            }

            particle_j = particle_j->next;
        }

        particle_i = particle_i->next_local;
    }

    int *n_added_array = (int *)malloc(params->number_of_processes * sizeof(int));
    MPI_Allgather(&n_added_local, 1, MPI_INT, n_added_array, 1, MPI_INT, comm_xslices);

    int *n_activated_array = (int *)malloc(params->number_of_processes * sizeof(int));
    MPI_Allgather(&n_activated_local, 1, MPI_INT, n_activated_array, 1, MPI_INT, comm_xslices);

    int n_added_total = 0.0;
    for (int i = 0; i < params->number_of_processes; i++)
    {
        n_added_total += n_added_array[i];
        n_added_array[i] *= sizeof(PreliminarySpring);
    }

    int n_activated_total = 0.0;
    for (int i = 0; i < params->number_of_processes; i++)
    {
        n_activated_total += n_activated_array[i];
    }

    PreliminarySpring *total_springs_to_add = (PreliminarySpring *)malloc(n_added_total * sizeof(PreliminarySpring));

    int *total_springs_to_activate = (int *)malloc(n_activated_total * sizeof(int));

    int *displs = (int *)malloc(params->number_of_processes * sizeof(int));

    displs[0] = 0;
    for (int i = 1; i < params->number_of_processes; i++)
    {
        displs[i] = displs[i - 1] + n_added_array[i - 1];
    }

    MPI_Allgatherv(local_springs_to_add, n_added_local * sizeof(PreliminarySpring), MPI_BYTE, total_springs_to_add, n_added_array, displs, MPI_BYTE, comm_xslices);

    displs[0] = 0;
    for (int i = 1; i < params->number_of_processes; i++)
    {
        displs[i] = displs[i - 1] + n_activated_array[i - 1];
    }

    MPI_Allgatherv(local_springs_to_activate, n_activated_local, MPI_INT, total_springs_to_activate, n_activated_array, displs, MPI_INT, comm_xslices);

    LSMParticle *pair[2] = {0};
    int num;

    for (int i = 0; i < n_activated_total; i++)
    {
        spring = lsm->spring_first;
        while (spring != NULL)
        {
            if (spring->id == total_springs_to_activate[i])
            {
                spring->active = 1;
                spring->particle_1->N_connections++;
                spring->particle_2->N_connections++;
            }

            spring = spring->next;
        }
    }

    for (int i = 0; i < n_added_total; i++)
    {
        particle = lsm->particle_first;

        id_1 = total_springs_to_add[i].particle_1_id;
        id_2 = total_springs_to_add[i].particle_2_id;

        num = 0;
        while (particle != NULL)
        {
            if ((particle->id == id_1) || (particle->id == id_2))
            {
                pair[num] = particle;
                num++;
                if (num == 2)
                    break;
            }

            particle = particle->next;
        }

        create_spring(lsm, pair[0], pair[1], total_springs_to_add[i].l_eq);
    }

    if (params->process_rank == 0)
    {
        sprintf(print_string, "\n        Total springs: %d, springs added: %d, springs activated: %d", lsm->n_springs, n_added_total, n_activated_total);
        printf("%-71s", print_string);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    duration_part = MPI_Wtime() - start_part;
    if (params->process_rank == 0)
    {
        printf("[%7.4fs]\033[0m", duration_part);
    }

    free(A);

    free(local_springs_to_add);
    free(local_springs_to_activate);

    free(n_added_array);
    free(n_activated_array);
    free(displs);

    free(total_springs_to_add);
    free(total_springs_to_activate);
}