#include <stdlib.h>
#include <math.h>

#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/lsm.h"

void create_particle(double x, double y, double z, double u, double v, double w, LSMBag *lsm)
{
    lsm->n_particles++;
    lsm->particles = (LSMParticle *)realloc(lsm->particles, lsm->n_particles * sizeof(LSMParticle));

    LSMParticle *particles = lsm->particles;
    int n_particles = lsm->n_particles;

    particles[n_particles - 1].x = x;
    particles[n_particles - 1].y = y;
    particles[n_particles - 1].z = z;
    particles[n_particles - 1].u = u;
    particles[n_particles - 1].v = v;
    particles[n_particles - 1].w = w;

    particles[n_particles - 1].Fx_ext = 0.0;
    particles[n_particles - 1].Fy_ext = 0.0;
    particles[n_particles - 1].Fz_ext = 0.0;
    particles[n_particles - 1].Fx_spring = 0.0;
    particles[n_particles - 1].Fy_spring = 0.0;
    particles[n_particles - 1].Fz_spring = 0.0;
    particles[n_particles - 1].Fx_damp = 0.0;
    particles[n_particles - 1].Fy_damp = 0.0;
    particles[n_particles - 1].Fz_damp = 0.0;
    particles[n_particles - 1].Fx_IBM = 0.0;
    particles[n_particles - 1].Fy_IBM = 0.0;
    particles[n_particles - 1].Fz_IBM = 0.0;
    particles[n_particles - 1].Fx = 0.0;
    particles[n_particles - 1].Fy = 0.0;
    particles[n_particles - 1].Fz = 0.0;

    particles[n_particles - 1].N_connections = 0;
    particles[n_particles - 1].chi = 1.0;
}

void create_spring(int p1, int p2, double l_eq, LSMBag *lsm)
{
    lsm->n_springs++;
    lsm->springs = (LSMSpring *)realloc(lsm->springs, lsm->n_springs * sizeof(LSMSpring));

    LSMParticle *particles = lsm->particles;
    LSMSpring *springs = lsm->springs;
    int n_springs = lsm->n_springs;

    springs[n_springs - 1].p1 = p1;
    springs[n_springs - 1].p2 = p2;
    springs[n_springs - 1].l_eq = l_eq;
    springs[n_springs - 1].S = 0.0;
    springs[n_springs - 1].active = 1;

    particles[p1].N_connections++;
    particles[p2].N_connections++;
}

void add_springs(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    double x, y, z, r;

    double r_max = params->r_max;

    int n_particles = lsm->n_particles;

    LSMParticle *particles = lsm->particles;

    for (int p1 = 0; p1 < n_particles; p1++)
    {
        for (int p2 = p1 + 1; p2 < n_particles; p2++)
        {
            x = particles[p1].x - particles[p2].x;
            y = particles[p1].y - particles[p2].y;
            z = particles[p1].z - particles[p2].z;

            r = sqrt(x * x + y * y + z * z);

            if (r < r_max)
                create_spring(p1, p2, r, lsm);
        }
    }
}

void compute_chi(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    LSMParticle *particles = lsm->particles;

    int n_particles = lsm->n_particles;

    for (int p = 0; p < n_particles; p++)
    {
        if (particles[p].N_connections >= params->N_connections_bulk)
            particles[p].chi = 0.0;
        else
            particles[p].chi = 1.0;
    }
}

void update_particle_positions(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double m = params->m;

    double *local_particle_positions = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_particle_positions = (double *)malloc(n_particles * 3 * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        local_particle_positions[(p - p_start) * 3] = particles[p].u + 0.5 / m * particles[p].Fx;
        local_particle_positions[(p - p_start) * 3 + 1] = particles[p].v + 0.5 / m * particles[p].Fy;
        local_particle_positions[(p - p_start) * 3 + 2] = particles[p].w + 0.5 / m * particles[p].Fz;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_positions, p_proc, data_type, total_particle_positions, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].x += total_particle_positions[p * 3];
        particles[p].y += total_particle_positions[p * 3 + 1];
        particles[p].z += total_particle_positions[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_positions);
    free(local_particle_positions);
}

void update_particle_predicted_velocities(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double m = params->m;

    double *local_particle_predicted_velocities = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_particle_predicted_velocities = (double *)malloc(n_particles * 3 * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        local_particle_predicted_velocities[(p - p_start) * 3] = particles[p].u + 1.0 / m * particles[p].Fx;
        local_particle_predicted_velocities[(p - p_start) * 3 + 1] = particles[p].v + 1.0 / m * particles[p].Fy;
        local_particle_predicted_velocities[(p - p_start) * 3 + 2] = particles[p].w + 1.0 / m * particles[p].Fz;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_predicted_velocities, p_proc, data_type, total_particle_predicted_velocities, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].u_predict = total_particle_predicted_velocities[p * 3];
        particles[p].v_predict = total_particle_predicted_velocities[p * 3 + 1];
        particles[p].w_predict = total_particle_predicted_velocities[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_predicted_velocities);
    free(local_particle_predicted_velocities);
}

// NOT PARALLELIZED
void evaluate_particle_external_forces(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    double m = params->m;
    double gx = params->gx;
    double gy = params->gy;
    double gz = params->gz;

    for (int p = 0; p < n_particles; p++)
    {
        particles[p].Fx_ext = m * gx;
        particles[p].Fy_ext = m * gy;
        particles[p].Fz_ext = m * gz;
    }
}

void evaluate_spring_forces(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    double rx, ry, rz, l, l_eq, Fx, Fy, Fz, S;
    int p1, p2;

    int n_springs = lsm->n_springs;
    LSMSpring *springs = lsm->springs;
    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int s_start = lsm->s_start;
    int s_end = lsm->s_end;
    int s_proc = lsm->s_proc;

    double k = params->k;

    double *local_spring_forces = (double *)malloc(s_proc * 4 * sizeof(double));
    double *total_spring_forces = (double *)malloc(n_springs * 4 * sizeof(double));

    for (int s = s_start; s < s_end; s++)
    {
        if (!springs[s].active)
        {
            local_spring_forces[(s - s_start) * 4] = 0.0;
            local_spring_forces[(s - s_start) * 4 + 1] = 0.0;
            local_spring_forces[(s - s_start) * 4 + 2] = 0.0;
            local_spring_forces[(s - s_start) * 4 + 3] = springs[s].S;
            continue;
        }

        p1 = springs[s].p1;
        p2 = springs[s].p2;

        rx = particles[p1].x - particles[p2].x;
        ry = particles[p1].y - particles[p2].y;
        rz = particles[p1].z - particles[p2].z;

        l = sqrt(rx * rx + ry * ry + rz * rz);
        l_eq = springs[s].l_eq;

        Fx = -k / (l_eq * l_eq) * (l - l_eq) * rx / l;
        Fy = -k / (l_eq * l_eq) * (l - l_eq) * ry / l;
        Fz = -k / (l_eq * l_eq) * (l - l_eq) * rz / l;
        S = (l - l_eq) / l_eq;

        local_spring_forces[(s - s_start) * 4] = Fx;
        local_spring_forces[(s - s_start) * 4 + 1] = Fy;
        local_spring_forces[(s - s_start) * 4 + 2] = Fz;
        local_spring_forces[(s - s_start) * 4 + 3] = S;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(4, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_spring_forces, s_proc, data_type, total_spring_forces, lsm->s_recvcounts, lsm->s_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].Fx_spring = 0.0;
        particles[p].Fy_spring = 0.0;
        particles[p].Fz_spring = 0.0;
    }

    for (int s = 0; s < n_springs; s++)
    {
        Fx = total_spring_forces[s * 4];
        Fy = total_spring_forces[s * 4 + 1];
        Fz = total_spring_forces[s * 4 + 2];

        p1 = springs[s].p1;
        p2 = springs[s].p2;

        particles[p1].Fx_spring += Fx;
        particles[p1].Fy_spring += Fy;
        particles[p1].Fz_spring += Fz;

        particles[p2].Fx_spring -= Fx;
        particles[p2].Fy_spring -= Fy;
        particles[p2].Fz_spring -= Fz;

        springs[s].S = total_spring_forces[s * 4 + 3];
    }

    MPI_Type_free(&data_type);
    free(total_spring_forces);
    free(local_spring_forces);
}

void evaluate_damping_forces(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    double rx, ry, rz, l, vr, Fx, Fy, Fz;
    int p1, p2;

    int n_springs = lsm->n_springs;
    LSMSpring *springs = lsm->springs;
    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int s_start = lsm->s_start;
    int s_end = lsm->s_end;
    int s_proc = lsm->s_proc;

    double c = params->c;

    double *local_damp_forces = (double *)malloc(s_proc * 3 * sizeof(double));
    double *total_damp_forces = (double *)malloc(n_springs * 3 * sizeof(double));

    for (int s = s_start; s < s_end; s++)
    {
        if (!springs[s].active)
        {
            local_damp_forces[(s - s_start) * 3] = 0.0;
            local_damp_forces[(s - s_start) * 3 + 1] = 0.0;
            local_damp_forces[(s - s_start) * 3 + 2] = 0.0;
            continue;
        }

        p1 = springs[s].p1;
        p2 = springs[s].p2;

        rx = particles[p1].x - particles[p2].x;
        ry = particles[p1].y - particles[p2].y;
        rz = particles[p1].z - particles[p2].z;

        l = sqrt(rx * rx + ry * ry + rz * rz);

        vr = (particles[p1].u_predict - particles[p2].u_predict) * rx / l + (particles[p1].v_predict - particles[p2].v_predict) * ry / l + (particles[p1].w_predict - particles[p2].w_predict) * rz / l;

        Fx = -c * vr * rx / l;
        Fy = -c * vr * ry / l;
        Fz = -c * vr * rz / l;

        local_damp_forces[(s - s_start) * 3] = Fx;
        local_damp_forces[(s - s_start) * 3 + 1] = Fy;
        local_damp_forces[(s - s_start) * 3 + 2] = Fz;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_damp_forces, s_proc, data_type, total_damp_forces, lsm->s_recvcounts, lsm->s_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].Fx_damp = 0.0;
        particles[p].Fy_damp = 0.0;
        particles[p].Fz_damp = 0.0;
    }

    for (int s = 0; s < n_springs; s++)
    {
        Fx = total_damp_forces[s * 3];
        Fy = total_damp_forces[s * 3 + 1];
        Fz = total_damp_forces[s * 3 + 2];

        p1 = springs[s].p1;
        p2 = springs[s].p2;

        particles[p1].Fx_damp += Fx;
        particles[p1].Fy_damp += Fy;
        particles[p1].Fz_damp += Fz;

        particles[p2].Fx_damp -= Fx;
        particles[p2].Fy_damp -= Fy;
        particles[p2].Fz_damp -= Fz;
    }

    MPI_Type_free(&data_type);
    free(total_damp_forces);
    free(local_damp_forces);
}

void update_particle_preliminary_velocities(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double m = params->m;

    double *local_particle_preliminary_velocities = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_particle_preliminary_velocities = (double *)malloc(n_particles * 3 * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        local_particle_preliminary_velocities[(p - p_start) * 3] = 0.5 / m * (particles[p].Fx + particles[p].Fx_spring + particles[p].Fx_damp + particles[p].Fx_ext);
        local_particle_preliminary_velocities[(p - p_start) * 3 + 1] = 0.5 / m * (particles[p].Fy + particles[p].Fy_spring + particles[p].Fy_damp + particles[p].Fy_ext);
        local_particle_preliminary_velocities[(p - p_start) * 3 + 2] = 0.5 / m * (particles[p].Fz + particles[p].Fz_spring + particles[p].Fz_damp + particles[p].Fz_ext);
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_preliminary_velocities, p_proc, data_type, total_particle_preliminary_velocities, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].u += total_particle_preliminary_velocities[p * 3];
        particles[p].v += total_particle_preliminary_velocities[p * 3 + 1];
        particles[p].w += total_particle_preliminary_velocities[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_preliminary_velocities);
    free(local_particle_preliminary_velocities);
}

void evaluate_particle_weights(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double *local_particle_weight = (double *)malloc(p_proc * sizeof(double));
    double *total_particle_weight = (double *)malloc(n_particles * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        particles[p].weight = 0.0;
        
    }

    free(local_particle_weight);
    free(total_particle_weight);
}

void evaluate_FSI_forces(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double rho_0 = params->rho_0;
    double m = params->m;
}