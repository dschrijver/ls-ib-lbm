#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../definitions.h"
#include "../include/communicate.h"
#include "../include/output.h"
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
    particles[n_particles - 1].chi = 1;
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
            particles[p].chi = 0;
        else
            particles[p].chi = 1;
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
        local_particle_positions[(p - p_start) * 3] = particles[p].x + particles[p].u + 0.5 / m * particles[p].Fx;
        local_particle_positions[(p - p_start) * 3 + 1] = particles[p].y + particles[p].v + 0.5 / m * particles[p].Fy;
        local_particle_positions[(p - p_start) * 3 + 2] = particles[p].z + particles[p].w + 0.5 / m * particles[p].Fz;

#ifdef XPERIODIC
        local_particle_positions[(p - p_start) * 3] = fmod(local_particle_positions[(p - p_start) * 3], (double)params->NX);
#endif
#ifdef YPERIODIC
        local_particle_positions[(p - p_start) * 3 + 1] = fmod(local_particle_positions[(p - p_start) * 3 + 1], (double)params->NY);
#endif
#ifdef ZPERIODIC
        local_particle_positions[(p - p_start) * 3 + 2] = fmod(local_particle_positions[(p - p_start) * 3 + 2], (double)params->NZ);
#endif
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_positions, p_proc, data_type, total_particle_positions, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].x = total_particle_positions[p * 3];
        particles[p].y = total_particle_positions[p * 3 + 1];
        particles[p].z = total_particle_positions[p * 3 + 2];
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
        local_particle_preliminary_velocities[(p - p_start) * 3] = particles[p].u + 0.5 / m * (particles[p].Fx + particles[p].Fx_spring + particles[p].Fx_damp + particles[p].Fx_ext);
        local_particle_preliminary_velocities[(p - p_start) * 3 + 1] = particles[p].v + 0.5 / m * (particles[p].Fy + particles[p].Fy_spring + particles[p].Fy_damp + particles[p].Fy_ext);
        local_particle_preliminary_velocities[(p - p_start) * 3 + 2] = particles[p].w + 0.5 / m * (particles[p].Fz + particles[p].Fz_spring + particles[p].Fz_damp + particles[p].Fz_ext);
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_preliminary_velocities, p_proc, data_type, total_particle_preliminary_velocities, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].u = total_particle_preliminary_velocities[p * 3];
        particles[p].v = total_particle_preliminary_velocities[p * 3 + 1];
        particles[p].w = total_particle_preliminary_velocities[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_preliminary_velocities);
    free(local_particle_preliminary_velocities);
}

void interpolate_fluid_velocity_to_particle(double *uf, double *vf, double *wf, SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    LSMParticle *particles = lsm->particles;

    MPI_Comm comm_xslices = params->comm_xslices;

    int i_start_loop, j_start_loop, k_start_loop;
    int i_end_loop, j_end_loop, k_end_loop;
    double x, y, z;
    double left, right;
    double rx, ry, rz;
    double uf_i, vf_i, wf_i;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int NY = params->NY;
    int NZ = params->NZ;

    int n_particles = lsm->n_particles;

    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    for (int p = 0; p < n_particles; p++)
    {
        left = (double)i_start;
        right = (double)i_end;
        if ((particles[p].x < left) || (particles[p].x >= right))
        {
            uf[p] = 0.0;
            vf[p] = 0.0;
            wf[p] = 0.0;
            continue;
        }

#ifdef LEFT_NEBB_VELOCITY
        i_start_loop = (int)particles[p].x;
#else
        i_start_loop = (int)(particles[p].x - 0.5);
#endif

#ifdef BOTTOM_NEBB_VELOCITY
        j_start_loop = (int)particles[p].y;
#else
        j_start_loop = (int)(particles[p].y - 0.5);
#endif

#ifdef BACK_NEBB_VELOCITY
        k_start_loop = (int)particles[p].z;
#else
        k_start_loop = (int)(particles[p].z - 0.5);
#endif

#ifdef XPERIODIC
        i_start_loop = max(i_start_loop, i_start - 1);
        i_end_loop = min(i_start_loop + 2, i_end + 1);
#else
        i_start_loop = max(i_start_loop, max(i_start - 1, 0));
        i_end_loop = min(i_start_loop + 2, min(i_end + 1, params->NX));
#endif

#ifdef YPERIODIC
        j_start_loop = max(j_start_loop, -1);
        j_end_loop = min(j_start_loop + 2, NY + 1);
#else
        j_start_loop = max(j_start_loop, 0);
        j_end_loop = min(j_start_loop + 2, NY);
#endif

#ifdef ZPERIODIC
        k_start_loop = max(k_start_loop, -1);
        k_end_loop = min(k_start_loop + 2, NZ + 1);
#else
        k_start_loop = max(k_start_loop, 0);
        k_end_loop = min(k_start_loop + 2, NZ);
#endif

        uf_i = 0.0;
        vf_i = 0.0;
        wf_i = 0.0;
        for (int i = i_start_loop; i < i_end_loop; i++)
        {
            for (int j = j_start_loop; j < j_end_loop; j++)
            {
                for (int k = k_start_loop; k < k_end_loop; k++)
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

                    rx = fabs(x - particles[p].x);
                    ry = fabs(y - particles[p].y);
                    rz = fabs(z - particles[p].z);

                    uf_i += u[INDEX(i, j, k)] * kernel(rx, ry, rz);
#ifdef YPERIODIC
                    vf_i += v[INDEX(i, mod(j, NY), k)] * kernel(rx, ry, rz);
#else
                    vf_i += v[INDEX(i, j, k)] * kernel(rx, ry, rz);
#endif
#ifdef ZPERIODIC
                    wf_i += w[INDEX(i, j, mod(k, NZ))] * kernel(rx, ry, rz);
#else
                    wf_i += w[INDEX(i, j, k)] * kernel(rx, ry, rz);
#endif
                }
            }
        }
        uf[p] = uf_i;
        vf[p] = vf_i;
        wf[p] = wf_i;
    }

    MPI_Allreduce(MPI_IN_PLACE, uf, n_particles, MPI_DOUBLE, MPI_SUM, comm_xslices);
    MPI_Allreduce(MPI_IN_PLACE, vf, n_particles, MPI_DOUBLE, MPI_SUM, comm_xslices);
    MPI_Allreduce(MPI_IN_PLACE, wf, n_particles, MPI_DOUBLE, MPI_SUM, comm_xslices);
}

void spread_FSI_forces_to_fluid(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    LSMParticle *particles = lsm->particles;

    double x, y, z;
    double rx, ry, rz;
    double Fx_IBM_i, Fy_IBM_i, Fz_IBM_i, Fx_rigid_i, Fy_rigid_i, Fz_rigid_i;

    int i_start_loop, j_start_loop, k_start_loop;
    int i_end_loop, j_end_loop, k_end_loop;
    int j_mod, k_mod;

    int n_particles = lsm->n_particles;

    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;
    int NX_proc = params->NX_proc;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *Fx_IBM_diff = fields->Fx_IBM_diff;
    double *Fy_IBM_diff = fields->Fy_IBM_diff;
    double *Fz_IBM_diff = fields->Fz_IBM_diff;

    double *Fx_rigid_diff = fields->Fx_rigid_diff;
    double *Fy_rigid_diff = fields->Fy_rigid_diff;
    double *Fz_rigid_diff = fields->Fz_rigid_diff;

    double *Fx_IBM = fields->Fx_IBM;
    double *Fy_IBM = fields->Fy_IBM;
    double *Fz_IBM = fields->Fz_IBM;

    double *Fx_rigid = fields->Fx_rigid;
    double *Fy_rigid = fields->Fy_rigid;
    double *Fz_rigid = fields->Fz_rigid;

    memset(&Fx_IBM_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fy_IBM_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fz_IBM_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fx_rigid_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fy_rigid_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fz_rigid_diff[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));

    for (int p = 0; p < n_particles; p++)
    {
#ifdef LEFT_NEBB_VELOCITY
        i_start_loop = (int)particles[p].x;
#else
        i_start_loop = (int)(particles[p].x - 0.5);
#endif

#ifdef BOTTOM_NEBB_VELOCITY
        j_start_loop = (int)particles[p].y;
#else
        j_start_loop = (int)(particles[p].y - 0.5);
#endif

#ifdef BACK_NEBB_VELOCITY
        k_start_loop = (int)particles[p].z;
#else
        k_start_loop = (int)(particles[p].z - 0.5);
#endif

#ifdef XPERIODIC
        if ((i_start == 0) && (i_start_loop == NX - 1))
        {
            i_start_loop = 0;
            i_end_loop = 1;
        }
        else if ((i_end == NX) && (i_start_loop == -1))
        {
            i_start_loop = NX - 1;
            i_end_loop = NX;
        }
        else
        {
            i_start_loop = max(i_start_loop, i_start);
            i_end_loop = min(i_start_loop + 2, i_end);
        }
#else
        i_start_loop = max(i_start_loop, i_start);
        i_end_loop = min(i_start_loop + 2, i_end);
#endif
        j_end_loop = j_start_loop + 2;
        k_end_loop = k_start_loop + 2;

        for (int i = i_start_loop; i < i_end_loop; i++)
        {
            for (int j = j_start_loop; j < j_end_loop; j++)
            {
                for (int k = k_start_loop; k < k_end_loop; k++)
                {
                    j_mod = mod(j, NY);
                    k_mod = mod(k, NZ);
#ifdef LEFT_NEBB_VELOCITY
                    if (i == 0)
                        continue;
#endif
#ifdef RIGHT_NEBB_VELOCITY
                    if (i == NX - 1)
                        continue;
#endif
#ifdef BOTTOM_NEBB_VELOCITY
                    if (j_mod == 0)
                        continue;
#endif
#ifdef TOP_NEBB_VELOCITY
                    if (j_mod == NY - 1)
                        continue;
#endif
#ifdef BACK_NEBB_VELOCITY
                    if (k_mod == 0)
                        continue;
#endif
#ifdef FRONT_NEBB_VELOCITY
                    if (k_mod == NZ - 1)
                        continue;
#endif

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
                    rx = fabs(x - particles[p].x);
                    ry = fabs(y - particles[p].y);
                    rz = fabs(z - particles[p].z);
#ifdef XPERIODIC
                    rx = fmin(rx, (double)NX - rx);
#endif

                    Fx_IBM_i = (double)particles[p].chi * particles[p].Fx_FSI_diff * kernel(rx, ry, rz);
                    Fy_IBM_i = (double)particles[p].chi * particles[p].Fy_FSI_diff * kernel(rx, ry, rz);
                    Fz_IBM_i = (double)particles[p].chi * particles[p].Fz_FSI_diff * kernel(rx, ry, rz);

                    Fx_rigid_i = (1.0 - (double)particles[p].chi) * particles[p].Fx_FSI_diff * kernel(rx, ry, rz);
                    Fy_rigid_i = (1.0 - (double)particles[p].chi) * particles[p].Fy_FSI_diff * kernel(rx, ry, rz);
                    Fz_rigid_i = (1.0 - (double)particles[p].chi) * particles[p].Fz_FSI_diff * kernel(rx, ry, rz);

                    Fx_IBM_diff[INDEX(i, j_mod, k_mod)] += Fx_IBM_i;
                    Fy_IBM_diff[INDEX(i, j_mod, k_mod)] += Fy_IBM_i;
                    Fz_IBM_diff[INDEX(i, j_mod, k_mod)] += Fz_IBM_i;

                    Fx_rigid_diff[INDEX(i, j_mod, k_mod)] += Fx_rigid_i;
                    Fy_rigid_diff[INDEX(i, j_mod, k_mod)] += Fy_rigid_i;
                    Fz_rigid_diff[INDEX(i, j_mod, k_mod)] += Fz_rigid_i;

                    Fx_IBM[INDEX(i, j_mod, k_mod)] += Fx_IBM_i;
                    Fy_IBM[INDEX(i, j_mod, k_mod)] += Fy_IBM_i;
                    Fz_IBM[INDEX(i, j_mod, k_mod)] += Fz_IBM_i;

                    Fx_rigid[INDEX(i, j_mod, k_mod)] += Fx_rigid_i;
                    Fy_rigid[INDEX(i, j_mod, k_mod)] += Fy_rigid_i;
                    Fz_rigid[INDEX(i, j_mod, k_mod)] += Fz_rigid_i;
                }
            }
        }
    }
}

void shift_fluid_velocity_FSI(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    double rho_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = fields->rho;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *Fx_IBM_diff = fields->Fx_IBM_diff;
    double *Fy_IBM_diff = fields->Fy_IBM_diff;
    double *Fz_IBM_diff = fields->Fz_IBM_diff;

    double *Fx_rigid_diff = fields->Fx_rigid_diff;
    double *Fy_rigid_diff = fields->Fy_rigid_diff;
    double *Fz_rigid_diff = fields->Fz_rigid_diff;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX(i, j, k)];

        u[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * (Fx_IBM_diff[INDEX(i, j, k)] + Fx_rigid_diff[INDEX(i, j, k)]);
        v[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * (Fy_IBM_diff[INDEX(i, j, k)] + Fy_rigid_diff[INDEX(i, j, k)]);
        w[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * (Fz_IBM_diff[INDEX(i, j, k)] + Fz_rigid_diff[INDEX(i, j, k)]);
    }
}

void evaluate_FSI_forces(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    char print_string[128];

    double local_error, total_error;

    MPI_Comm comm_xslices = params->comm_xslices;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double rho_0 = params->rho_0;
    double m = params->m;

    double tol = params->tol;
    int max_iters = params->max_iters;

    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;
    int i_start = params->i_start;

    double *Fx_IBM = fields->Fx_IBM;
    double *Fy_IBM = fields->Fy_IBM;
    double *Fz_IBM = fields->Fz_IBM;

    double *Fx_rigid = fields->Fx_rigid;
    double *Fy_rigid = fields->Fy_rigid;
    double *Fz_rigid = fields->Fz_rigid;

    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    double *uf = (double *)malloc(n_particles * sizeof(double));
    double *vf = (double *)malloc(n_particles * sizeof(double));
    double *wf = (double *)malloc(n_particles * sizeof(double));

    double *local_F_FSI_diff = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_F_FSI_diff = (double *)malloc(n_particles * 3 * sizeof(double));

    communicate_fields(sim);
    interpolate_fluid_velocity_to_particle(uf, vf, wf, sim);

    for (int p = p_start; p < p_end; p++)
    {
        if (particles[p].chi)
        {
            local_F_FSI_diff[3 * (p - p_start)] = 2.0 * m * rho_0 / (m + rho_0) * (particles[p].u - uf[p]);
            local_F_FSI_diff[3 * (p - p_start) + 1] = 2.0 * m * rho_0 / (m + rho_0) * (particles[p].v - vf[p]);
            local_F_FSI_diff[3 * (p - p_start) + 2] = 2.0 * m * rho_0 / (m + rho_0) * (particles[p].w - wf[p]);
        }
        else
        {
            local_F_FSI_diff[3 * (p - p_start)] = 2.0 * rho_0 * (particles[p].u - uf[p]);
            local_F_FSI_diff[3 * (p - p_start) + 1] = 2.0 * rho_0 * (particles[p].v - vf[p]);
            local_F_FSI_diff[3 * (p - p_start) + 2] = 2.0 * rho_0 * (particles[p].w - wf[p]);
        }
    }

    // Communicate
    MPI_Allgatherv(local_F_FSI_diff, p_proc, data_type, total_F_FSI_diff, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].Fx_FSI_diff = total_F_FSI_diff[p * 3];
        particles[p].Fy_FSI_diff = total_F_FSI_diff[p * 3 + 1];
        particles[p].Fz_FSI_diff = total_F_FSI_diff[p * 3 + 2];

        if (particles[p].chi)
        {
            particles[p].Fx_IBM = -total_F_FSI_diff[p * 3];
            particles[p].Fy_IBM = -total_F_FSI_diff[p * 3 + 1];
            particles[p].Fz_IBM = -total_F_FSI_diff[p * 3 + 2];
        }
        else
        {
            particles[p].Fx_IBM = 0.0;
            particles[p].Fy_IBM = 0.0;
            particles[p].Fz_IBM = 0.0;
        }
    }

    memset(&Fx_IBM[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fy_IBM[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fz_IBM[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fx_rigid[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fy_rigid[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));
    memset(&Fz_rigid[INDEX(i_start, 0, 0)], 0, NX_proc * NY * NZ * sizeof(double));

    spread_FSI_forces_to_fluid(sim);
    shift_fluid_velocity_FSI(sim);

    int iters = 1;
    while (iters < max_iters)
    {
        local_error = 0.0;
        for (int p = p_start; p < p_end; p++)
        {
            local_error += particles[p].Fx_FSI_diff * particles[p].Fx_FSI_diff + particles[p].Fy_FSI_diff * particles[p].Fy_FSI_diff + particles[p].Fz_FSI_diff * particles[p].Fz_FSI_diff;
        }
        MPI_Allreduce(&local_error, &total_error, 1, MPI_DOUBLE, MPI_SUM, comm_xslices);
        if (sqrt(total_error) / (double)n_particles < tol)
        {
            break;
        }

        communicate_fields(sim);
        interpolate_fluid_velocity_to_particle(uf, vf, wf, sim);

        for (int p = p_start; p < p_end; p++)
        {
            if (particles[p].chi)
            {
                local_F_FSI_diff[3 * (p - p_start)] = 2.0 * m * (particles[p].u - uf[p]) + particles[p].Fx_IBM;
                local_F_FSI_diff[3 * (p - p_start) + 1] = 2.0 * m * (particles[p].v - vf[p]) + particles[p].Fy_IBM;
                local_F_FSI_diff[3 * (p - p_start) + 2] = 2.0 * m * (particles[p].w - wf[p]) + particles[p].Fz_IBM;
            }
            else
            {
                local_F_FSI_diff[3 * (p - p_start)] = 2.0 * rho_0 * (particles[p].u - uf[p]);
                local_F_FSI_diff[3 * (p - p_start) + 1] = 2.0 * rho_0 * (particles[p].v - vf[p]);
                local_F_FSI_diff[3 * (p - p_start) + 2] = 2.0 * rho_0 * (particles[p].w - wf[p]);
            }
        }

        // Communicate
        MPI_Allgatherv(local_F_FSI_diff, p_proc, data_type, total_F_FSI_diff, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

        // Setting values
        for (int p = 0; p < n_particles; p++)
        {
            particles[p].Fx_FSI_diff = total_F_FSI_diff[p * 3];
            particles[p].Fy_FSI_diff = total_F_FSI_diff[p * 3 + 1];
            particles[p].Fz_FSI_diff = total_F_FSI_diff[p * 3 + 2];

            if (particles[p].chi)
            {
                particles[p].Fx_IBM -= total_F_FSI_diff[p * 3];
                particles[p].Fy_IBM -= total_F_FSI_diff[p * 3 + 1];
                particles[p].Fz_IBM -= total_F_FSI_diff[p * 3 + 2];
            }
        }

        spread_FSI_forces_to_fluid(sim);
        shift_fluid_velocity_FSI(sim);

        iters++;
    }

    if (params->process_rank == 0)
    {
        sprintf(print_string, "\n\033[2;37m    > Number of iterations: %d\033[0m", iters);
        printf("%-82s", print_string);
    }

    if (iters == max_iters)
    {
        if (params->process_rank == 0)
        {
            printf("\n\nFSI solver did not converge!\n");
            fflush(stdout);
        }
        output_data(sim);
        MPI_Finalize();
        exit(1);
    }

    MPI_Type_free(&data_type);
    free(local_F_FSI_diff);
    free(total_F_FSI_diff);

    free(uf);
    free(vf);
    free(wf);
}

void compute_particle_final_velocity(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double m = params->m;

    double *local_particle_final_velocities = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_particle_final_velocities = (double *)malloc(n_particles * 3 * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        local_particle_final_velocities[(p - p_start) * 3] = particles[p].u + 0.5 / m * particles[p].Fx_IBM;
        local_particle_final_velocities[(p - p_start) * 3 + 1] = particles[p].v + 0.5 / m * particles[p].Fy_IBM;
        local_particle_final_velocities[(p - p_start) * 3 + 2] = particles[p].w + 0.5 / m * particles[p].Fz_IBM;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_final_velocities, p_proc, data_type, total_particle_final_velocities, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].u = total_particle_final_velocities[p * 3];
        particles[p].v = total_particle_final_velocities[p * 3 + 1];
        particles[p].w = total_particle_final_velocities[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_final_velocities);
    free(local_particle_final_velocities);
}

void compute_particle_total_force(SimulationBag *sim)
{
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    LSMParticle *particles = lsm->particles;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;
    int p_proc = lsm->p_proc;

    double *local_particle_total_forces = (double *)malloc(p_proc * 3 * sizeof(double));
    double *total_particle_total_forces = (double *)malloc(n_particles * 3 * sizeof(double));

    for (int p = p_start; p < p_end; p++)
    {
        local_particle_total_forces[(p - p_start) * 3] = particles[p].Fx_spring + particles[p].Fx_damp + particles[p].Fx_IBM + particles[p].Fx_ext;
        local_particle_total_forces[(p - p_start) * 3 + 1] = particles[p].Fy_spring + particles[p].Fy_damp + particles[p].Fy_IBM + particles[p].Fy_ext;
        local_particle_total_forces[(p - p_start) * 3 + 2] = particles[p].Fz_spring + particles[p].Fz_damp + particles[p].Fz_IBM + particles[p].Fz_ext;
    }

    // Communicate
    MPI_Datatype data_type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &data_type);
    MPI_Type_commit(&data_type);

    MPI_Allgatherv(local_particle_total_forces, p_proc, data_type, total_particle_total_forces, lsm->p_recvcounts, lsm->p_displs, data_type, params->comm_xslices);

    // Setting values
    for (int p = 0; p < n_particles; p++)
    {
        particles[p].Fx = total_particle_total_forces[p * 3];
        particles[p].Fy = total_particle_total_forces[p * 3 + 1];
        particles[p].Fz = total_particle_total_forces[p * 3 + 2];
    }

    MPI_Type_free(&data_type);
    free(total_particle_total_forces);
    free(local_particle_total_forces);
}