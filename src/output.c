#include <hdf5.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/output.h"

void output_data(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    char filename[32];

    int t = params->t;

    double *rho = fields->rho;
    double *pressure = fields->pressure;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *Fx_grav = fields->Fx_grav;
    double *Fy_grav = fields->Fy_grav;
    double *Fz_grav = fields->Fz_grav;

    // Create file
    sprintf(filename, "data_%d.h5", params->n_output);
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    // Create groups
    hid_t gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
    hid_t lsm_group_id = H5Gcreate2(file_id, "/LSM", H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    H5Pclose(gcpl_id);

    // Write time
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t dset_scalar = H5Dcreate2(file_id, "t", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    H5Dclose(dset_scalar);
    H5Sclose(scalar_space);

    output_field(rho, "rho", file_id, sim);
    output_field(pressure, "pressure", file_id, sim);
    output_field(u, "u", file_id, sim);
    output_field(v, "v", file_id, sim);
    output_field(w, "w", file_id, sim);

    output_field(Fx, "Fx", file_id, sim);
    output_field(Fy, "Fy", file_id, sim);
    output_field(Fz, "Fz", file_id, sim);

    output_field(Fx_grav, "Fx_grav", file_id, sim);
    output_field(Fy_grav, "Fy_grav", file_id, sim);
    output_field(Fz_grav, "Fz_grav", file_id, sim);

    // LSM group
    if (sim->lsm->n_particles > 0) {
        output_lsm(lsm_group_id, sim);
    }

    // Close groups
    H5Gclose(lsm_group_id);

    // Close file
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    params->n_output++;
}

void output_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int i_start = params->i_start;
    int NX_proc = params->NX_proc;
    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    hid_t filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[3] = {NX_proc + 4, NY, NZ};
    hid_t memspace = H5Screate_simple(3, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t start_file[3] = {i_start, 0, 0};
    hsize_t count[3] = {NX_proc, NY, NZ};
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Process hyperslab
    hsize_t start_proc[3] = {2, 0, 0};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Sclose(filespace);
}

void output_lsm(hid_t loc_id, SimulationBag *sim) {
    LSMBag *lsm = sim->lsm;
    ParamBag *params = sim->params;

    int n_particles = lsm->n_particles;
    int n_springs = lsm->n_springs;

    // Write number of particles
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t dset_scalar = H5Dcreate2(loc_id, "n_particles", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_particles);
    H5Dclose(dset_scalar);
    H5Sclose(scalar_space);

    // Write number of particles
    scalar_space = H5Screate(H5S_SCALAR);
    dset_scalar = H5Dcreate2(loc_id, "n_springs", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_springs);
    H5Dclose(dset_scalar);
    H5Sclose(scalar_space);

    // Create groups
    hid_t gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
    hid_t particle_group_id = H5Gcreate2(loc_id, "Particles", H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    hid_t spring_group_id = H5Gcreate2(loc_id, "Springs", H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    H5Pclose(gcpl_id);

    // Find local particles
    int particle_start = (float)(params->process_coords[0])/(float)params->number_of_processes*n_particles;
    int particle_end = (float)(params->process_coords[0]+1)/(float)params->number_of_processes*n_particles;
    int particles_proc = particle_end - particle_start;

    int malloc_size = particles_proc*3*sizeof(double);

    double *r = (double *)malloc(malloc_size);
    double *v = (double *)malloc(malloc_size);
    double *F = (double *)malloc(malloc_size);
    double *F_ext = (double *)malloc(malloc_size);
    double *F_spring = (double *)malloc(malloc_size);
    double *F_damp = (double *)malloc(malloc_size);
    double *F_FSI = (double *)malloc(malloc_size);

    int *N_connections = (int *)malloc(particles_proc*sizeof(int));
    double *chi = (double *)malloc(particles_proc*sizeof(double));

    int p = 0;
    int p_local = 0;
    LSMParticle *particle = lsm->particle_first;
    while (particle != NULL) {
        particle->id = p;

        if ((particle->id >= particle_start) && (particle->id < particle_end)) {

            r[3*p_local] = particle->x;
            r[3*p_local + 1] = particle->y;
            r[3*p_local + 2] = particle->z;

            v[3*p_local] = particle->u;
            v[3*p_local + 1] = particle->v;
            v[3*p_local + 2] = particle->w;

            F[3*p_local] = particle->Fx;
            F[3*p_local + 1] = particle->Fy;
            F[3*p_local + 2] = particle->Fz;

            F_ext[3*p_local] = particle->Fx_ext;
            F_ext[3*p_local + 1] = particle->Fy_ext;
            F_ext[3*p_local + 2] = particle->Fz_ext;

            F_spring[3*p_local] = particle->Fx_spring;
            F_spring[3*p_local + 1] = particle->Fy_spring;
            F_spring[3*p_local + 2] = particle->Fz_spring;

            F_damp[3*p_local] = particle->Fx_damp;
            F_damp[3*p_local + 1] = particle->Fy_damp;
            F_damp[3*p_local + 2] = particle->Fz_damp;

            F_FSI[3*p_local] = particle->Fx_FSI;
            F_FSI[3*p_local + 1] = particle->Fy_FSI;
            F_FSI[3*p_local + 2] = particle->Fz_FSI;
            
            N_connections[p_local] = particle->N_connections;
            chi[p_local] = particle->chi;

            p_local++;
        }

        p++;
        particle = particle->next;
    }

    int spring_start = (float)(params->process_coords[0])/(float)params->number_of_processes*n_springs;
    int spring_end = (float)(params->process_coords[0]+1)/(float)params->number_of_processes*n_springs;
    int springs_proc = spring_end - spring_start;

    int *particle_1 = (int *)malloc(springs_proc*sizeof(int));
    int *particle_2 = (int *)malloc(springs_proc*sizeof(int));

    double *l_eq = (double *)malloc(springs_proc*sizeof(double));
    double *S = (double *)malloc(springs_proc*sizeof(double));

    int *active = (int *)malloc(springs_proc*sizeof(int));

    int s = 0;
    int s_local = 0;
    LSMSpring *spring = lsm->spring_first;
    while (spring != NULL) {

        if ((s >= spring_start) && (s < spring_end)) {
            particle_1[s_local] = spring->particle_1->id;
            particle_2[s_local] = spring->particle_2->id;

            l_eq[s_local] = spring->l_eq;
            S[s_local] = spring->S;

            active[s_local] = spring->active;

            s_local++;
        }

        s++;
        spring = spring->next;
    }

    output_lsm_field(r, "r", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(v, "v", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(F, "F", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(F_ext, "F_ext", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(F_spring, "F_spring", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(F_damp, "F_damp", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field(F_FSI, "F_FSI", particle_start, particle_end, particle_group_id, sim);
    output_lsm_field_1d(N_connections, "N_connections", particle_start, particle_end, n_particles, H5T_NATIVE_INT, particle_group_id);
    output_lsm_field_1d(chi, "chi", particle_start, particle_end, n_particles, H5T_NATIVE_DOUBLE, particle_group_id);

    output_lsm_field_1d(particle_1, "particle_1", spring_start, spring_end, n_springs, H5T_NATIVE_INT, spring_group_id);
    output_lsm_field_1d(particle_2, "particle_2", spring_start, spring_end, n_springs, H5T_NATIVE_INT, spring_group_id);
    output_lsm_field_1d(l_eq, "l_eq", spring_start, spring_end, n_springs, H5T_NATIVE_DOUBLE, spring_group_id);
    output_lsm_field_1d(S, "S", spring_start, spring_end, n_springs, H5T_NATIVE_DOUBLE, spring_group_id);
    output_lsm_field_1d(active, "active", spring_start, spring_end, n_springs, H5T_NATIVE_INT, spring_group_id);

    // Close groups
    H5Gclose(particle_group_id);
    H5Gclose(spring_group_id);

    free(r);
    free(v);
    free(F);
    free(F_ext);
    free(F_spring);
    free(F_damp);
    free(F_FSI);
    free(N_connections);
    free(chi);

    free(particle_1);
    free(particle_2);
    free(l_eq);
    free(S);
    free(active);
}

void output_lsm_field(double* field, char *fieldname, int start, int end, hid_t loc_id, SimulationBag *sim) {
    LSMBag *lsm = sim->lsm;

    int n_particles = lsm->n_particles;

     // Space occupied in file
    hsize_t dims_file[2] = {n_particles, 3};
    hid_t filespace = H5Screate_simple(2, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[2] = {end-start, 3};
    hid_t memspace  = H5Screate_simple(2, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t start_file[2] = {start, 0};
    hsize_t count[2] = {end-start, 3};   
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Sclose(filespace);
}

void output_lsm_field_1d(void* field, char *fieldname, int start, int end, int total, hid_t data_type, hid_t loc_id) {
     // Space occupied in file
    hsize_t dims_file[1] = {total};
    hid_t filespace = H5Screate_simple(1, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[1] = {end-start};
    hid_t memspace  = H5Screate_simple(1, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, data_type, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t start_file[1] = {start};
    hsize_t count[1] = {end-start};   
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, data_type, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Sclose(filespace);
}