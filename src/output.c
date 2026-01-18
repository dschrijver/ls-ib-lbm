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
    
    double *Fx_IBM = fields->Fx_IBM;
    double *Fy_IBM = fields->Fy_IBM;
    double *Fz_IBM = fields->Fz_IBM;

    double *Fx_rigid = fields->Fx_rigid;
    double *Fy_rigid = fields->Fy_rigid;
    double *Fz_rigid = fields->Fz_rigid;

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

    output_field(Fx_IBM, "Fx_IBM", file_id, sim);
    output_field(Fy_IBM, "Fy_IBM", file_id, sim);
    output_field(Fz_IBM, "Fz_IBM", file_id, sim);

    output_field(Fx_rigid, "Fx_rigid", file_id, sim);
    output_field(Fy_rigid, "Fy_rigid", file_id, sim);
    output_field(Fz_rigid, "Fz_rigid", file_id, sim);

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

    int n_particles = lsm->n_particles;
    int n_springs = lsm->n_springs;

    int p_start = lsm->p_start;
    int p_end = lsm->p_end;

    int s_start = lsm->s_start;
    int s_end = lsm->s_end;

    hid_t particle_type = H5Tcreate(H5T_COMPOUND, sizeof(LSMParticle));

    H5Tinsert(particle_type, "x", HOFFSET(LSMParticle, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "y", HOFFSET(LSMParticle, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "z", HOFFSET(LSMParticle, z), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "u", HOFFSET(LSMParticle, u), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "v", HOFFSET(LSMParticle, v), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "w", HOFFSET(LSMParticle, w), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "u_predict", HOFFSET(LSMParticle, u_predict), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "v_predict", HOFFSET(LSMParticle, v_predict), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "w_predict", HOFFSET(LSMParticle, w_predict), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "u_fluid", HOFFSET(LSMParticle, u_fluid), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "v_fluid", HOFFSET(LSMParticle, v_fluid), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "w_fluid", HOFFSET(LSMParticle, w_fluid), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "Fx_ext", HOFFSET(LSMParticle, Fx_ext), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fy_ext", HOFFSET(LSMParticle, Fy_ext), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fz_ext", HOFFSET(LSMParticle, Fz_ext), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "Fx_spring", HOFFSET(LSMParticle, Fx_spring), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fy_spring", HOFFSET(LSMParticle, Fy_spring), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fz_spring", HOFFSET(LSMParticle, Fz_spring), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "Fx_damp", HOFFSET(LSMParticle, Fx_damp), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fy_damp", HOFFSET(LSMParticle, Fy_damp), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fz_damp", HOFFSET(LSMParticle, Fz_damp), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "Fx_IBM", HOFFSET(LSMParticle, Fx_IBM), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fy_IBM", HOFFSET(LSMParticle, Fy_IBM), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fz_IBM", HOFFSET(LSMParticle, Fz_IBM), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "Fx", HOFFSET(LSMParticle, Fx), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fy", HOFFSET(LSMParticle, Fy), H5T_NATIVE_DOUBLE);
    H5Tinsert(particle_type, "Fz", HOFFSET(LSMParticle, Fz), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "weight", HOFFSET(LSMParticle, weight), H5T_NATIVE_DOUBLE);

    H5Tinsert(particle_type, "N_connections", HOFFSET(LSMParticle, N_connections), H5T_NATIVE_INT);
    H5Tinsert(particle_type, "chi", HOFFSET(LSMParticle, chi), H5T_NATIVE_DOUBLE);

    hid_t spring_type = H5Tcreate(H5T_COMPOUND, sizeof(LSMSpring));

    H5Tinsert(spring_type, "p1", HOFFSET(LSMSpring, p1), H5T_NATIVE_INT);
    H5Tinsert(spring_type, "p2", HOFFSET(LSMSpring, p2), H5T_NATIVE_INT);

    H5Tinsert(spring_type, "l_eq", HOFFSET(LSMSpring, l_eq), H5T_NATIVE_DOUBLE);
    H5Tinsert(spring_type, "S", HOFFSET(LSMSpring, S), H5T_NATIVE_DOUBLE);

    H5Tinsert(spring_type, "active", HOFFSET(LSMSpring, active), H5T_NATIVE_INT);
    
    output_lsm_field(&lsm->particles[p_start], "Particles", particle_type, p_start, p_end, n_particles, loc_id);
    output_lsm_field(&lsm->springs[s_start], "Springs", spring_type, s_start, s_end, n_springs, loc_id);

    H5Tclose(particle_type);
    H5Tclose(spring_type);
}

void output_lsm_field(void* field, char *fieldname, hid_t field_type, int start, int end, int total, hid_t loc_id) 
{
     // Space occupied in file
    hsize_t dims_file[1] = {total};
    hid_t filespace = H5Screate_simple(1, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[1] = {end-start};
    hid_t memspace  = H5Screate_simple(1, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, field_type, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t file_start[1] = {start};
    hsize_t count[1] = {end-start};
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, field_type, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Sclose(filespace);
}