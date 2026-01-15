#include <hdf5.h>

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