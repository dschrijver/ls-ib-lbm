#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/forcing.h"

void evaluate_gravity_force(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    double rho_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double gx = params->gx;
    double gy = params->gy;
    double gz = params->gz;

    double *rho = fields->rho;
    double *Fx_grav = fields->Fx_grav;
    double *Fy_grav = fields->Fy_grav;
    double *Fz_grav = fields->Fz_grav;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX(i, j, k)];

        Fx_grav[INDEX(i, j, k)] = rho_i * gx;
        Fy_grav[INDEX(i, j, k)] = rho_i * gy;
        Fz_grav[INDEX(i, j, k)] = rho_i * gz;
    }
}

void evaluate_total_force(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

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

    FOR_DOMAIN
    {
        Fx[INDEX(i, j, k)] = Fx_grav[INDEX(i, j, k)] + Fx_IBM[INDEX(i, j, k)] + Fx_rigid[INDEX(i, j, k)];
        Fy[INDEX(i, j, k)] = Fy_grav[INDEX(i, j, k)] + Fy_IBM[INDEX(i, j, k)] + Fy_rigid[INDEX(i, j, k)];
        Fz[INDEX(i, j, k)] = Fz_grav[INDEX(i, j, k)] + Fz_IBM[INDEX(i, j, k)] + Fz_rigid[INDEX(i, j, k)];
    }
}