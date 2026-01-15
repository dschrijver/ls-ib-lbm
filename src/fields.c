#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/fields.h"

void extract_moments(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    FieldBag *fields = sim->fields;
    Stencil *stencil = sim->stencil;

    double rho_i, u_i, v_i, w_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double cs2 = stencil->cs2;

    double *rho = fields->rho;
    double *pressure = fields->pressure;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *f1 = dists->f1;

    FOR_DOMAIN
    {
        rho_i = 0.0;
        u_i = 0.0;
        v_i = 0.0;
        w_i = 0.0;

        for (int p = 0; p < NP; p++)
        {
            rho_i += f1[INDEX_F(i, j, k, p)];
            u_i += f1[INDEX_F(i, j, k, p)] * (double)cx[p];
            v_i += f1[INDEX_F(i, j, k, p)] * (double)cy[p];
            w_i += f1[INDEX_F(i, j, k, p)] * (double)cz[p];
        }

        rho[INDEX(i, j, k)] = rho_i;

        pressure[INDEX(i, j, k)] = rho[INDEX(i, j, k)] * cs2;

        u[INDEX(i, j, k)] = u_i / rho_i;
        v[INDEX(i, j, k)] = v_i / rho_i;
        w[INDEX(i, j, k)] = w_i / rho_i;
    }
}

void update_final_velocity(SimulationBag *sim)
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

    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX(i, j, k)];

        u[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * Fx[INDEX(i, j, k)];
        v[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * Fy[INDEX(i, j, k)];
        w[INDEX(i, j, k)] += 1.0 / (2.0 * rho_i) * Fz[INDEX(i, j, k)];
    }
}