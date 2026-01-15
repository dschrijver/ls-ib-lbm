#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/collide.h"

void collide_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    FieldBag *fields = sim->fields;
    Stencil *stencil = sim->stencil;

    double rho_i;
    double u_i, v_i, w_i, u2, uc;
    double feq, S;
    double Fx_i, Fy_i, Fz_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double tau = params->tau;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *rho = fields->rho;
    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX(i, j, k)];
        u_i = u[INDEX(i, j, k)];
        v_i = v[INDEX(i, j, k)];
        w_i = w[INDEX(i, j, k)];
        u2 = u_i * u_i + v_i * v_i + w_i * w_i;

        rho_i = rho[INDEX(i, j, k)];
        Fx_i = Fx[INDEX(i, j, k)];
        Fy_i = Fy[INDEX(i, j, k)];
        Fz_i = Fz[INDEX(i, j, k)];

        for (int p = 0; p < NP; p++)
        {
            uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

            feq = rho_i * wp[p] * (1.0 + uc / cs2 + (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
            
            S = wp[p] * ((((double)cx[p] - u_i) / cs2 + uc / (cs2 * cs2) * (double)cx[p]) * Fx_i + (((double)cy[p] - v_i) / cs2 + uc / (cs2 * cs2) * (double)cy[p]) * Fy_i + (((double)cz[p] - w_i) / cs2 + uc / (cs2 * cs2) * (double)cz[p]) * Fz_i);
            f2[INDEX_F(i, j, k, p)] = f1[INDEX_F(i, j, k, p)] - 1.0 / tau * (f1[INDEX_F(i, j, k, p)] - feq) + (1.0 - 1.0 / (2.0 * tau)) * S;
        }
    }
}