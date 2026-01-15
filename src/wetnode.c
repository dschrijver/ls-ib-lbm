#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/utils.h"
#include "../include/forcing.h"
#include "../include/wetnode.h"

#ifdef WETNODE
void wetnode_boundary_conditions(SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    // Mass conservation after streaming
#if defined(LEFT_NEBB_VELOCITY)
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_mass_conservation_streaming(0, j, k, 1, 0, 0, sim);
            }
        }
    }
#endif

#if defined(RIGHT_NEBB_VELOCITY)
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_mass_conservation_streaming(params->NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#if defined(BOTTOM_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#if defined(TOP_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#if defined(BACK_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#if defined(FRONT_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

    // Set velocities
#ifdef LEFT_NEBB_VELOCITY
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->fields->u[INDEX(0, j, k)] = LEFT_U_VELOCITY;
                sim->fields->v[INDEX(0, j, k)] = LEFT_V_VELOCITY;
                sim->fields->w[INDEX(0, j, k)] = LEFT_W_VELOCITY;
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_VELOCITY
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->fields->u[INDEX(params->NX - 1, j, k)] = RIGHT_U_VELOCITY;
                sim->fields->v[INDEX(params->NX - 1, j, k)] = RIGHT_V_VELOCITY;
                sim->fields->w[INDEX(params->NX - 1, j, k)] = RIGHT_W_VELOCITY;
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->fields->u[INDEX(i, 0, k)] = BOTTOM_U_VELOCITY;
            sim->fields->v[INDEX(i, 0, k)] = BOTTOM_V_VELOCITY;
            sim->fields->w[INDEX(i, 0, k)] = BOTTOM_W_VELOCITY;
        }
    }
#endif

#ifdef TOP_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->fields->u[INDEX(i, NY - 1, k)] = TOP_U_VELOCITY;
            sim->fields->v[INDEX(i, NY - 1, k)] = TOP_V_VELOCITY;
            sim->fields->w[INDEX(i, NY - 1, k)] = TOP_W_VELOCITY;
        }
    }
#endif

#ifdef BACK_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->fields->u[INDEX(i, j, 0)] = BACK_U_VELOCITY;
            sim->fields->v[INDEX(i, j, 0)] = BACK_V_VELOCITY;
            sim->fields->w[INDEX(i, j, 0)] = BACK_W_VELOCITY;
        }
    }
#endif

#ifdef FRONT_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->fields->u[INDEX(i, j, NZ - 1)] = FRONT_U_VELOCITY;
            sim->fields->v[INDEX(i, j, NZ - 1)] = FRONT_V_VELOCITY;
            sim->fields->w[INDEX(i, j, NZ - 1)] = FRONT_W_VELOCITY;
        }
    }
#endif

    // Compute densities
#ifdef LEFT_NEBB_VELOCITY
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_density(0, j, k, 1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_VELOCITY
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_density(params->NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#ifdef TOP_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#ifdef BACK_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_density(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_density(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

    // Evaluate forces
    evaluate_gravity_force(sim);
    evaluate_total_force(sim);

    // Non-equilibrium bounce-back
#if defined(LEFT_NEBB_VELOCITY)
    if (i_start == 0)
    {
        non_equilibrium_bounce_back_x(0, 1, sim);
    }
#endif

#if defined(RIGHT_NEBB_VELOCITY)
    if (i_end == params->NX)
    {
        non_equilibrium_bounce_back_x(params->NX - 1, -1, sim);
    }
#endif

#if defined(BOTTOM_NEBB_VELOCITY)
    non_equilibrium_bounce_back_y(0, 1, sim);
#endif

#if defined(TOP_NEBB_VELOCITY)
    non_equilibrium_bounce_back_y(NY - 1, -1, sim);
#endif

#if defined(BACK_NEBB_VELOCITY)
    non_equilibrium_bounce_back_z(0, 1, sim);
#endif

#if defined(FRONT_NEBB_VELOCITY)
    non_equilibrium_bounce_back_z(NZ - 1, -1, sim);
#endif

}

void wetnode_mass_conservation_streaming(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;

    double cn;

    int i_start = params->i_start;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    for (int p = 0; p < NP; p++)
    {
        cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;
        if (cn < 0)
        {
            f1[INDEX_F(i, j, k, 0)] += f2[INDEX_F(i, j, k, p)] - f1[INDEX_F(i, j, k, p)];
        }
    }
}

void wetnode_compute_density(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i;
    double un;
    int cn;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    int i_start = params->i_start;

    double *rho = fields->rho;

    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    rho_i = f2[INDEX_F(i, j, k, 0)];

    for (int p = 1; p < NP; p++)
    {
        cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

        if (cn < 0)
        {
            rho_i += f1[INDEX_F(i, j, k, p)] + f2[INDEX_F(i, j, k, p)];
        }
        else if (cn == 0)
        {
            rho_i += f1[INDEX_F(i, j, k, p)];
        }
    }

    un = u[INDEX(i, j, k)] * (double)nx + v[INDEX(i, j, k)] * (double)ny + w[INDEX(i, j, k)] * (double)nz;

    rho[INDEX(i, j, k)] = rho_i / (1.0 - un);
}

void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;

    double *rho = fields->rho;
    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *f1 = dists->f1;

    double C_norm = stencil->C_norm;
    double C_par = stencil->C_par;

    for (int j = 0; j < NY; j++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_i = rho[INDEX(i, j, k)];

            if (rho_i == 0.0)
            {
                for (int p = 0; p < NP; p++)
                {
                    f1[INDEX_F(i, j, k, p)] = 0.0;
                }
                continue;
            }

            u_i = u[INDEX(i, j, k)];
            v_i = v[INDEX(i, j, k)];
            w_i = w[INDEX(i, j, k)];

            Fx_i = Fx[INDEX(i, j, k)];
            Fy_i = Fy[INDEX(i, j, k)];
            Fz_i = Fz[INDEX(i, j, k)];

            Nx = 0.5 * Fx_i / C_norm;
            Ny = 0.5 * Fy_i - rho_i * v_i;
            Nz = 0.5 * Fz_i - rho_i * w_i;

            for (int p = 1; p < NP; p++)
            {
                if (cx[p] == 0)
                {
                    Ny += f1[INDEX_F(i, j, k, p)] * (double)cy[p];
                    Nz += f1[INDEX_F(i, j, k, p)] * (double)cz[p];
                }

                else if (cx[p] * nx > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Ny += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cy[p];
                    Nz += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cz[p];
                }
            }

            Ny /= C_par;
            Nz /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cx[p] * nx > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                    f1[INDEX_F(i, j, k, p)] = f1[INDEX_F(i, j, k, p_bounceback[p])] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                }
            }

            f1[INDEX_F(i, j, k, 0)] += 0.5 * Fx_i * nx;
        }
    }
}

void non_equilibrium_bounce_back_y(int j, int ny, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = fields->rho;
    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *f1 = dists->f1;

    double C_norm = stencil->C_norm;
    double C_par = stencil->C_par;

    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_i = rho[INDEX(i, j, k)];

            if (rho_i == 0.0)
            {
                for (int p = 0; p < NP; p++)
                {
                    f1[INDEX_F(i, j, k, p)] = 0.0;
                }
                continue;
            }

            u_i = u[INDEX(i, j, k)];
            v_i = v[INDEX(i, j, k)];
            w_i = w[INDEX(i, j, k)];

            Fx_i = Fx[INDEX(i, j, k)];
            Fy_i = Fy[INDEX(i, j, k)];
            Fz_i = Fz[INDEX(i, j, k)];

            Nx = 0.5 * Fx_i - rho_i * u_i;
            Ny = 0.5 * Fy_i / C_norm;
            Nz = 0.5 * Fz_i - rho_i * w_i;

            for (int p = 1; p < NP; p++)
            {
                if (cy[p] == 0)
                {
                    Nx += f1[INDEX_F(i, j, k, p)] * (double)cx[p];
                    Nz += f1[INDEX_F(i, j, k, p)] * (double)cz[p];
                }

                else if (cy[p] * ny > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Nx += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cx[p];
                    Nz += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cz[p];
                }
            }

            Nx /= C_par;
            Nz /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cy[p] * ny > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                    f1[INDEX_F(i, j, k, p)] = f1[INDEX_F(i, j, k, p_bounceback[p])] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                }
            }

            f1[INDEX_F(i, j, k, 0)] += 0.5 * Fy_i * ny;
        }
    }
}

void non_equilibrium_bounce_back_z(int k, int nz, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    FieldBag *fields = sim->fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = fields->rho;
    double *Fx = fields->Fx;
    double *Fy = fields->Fy;
    double *Fz = fields->Fz;

    double *u = fields->u;
    double *v = fields->v;
    double *w = fields->w;

    double *f1 = dists->f1;

    double C_norm = stencil->C_norm;
    double C_par = stencil->C_par;

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            rho_i = rho[INDEX(i, j, k)];

            if (rho_i == 0.0)
            {
                for (int p = 0; p < NP; p++)
                {
                    f1[INDEX_F(i, j, k, p)] = 0.0;
                }
                continue;
            }

            u_i = u[INDEX(i, j, k)];
            v_i = v[INDEX(i, j, k)];
            w_i = w[INDEX(i, j, k)];

            Fx_i = Fx[INDEX(i, j, k)];
            Fy_i = Fy[INDEX(i, j, k)];
            Fz_i = Fz[INDEX(i, j, k)];

            Nx = 0.5 * Fx_i - rho_i * u_i;
            Ny = 0.5 * Fy_i - rho_i * v_i;
            Nz = 0.5 * Fz_i / C_norm;

            for (int p = 1; p < NP; p++)
            {
                if (cz[p] == 0)
                {
                    Nx += f1[INDEX_F(i, j, k, p)] * (double)cx[p];
                    Ny += f1[INDEX_F(i, j, k, p)] * (double)cy[p];
                }

                else if (cz[p] * nz > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Nx += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cx[p];
                    Ny += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cy[p];
                }
            }

            Nx /= C_par;
            Ny /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cz[p] * nz > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                    f1[INDEX_F(i, j, k, p)] = f1[INDEX_F(i, j, k, p_bounceback[p])] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                }
            }

            f1[INDEX_F(i, j, k, 0)] += 0.5 * Fz_i * nz;
        }
    }
}

#endif