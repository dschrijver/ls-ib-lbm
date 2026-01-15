#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../definitions.h"
#include "../include/stream.h"

void stream_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            ic = i - cx[p];
            jc = j - cy[p];
            kc = k - cz[p];

#ifdef YPERIODIC
            jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
            kc = mod(kc, NZ);
#endif

#ifdef LEFT_BOUNCEBACK_NOSLIP
            if (ic < 0)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

#ifdef RIGHT_BOUNCEBACK_NOSLIP
            if (ic > params->NX - 1)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

#ifdef BOTTOM_BOUNCEBACK_NOSLIP
            if (jc < 0)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

#ifdef TOP_BOUNCEBACK_NOSLIP
            if (jc > NY - 1)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

#ifdef BACK_BOUNCEBACK_NOSLIP
            if (kc < 0)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

#ifdef FRONT_BOUNCEBACK_NOSLIP
            if (kc > NZ - 1)
            {
                f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p])];
                continue;
            }
#endif

            f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(ic, jc, kc, p)];
        }
    }
}