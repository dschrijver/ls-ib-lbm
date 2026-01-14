#ifndef STENCIL_H
#define STENCIL_H

#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "memory.h"

static inline void initialize_stencil(SimulationBag *sim)
{
    Stencil *stencil = sim->stencil;

    stencil->NP = 19;

    allocate_stencil(sim);

    stencil->cs2 = 1.0 / 3.0;

    int cx[] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0};
    int cy[] = {0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1};
    int cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1};
    double wp[] = {1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    int p_bounceback[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};

    memcpy(stencil->cx, cx, stencil->NP * sizeof(int));
    memcpy(stencil->cy, cy, stencil->NP * sizeof(int));
    memcpy(stencil->cz, cz, stencil->NP * sizeof(int));
    memcpy(stencil->wp, wp, stencil->NP * sizeof(double));
    memcpy(stencil->p_bounceback, p_bounceback, stencil->NP * sizeof(int));

    stencil->C_norm = 1.0 + 2.0 * sqrt(2.0);
    stencil->C_par = sqrt(2.0);
}

#endif