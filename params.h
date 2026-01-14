#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 250000;
    params->NSTORE = 2500;
    params->NLOG = 250;
    params->NX = 101;
    params->NY = 65;
    params->NZ = 1;

    // Relaxation time
    params->tau = 1.0;

    // Starting density
    params->rho_0 = 1.0;

    params->gx = 1e-6;
    params->gy = 0.0;
    params->gz = 0.0;
}

#endif