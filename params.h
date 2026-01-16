#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 250000;
    params->NSTORE = 2500;
    params->NLOG = 250;
    params->NX = 100;
    params->NY = 20;
    params->NZ = 1;

    // Relaxation time
    params->tau = 1.0;

    // Starting density
    params->rho_0 = 1.0;

    // Gravity 
    params->gx = 1e-6;
    params->gy = 0.0;
    params->gz = 0.0;

    // LSM
    params->N_connections_bulk = 18;
}

#endif