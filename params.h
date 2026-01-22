#ifndef PARAMS_H
#define PARAMS_H

#include <math.h>

#include "include/datatypes.h"
#include "include/utils.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 100;
    params->NSTORE = 1;
    params->NLOG = 1;
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
    params->N_connections_bulk = 8;
    params->r_max = sqrt((double)(params->N_connections_bulk + 1)/DS_PI); // 2D
    // params->r_max = pow((double)(params->N_connections_bulk + 1) / (4.0 / 3.0 * DS_PI), 1.0 / 3.0); // 3D
    params->m = 1.0;
    params->k = 1e-3;
    params->c = 5e-2;
    params->tol = 1e-6;
    params->max_iters = 1000;
}

#endif