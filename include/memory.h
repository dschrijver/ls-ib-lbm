#ifndef MEMORY_H
#define MEMORY_H

#include "datatypes.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil, LSMBag **lsm);
void allocate_stencil(SimulationBag *sim);
void allocate_distributions(SimulationBag *sim);
void allocate_fields(SimulationBag *sim);

#endif