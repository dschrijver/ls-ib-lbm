#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "datatypes.h"

void initialize_MPI(ParamBag *params);
void initialize_fields(SimulationBag *sim);
void initialize_distribution(SimulationBag *sim);

#endif