#ifndef OUTPUT_H
#define OUTPUT_H

#include <hdf5.h>
#include "datatypes.h"

void output_data(SimulationBag *sim);
void output_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim);

#endif