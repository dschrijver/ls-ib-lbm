#ifndef OUTPUT_H
#define OUTPUT_H

#include <hdf5.h>
#include "datatypes.h"

void output_data(SimulationBag *sim);
void output_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim);
void output_lsm(hid_t loc_id, SimulationBag *sim);
void output_lsm_field(void* field, char *fieldname, hid_t field_type, int start, int end, int total, hid_t loc_id);

#endif