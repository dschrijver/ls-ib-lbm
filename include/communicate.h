#ifndef COMMUNICATE_H
#define COMMUNICATE_H

#include "datatypes.h"

void communicate_dists(SimulationBag *sim);
void communicate_dist(double *dist, double *send_buffer, double *recv_buffer, int tag, SimulationBag *sim);

#endif