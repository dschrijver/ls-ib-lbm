#ifndef LSM_H
#define LSM_H

#include "datatypes.h"

void create_particle(double x, double y, double z, double u, double v, double w, LSMBag *lsm);
void create_spring(int p1, int p2, double l_eq, LSMBag *lsm);
void add_springs(SimulationBag *sim);
void compute_chi(SimulationBag *sim);

void update_particle_positions(SimulationBag *sim);
void update_particle_predicted_velocities(SimulationBag *sim);

void evaluate_particle_external_forces(SimulationBag *sim);
void evaluate_spring_forces(SimulationBag *sim);
void evaluate_damping_forces(SimulationBag *sim);
void update_particle_preliminary_velocities(SimulationBag *sim);

void interpolate_fluid_velocity_to_particle(double *uf, double *vf, double *wf, SimulationBag *sim);
void spread_FSI_forces_to_fluid(SimulationBag *sim);
void shift_fluid_velocity_FSI(SimulationBag *sim);
void evaluate_FSI_forces(SimulationBag *sim);

void compute_particle_final_velocity(SimulationBag *sim);
void compute_particle_total_force(SimulationBag *sim);

#endif