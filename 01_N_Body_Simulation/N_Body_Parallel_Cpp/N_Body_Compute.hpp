#ifndef N_BODY_COMPUTE
#define N_BODY_COMPUTE

#include <boost/thread/barrier.hpp>

#include "Particle_Set.hpp"
#include "Computation_Info.hpp"

auto n_body_velocity_calc(unsigned index, Particle_Set &particles, const float deltaTime) -> void;

auto n_body_position_calc(unsigned index, Particle_Set &particles, const float deltaTime) -> void;

auto n_body_thread(unsigned begin, unsigned end, Particle_Set &particles, const Computation_Info &info, boost::barrier &sync_point) -> void;

auto n_body_computation_dispatcher(Particle_Set &particles, const Computation_Info &info) -> void;

#endif
