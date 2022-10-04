#include "N_Body_Compute.hpp"

#include <iostream>
#include <omp.h>

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "Particle_Set.hpp"
#include "Computation_Info.hpp"

auto n_body_velocity_calc(unsigned index, Particle_Set &particles, const float deltaTime) -> void
{
	// Aliases for our variables
	std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;
    const float mass = particles.mass;
    const unsigned int numParticles = particles.numParticles;

    unsigned i = index;

    /// Gravitational constant
    const float G = 6.6743e-11f;

    glm::vec3 force(0.0f, 0.0f, 0.0f);

    /// Loop to calculate a particle force due to all the other particles
    for (unsigned j = 0; j < numParticles; ++j)
    {
        /// If indexes i & j refer to the same body, do nothing
        if (i == j)
            continue;

        float distance2 = glm::distance2(positions[i], positions[j]);
        glm::vec3 direction = glm::normalize(positions[j] - positions[i]);

        force += G * mass * mass * direction / distance2;
    }

    /// f=ma => a=f/m
    glm::vec3 acceleration = force / mass;

    /// v=v0+a*dt
    velocities[i] += acceleration * deltaTime;


}

auto n_body_position_calc(unsigned index, Particle_Set &particles, const float deltaTime) -> void
{
	/// Aliases for our variables
	std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;

    unsigned i = index;

    positions[i] += velocities[i] * deltaTime;

}

auto n_body_computation_dispatcher(Particle_Set &particles, const Computation_Info &info) -> void
{
    /// Loop for every iteration
	for (unsigned iter = 0; iter < info.numIterations; ++iter)
    {

        /** 'pragma omp parallel for' indicates that the next for loop is going to
            split the work equally in 'num_threads' (max_threads if not specified).
            First we start with the velocity calculations and then do the positions.
        */
        #pragma omp parallel for num_threads(info.numThreads)
        for (unsigned p = 0; p < particles.numParticles; ++p)
        {
            n_body_velocity_calc(p, particles, info.timeStep);
        }

        /** By default OpenMP sets a barrier after the parallel for loop and waits
            for all threads to finish their work.
        */

        #pragma omp parallel for num_threads(info.numThreads)
        for (unsigned p = 0; p < particles.numParticles; ++p)
        {
            n_body_position_calc(p, particles, info.timeStep);
        }

    }

}
