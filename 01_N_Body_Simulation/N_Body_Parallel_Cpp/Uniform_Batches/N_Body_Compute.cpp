#include "N_Body_Compute.hpp"
#include <iostream>
#include <thread>
#include <boost/thread/barrier.hpp>

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

auto n_body_thread(unsigned begin, unsigned end, Particle_Set &particles, const Computation_Info &info, boost::barrier &sync_point) -> void
{

    for (unsigned i = 0; i < info.numIterations; ++i)
    {
        /// Dispatch the velocity calculation
        for (unsigned j = begin; j < end; ++j)
        {
            n_body_velocity_calc(j, particles, info.timeStep);
        }

        /// Sync threads for memory coherence
        sync_point.wait();

        /// Dispatch the position calculation
        for (unsigned j = begin; j < end; ++j)
        {
            n_body_position_calc(j, particles, info.timeStep);
        }

        /// Sync again
        sync_point.wait();
    }

    return;

}

auto n_body_computation_dispatcher(Particle_Set &particles, const Computation_Info &info) -> void
{

    /// Thread dispatch for updating particles speed
	unsigned step_index = particles.numParticles / info.numThreads;
	unsigned remainder = particles.numParticles % info.numThreads;

	/// Barrier for syncing all threads
	boost::barrier sync_point{info.numThreads};

	std::vector<std::jthread> threads;

	for (unsigned curr_index = 0; curr_index < particles.numParticles; curr_index += step_index)
    {
        if (remainder != 0)
        {
            threads.emplace_back(n_body_thread, curr_index, curr_index + step_index + 1, std::ref(particles), std::cref(info), std::ref(sync_point));
            ++curr_index;
            --remainder;
        }
        else
            threads.emplace_back(n_body_thread, curr_index, curr_index + step_index, std::ref(particles), std::cref(info), std::ref(sync_point));
    }

}
