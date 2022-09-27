#ifndef N_BODY_COMPUTE
#define N_BODY_COMPUTE

#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "Particle_Set.hpp"
#include "Computation_Info.hpp"

std::mutex m;
std::condition_variable cv;
unsigned threads_done = 0;

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
	// Aliases for our variables
	std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;

    unsigned i = index;

    positions[i] += velocities[i] * deltaTime;

}

auto n_body_thread(unsigned begin, unsigned end, Particle_Set &particles, const Computation_Info &info) -> void
{

    for (int i = 0; i < info.numIterations; ++i)
    {
        /// Dispatch the velocity calculation
        for (int j = begin; j < end; ++j)
        {
            n_body_velocity_calc(j, particles, info.timeStep);
        }

        /// Sync threads for memory coherence
        std::unique_lock<std::mutex> lk(m);
        ++threads_done;

        if (threads_done == info.numThreads)
        {
            lk.unlock();
            cv.notify_all();
            lk.lock();
        }
        else
            cv.wait(lk, [&info]{ return threads_done == info.numThreads;});

        --threads_done;
        lk.unlock();


        /// Dispatch the position calculation
        for (int j = begin; j < end; ++j)
        {
            n_body_position_calc(j, particles, info.timeStep);
        }

        /// Sync again
        lk.lock();
        ++threads_done;

        if (threads_done == info.numThreads)
        {
            std::cout << particles.positions[0].x << ' ' << particles.positions[0].y << ' ' << particles.positions[0].z << '\n';
            lk.unlock();
            cv.notify_all();
            lk.lock();
        }
        else
            cv.wait(lk, [&info]{ return threads_done == info.numThreads;});

        --threads_done;
        lk.unlock();

    }

    return;

}

auto n_body_computation_dispatcher(Particle_Set &particles, const Computation_Info &info) -> void
{

    /// Thread dispatch for updating particles speed
	unsigned curr_index = 0;
	unsigned step_index = particles.numParticles / info.numThreads;

	while (curr_index < particles.numParticles - step_index)
	{
		std::jthread thr{n_body_thread, curr_index, curr_index + step_index, std::ref(particles), std::cref(info)};
		curr_index += step_index;
	}
	std::jthread thr{n_body_thread, curr_index, particles.numParticles, std::ref(particles), std::cref(info)};

    //std::cout << positions[0].x << ' ' << positions[0].y << ' ' << positions[0].z << '\n';

}

struct N_Body_Compute
{



};

#endif
