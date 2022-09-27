#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <thread>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "CLI/CLI.hpp"

struct Particle_Set
{

    /**
    A particle has a position and velocity
	All particles shares mass in this scenario
	Positions are randomly selected with uniform distribution
    */
    Particle_Set(unsigned int _numParticles, float _mass = 1.0f) : positions(_numParticles), velocities(_numParticles), numParticles(_numParticles), mass(_mass)
    {
        for (auto& pos : positions)
        {
            std::random_device rd;
            std::mt19937 e2(rd());
            std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

            pos.x = dist(e2);
            pos.y = dist(e2);
			pos.z = dist(e2);
        }
    }

    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> velocities;

    const unsigned int numParticles;
    const float mass;

};

struct Computation_Info
{
	
	Computation_Info(float _timeStep, unsigned _numIterations, unsigned _numThreads = 1) :
					 timeStep(_timeStep), numIterations(_numIterations), numThreads(_numThreads) { }
	
	float timeStep;
	unsigned numIterations;
	unsigned numThreads;
}

auto n_body_velocity_calc(unsigned begin, unsigned end, Particle_Set &particles, const float deltaTime, ) -> void
{
	// Aliases for our variables
	std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;
    const float mass = particles.mass;
    const unsigned int numParticles = particles.numParticles;

    /// Gravitational constant
    const float G = 6.6743e-11f;

	for (unsigned i = begin; i < end; ++i)
	{
		glm::vec3 force(0.0f, 0.0f, 0.0f);
		
		/// Loop to calculate a particle's force due to all the other particles
		for (unsigned j = 0; j < numParticles; ++j)
		{
			/// If indexes i & j refer to the same body, do nothing
			if (i == j)
				continue;

			float distance2 = glm::distance2(positions[i], positions[j]);
			glm::vec3 direction = glm::normalize(positions[j] - positions[i]);

			force += G * mass * mass * direction / distance2;
			
		}
	}
	
	// f=ma => a=f/m
	glm::vec3 acceleration = force / mass;

	// v=v0+a*dt
	glm::vec3 velocity = velocities[i] + acceleration * deltaTime;
}

auto n_body_position_calc(unsigned begin, unsigned end, Particle_Set &particles, const float deltaTime) -> void
{
	// Aliases for our variables
	std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;
	
	unsigned int i = particleIndex;
	
	glm::vec3 position = positions[i] + velocities[i] * deltaTime;
	
	return position;
}

auto n_body_computation(Particle_Set &particles, const Computation_Info &info) -> void
{

    /// Thread dispatch for updating particles speed
	unsigned curr_index = 0;
	unsigned step_index = particles.numParticles / info.numThreads;
	
	while (curr_index < particles.numParticles - step_index)
	{
		n_body_velocity_calc(curr_index, curr_index + step_index, particles, info.timeStep);
		curr_index += step_index;
	}
	n_body_velocity_calc(curr_index, particles.numParticles, particles, info.timeStep);
	
	

    //std::cout << positions[0].x << ' ' << positions[0].y << ' ' << positions[0].z << '\n';

}

auto get_computation_info(int argc, char **argv) -> std::tuple<Particle_Set, Computation_Info>
{
	/// Parse arguments to form the particle set and the computation info
	CLI::App app{"N-Body serial simulation"};
	
	unsigned numParticles = 0;
	float mass = 1.0f;
	float dt = 1.0f;
	unsigned iterations = 0;
	unsigned threads = std::thread::hardware_concurrency();
	
    app.add_option("--particles", numParticles, "Particle count") -> required();
	app.add_option("--mass", mass, "Mass of each particle (def=1.0[kg])");
	app.add_option("--dt", dt, "Timestep between each iteration (def=1.0[s])");
	app.add_option("--iterations", iterations, "Number of iterations of the simulation") -> required();
	app.add_option("--threads", threads, "Number of simultaneous threads (def=max_pc_threads)");

    CLI11_PARSE(app, argc, argv);

    std::cout << "Particle count: " << numParticles << " particles" << "\n";
	std::cout << "Particle mass: " << mass << " [kg]" << "\n";
	std::cout << "Timestep: " << dt << " [s]" << "\n";
	std::cout << "Iterations: " << iterations << " steps" << "\n";
	std::cout << "Threads: " << threads << "\n";
	
	/// Form the particle set and the computation info
	Particle_Set my_set{numParticles, mass};
	Computation_Info info{dt, iterations, threads};
	
	return {my_set, info};
	
	
}

int main(int argc, char **argv)
{
	
	auto [particle_set, computation_info] = get_computation_info(argc, argv);

    n_body_computation(particle_set, computation_info);

}
