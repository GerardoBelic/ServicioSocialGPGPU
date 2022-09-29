#include <iostream>
#include <vector>
#include <random>
#include <cmath>

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
    Particle_Set(unsigned int _numParticles, float _mass = 1.0f) : positions(_numParticles), velocities(_numParticles), numParticle_Set(_numParticles), mass(_mass)
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

    const unsigned int numParticle_Set;
    const float mass;

};

struct Computation_Info
{

	Computation_Info(float _timeStep, unsigned _numIterations) :
					 timeStep(_timeStep), numIterations(_numIterations) { }

	float timeStep;
	unsigned numIterations;
};

auto n_body_iteration(Particle_Set &particles, const float deltaTime) -> void
{
	/// Aliases for our data variables
    std::vector<glm::vec3> &positions = particles.positions;
    std::vector<glm::vec3> &velocities = particles.velocities;
    const float mass = particles.mass;
    const unsigned int numParticles = particles.numParticles;

    /// Gravitational constant
    const float G = 6.6743e-11f;

    /// Ciclo para calcular la nueva velocidad de cada particula calculando su vector de force y aceleración
    for (unsigned int i = 0; i < numParticles; ++i)
    {

        glm::vec3 force(0.0f, 0.0f, 0.0f);

        /// Loop to calculate a particle's force due to all the other particles
        for (unsigned int j = 0; j < numParticles; ++j)
        {

            /// If indexes i & j refer to the same body, do nothing
            if (i == j)
                continue;

            float distance2 = glm::distance2(positions[i], positions[j]);
            glm::vec3 direction = glm::normalize(positions[j] - positions[i]);

            force += G * mass * mass * direction / distance2;
            
        }
        
		// f=ma => a=f/m
        glm::vec3 acceleration = force / mass;

		// v=v0+a*dt
        velocities[i] += acceleration * deltaTime;

    }

    /// Update positions with new velocities
    for (int i = 0; i < numParticles; ++i)
    {

        positions[i] += velocities[i] * deltaTime;

    }

    //std::cout << positions[0].x << ' ' << positions[0].y << ' ' << positions[0].z << '\n';

}

auto parse_arguments(int argc, char **argv) -> Particle_Set
{
	
	/// Parse arguments to form the particle set and the computation info
	CLI::App app{"N-Body serial simulation"};

	unsigned numParticles = 0;
	float mass = 1.0f;
	float dt = 1.0f;
	unsigned iterations = 0;

    app.add_option("--particles", numParticles, "Particle count") -> required();
	app.add_option("--mass", mass, "Mass of each particle (def=1.0[kg])");
	app.add_option("--dt", dt, "Timestep between each iteration (def=1.0[s])");
	app.add_option("--iterations", iterations, "Number of iterations of the simulation") -> required();

    //CLI11_PARSE(app, argc, argv);
    app.parse(argc, argv);

    std::cout << "Particle count: " << numParticles << " particles" << "\n";
	std::cout << "Particle mass: " << mass << " [kg]" << "\n";
	std::cout << "Timestep: " << dt << " [s]" << "\n";
	std::cout << "Iterations: " << iterations << " steps" << "\n";

	/// Form the particle set and the computation info
	Particle_Set my_set{numParticles, mass};
	Computation_Info info{dt, iterations}

	return {my_set, info};
	
	
}

int main(int argc, char **argv)
{
	
	/// Instance of our particle set and parameters for solving the sim
	auto [particle_set, computation_info] = parse_arguments(argc, argv);

    for (unsigned i = 0; i < computation_info.numIterations; ++i)
        n_body_iteration(particle_set, computation_info.timeStep);

}
