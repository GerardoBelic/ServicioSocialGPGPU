#include <iostream>
#include <vector>

#include <cmath>
#include <thread>

#include "CLI/CLI.hpp"

#include "Particle_Set.hpp"
#include "Computation_Info.hpp"
#include "N_Body_Compute.hpp"


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

    //CLI11_PARSE(app, argc, argv);
    app.parse(argc, argv);

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

	auto [particle_set, computation_info] = get_computation_info(argc, argv);//std::cout << particle_set.positions[0].x << ' ' << particle_set.positions[0].y << ' ' << particle_set.positions[0].z << '\n';

    n_body_computation_dispatcher(particle_set, computation_info);//std::cout << particle_set.positions[0].x << ' ' << particle_set.positions[0].y << ' ' << particle_set.positions[0].z << '\n';

}
