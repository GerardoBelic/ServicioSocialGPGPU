#include <iostream>
#include <vector>
#include <thread>

#include "CLI/CLI.hpp"

#include "Particle_Set.hpp"
#include "Computation_Info.hpp"
#include "N_Body_Compute.hpp"


auto parse_arguments(int argc, char **argv) -> std::tuple<Particle_Set, Computation_Info>
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
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    std::cout << "Compute info:" << "\n";
    std::cout << "\tParticle count: " << numParticles << " particles" << "\n";
	std::cout << "\tParticle mass: " << mass << " [kg]" << "\n";
	std::cout << "\tTimestep: " << dt << " [s]" << "\n";
	std::cout << "\tIterations: " << iterations << " steps" << "\n";
	std::cout << "\tThreads: " << threads << "\n\n";

	/// Form the particle set and the computation info
	Particle_Set my_set{numParticles, mass};
	Computation_Info info{dt, iterations, threads};

	return {my_set, info};

}

int main(int argc, char **argv)
{

	/// Get computation info
	auto [particle_set, computation_info] = parse_arguments(argc, argv);

	/// Start measuring time
	auto start = std::chrono::steady_clock::now();

	/// Start the computation
    n_body_computation_dispatcher(particle_set, computation_info);

    /// Finish measuring time
    auto finish = std::chrono::steady_clock::now();

    double duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "Compute elapsed time: " << "\n";
    std::cout << "\t" << duration / 1e3 << " [ms] (" << duration / 1e6 << " [s]) " << "\n";

    std::cout << "Average time per iteration: " << "\n";
    std::cout << "\t" << duration / computation_info.numIterations << " [us]" << "\n";

}
