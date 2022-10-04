#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"

#include "CLI/CLI.hpp"

struct Computation_Info
{

	Computation_Info(float _timeStep, unsigned _numIterations, unsigned _workGroupSize) :
					 timeStep(_timeStep), numIterations(_numIterations), workGroupSize(_workGroupSize) { }

	float timeStep;
	unsigned numIterations;
	
	unsigned workGroupSize;
	
};

struct Particle_Set
{

    /**
	A particle must store its currect position and velocty.
    */
    Particle_Set(unsigned _numParticles, float _mass = 1.0f) : positions(_numParticles), velocities(_numParticles), numParticles(_numParticles), mass(_mass)
    {
        std::random_device rd;
        std::mt19937 e2(rd());
        std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

        for (auto& pos : positions)
        {
            pos.x = dist(e2);
            pos.y = dist(e2);
            pos.z = dist(e2);
        }
    }

    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> velocities;

    const int numParticles;
    const float mass;

};

__global__ void iteracionNCuerpos(glm::vec3 *posiciones, glm::vec3 *velocidades,
								  int numParticulas, float masa, float deltaTiempo)
{

	glm::vec3 fuerza(0.0f, 0.0f, 0.0f);

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= numParticulas)
        return;

    float G = 6.6743e-11f;

    for (int j = 0; j < numParticulas; ++j)
    {

        if (i == j)
            continue;

		float distancia2 = glm::distance2(posiciones[i], posiciones[j]);
        glm::vec3 direccion = glm::normalize(posiciones[j] - posiciones[i]);

        fuerza += G * masa * masa * direccion / distancia2;
        

    }
	
	glm::vec3 aceleracion = fuerza / masa;

    velocidades[i] += aceleracion * deltaTiempo;
	
}

__global__ void integracionNCuerpos(glm::vec3 *posiciones, glm::vec3 *velocidades,
								  int numParticulas, float deltaTiempo)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i >= numParticulas)
        return;
		
	posiciones[i] += velocidades[i] * deltaTiempo;
}

auto parse_arguments(int argc, char **argv) -> std::tuple<Particle_Set, Computation_Info>
{
	/// Parse arguments to form the particle set and the computation info
	CLI::App app{"N-Body serial simulation"};

	unsigned numParticles = 0;
	float mass = 1.0f;
	float dt = 1.0f;
	unsigned iterations = 0;
	unsigned workgroupSize = 32;

    app.add_option("--particles", numParticles, "Particle count") -> required();
	app.add_option("--mass", mass, "Mass of each particle (def=1.0[kg])");
	app.add_option("--dt", dt, "Timestep between each iteration (def=1.0[s])");
	app.add_option("--iterations", iterations, "Number of iterations of the simulation") -> required();
	app.add_option("--workgroup_size", workgroupSize, "Number of threads per workgroup (def=32)");

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
	std::cout << "\tWorkgroup size: " << workgroupSize << " threads" << "\n\n";

	/// Form the particle set and the computation info
	Particle_Set my_set{numParticles, mass};
	Computation_Info info{dt, iterations, workgroupSize};

	return {my_set, info};

}

int main(int argc, char **argv)
{

	auto [particle_set, computation_info] = parse_arguments(argc, argv);
	
	/// Variables parsed from program arguments
	unsigned numParticles = particle_set.numParticles;
	float mass = particle_set.mass;
	float timeStep = computation_info.timeStep;
	unsigned numIterations = computation_info.numIterations;
	unsigned workgroupSize = computation_info.workGroupSize;
	
	/// Kernel variables for dispatching
	dim3 threadsPerBlock = workgroupSize;
	dim3 blocksPerGrid = ((numParticles + workgroupSize - 1) / workgroupSize);
	
	glm::vec3* d_positions = nullptr;
	glm::vec3* d_velocities = nullptr;
	
	unsigned size_memory = numParticles * sizeof(glm::vec3);
	
	/// Start measuring time
	auto start = std::chrono::steady_clock::now();
	
	cudaMalloc((void **)&d_positions, size_memory);
	cudaMalloc((void **)&d_velocities, size_memory);
	
	cudaMemcpy(d_positions, &particle_set.positions[0], size_memory, cudaMemcpyHostToDevice);
	cudaMemcpy(d_velocities, &particle_set.velocities[0], size_memory, cudaMemcpyHostToDevice);
	
	auto memcpy_to_device = std::chrono::steady_clock::now();
	
	for (unsigned iter = 0; iter < numIterations; ++iter)
	{
		iteracionNCuerpos<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, numParticles, mass, timeStep);
		integracionNCuerpos<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, numParticles, timeStep);
	}
	
	auto kernel_compute = std::chrono::steady_clock::now();
	
	cudaMemcpy(&particle_set.positions[0], d_positions, size_memory, cudaMemcpyDeviceToHost);
	cudaMemcpy(&particle_set.velocities[0], d_velocities, size_memory, cudaMemcpyDeviceToHost);
	
	/// Finish measuring time
    auto finish = std::chrono::steady_clock::now();
	
	double duration_host_device = std::chrono::duration_cast<std::chrono::microseconds>(memcpy_to_device - start).count();
	double duration_compute = std::chrono::duration_cast<std::chrono::microseconds>(kernel_compute - memcpy_to_device).count();
	double duration_device_host = std::chrono::duration_cast<std::chrono::microseconds>(finish - kernel_compute).count();

    std::cout << "Time to copy memory from host to device: " << "\n";
    std::cout << "\t" << duration_host_device << " [us]" << "\n";

    std::cout << "Compute elapsed time: " << "\n";
    std::cout << "\t" << duration_compute << " [us]" << "\n";
	
	std::cout << "Time to copy memory from device back to host: " << "\n";
    std::cout << "\t" << duration_device_host << " [us]" << "\n";
	
	std::cout << "Total time: " << "\n";
    std::cout << "\t" << duration_host_device + duration_compute + duration_device_host << " [us]" << "\n";

}