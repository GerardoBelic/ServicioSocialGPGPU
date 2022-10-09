#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>

#include <cstdio>

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
    Particle_Set(unsigned _numParticles, float _mass = 1.0f) : numParticles(_numParticles), mass(_mass)
    {
		/// Allocate in pinned memory
		unsigned size_memory = numParticles * sizeof(glm::vec3);

		cudaMallocHost((void **)&positions, size_memory);
		cudaMallocHost((void **)&velocities, size_memory);
		
        std::random_device rd;
        std::mt19937 e2(rd());
        std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

        for (unsigned i = 0; i < numParticles; ++i)
        {
            (positions + i)->x = dist(e2);
            (positions + i)->y = dist(e2);
            (positions + i)->z = dist(e2);
			
			(velocities + i)->x = 0.0f;
            (velocities + i)->y = 0.0f;
            (velocities + i)->z = 0.0f;
        }
    }
	
	~Particle_Set()
	{
		//cudaFreeHost(positions);
		//cudaFreeHost(velocities);
	}

    glm::vec3* positions = nullptr;
    glm::vec3* velocities = nullptr;

    const int numParticles;
    const float mass;

};

__global__ void n_body_vel_calc(glm::vec3 *positions, glm::vec3 *velocities,
								  unsigned numParticles, float mass, float deltaTime)
{
	
	/// Shared memory between a thread group
	/// The size is equal to the number of threads in a group;
	/// this implies that size==blockDim
	__shared__ glm::vec3 temp_tile[256];
	

    unsigned i = blockDim.x * blockIdx.x + threadIdx.x;

    //if (i >= numParticles)
        //return;

    float G = 6.6743e-11f;
	
	glm::vec3 cur_position = positions[i];
	glm::vec3 force(0.0f, 0.0f, 0.0f);

    for (unsigned tile = 0; tile < numParticles; tile += 256)
    {
		temp_tile[threadIdx.x] = positions[tile + threadIdx.x];
		__syncthreads();

		for (unsigned j = 0; j < 256; ++j)
		{
			//if (i == j || ((tile == (numParticles - 1) / blockDim.x) && numParticles % j < ))
			//if (i == j || (tile + j) >= numParticles)
				//continue;
			if (i == j)
				continue;

			float distance2 = glm::distance2(cur_position, temp_tile[j]);
			glm::vec3 direction = glm::normalize(temp_tile[j] - cur_position);

			force += G * mass * mass * direction / distance2;
		}
		
		//__syncthreads();
        
    }
	
	glm::vec3 acceleration = force / mass;

    velocities[i] += acceleration * deltaTime;
	
}

__global__ void n_body_pos_calc(glm::vec3 *positions, glm::vec3 *velocities,
								  unsigned numParticles, float deltaTime)
{
	unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i >= numParticles)
        return;
		
	positions[i] += velocities[i] * deltaTime;
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
	//unsigned workgroupSize = computation_info.workGroupSize;
	unsigned workgroupSize = 256;
	
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
	
	/// Wait for cudaMalloc to finish
	cudaDeviceSynchronize();
	
	/// Wait for malloc to finish
	auto device_malloc = std::chrono::steady_clock::now();
	
	cudaMemcpy(d_positions, particle_set.positions, size_memory, cudaMemcpyHostToDevice);
	cudaMemcpy(d_velocities, particle_set.velocities, size_memory, cudaMemcpyHostToDevice);
	
	/// Wait for memcpy to finish
	cudaDeviceSynchronize();
	
	auto memcpy_to_device = std::chrono::steady_clock::now();
	
	for (unsigned iter = 0; iter < numIterations; ++iter)
	{
		//n_body_vel_calc<<<blocksPerGrid, threadsPerBlock, workgroupSize * sizeof(glm::vec3)>>>(d_positions, d_velocities, numParticles, mass, timeStep);
		//n_body_pos_calc<<<blocksPerGrid, threadsPerBlock, workgroupSize * sizeof(glm::vec3)>>>(d_positions, d_velocities, numParticles, timeStep);
		n_body_vel_calc<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, numParticles, mass, timeStep);
		n_body_pos_calc<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, numParticles, timeStep);
	}
	
	/// Wait for computation to finish
	cudaDeviceSynchronize();
	
	auto kernel_compute = std::chrono::steady_clock::now();
	
	cudaMemcpy(&particle_set.positions[0], d_positions, size_memory, cudaMemcpyDeviceToHost);
	cudaMemcpy(&particle_set.velocities[0], d_velocities, size_memory, cudaMemcpyDeviceToHost);
	
	/// Wait for memcpy to finish
	cudaDeviceSynchronize();
	
	/// Finish measuring time
    auto finish = std::chrono::steady_clock::now();
	
	double duration_cuda_malloc = std::chrono::duration_cast<std::chrono::microseconds>(device_malloc - start).count();
	double duration_host_device = std::chrono::duration_cast<std::chrono::microseconds>(memcpy_to_device - device_malloc).count();
	double duration_compute = std::chrono::duration_cast<std::chrono::microseconds>(kernel_compute - memcpy_to_device).count();
	double duration_device_host = std::chrono::duration_cast<std::chrono::microseconds>(finish - kernel_compute).count();

	std::cout << "cudaMalloc elapsed time: " << "\n";
    std::cout << "\t" << duration_cuda_malloc << " [us]" << "\n";

    std::cout << "Time to copy memory from host to device: " << "\n";
    std::cout << "\t" << duration_host_device << " [us]" << "\n";

    std::cout << "Compute elapsed time: " << "\n";
    std::cout << "\t" << duration_compute << " [us] (" << duration_compute/1e3 << " [ms]) (" << duration_compute/1e6 << " [s])" << "\n";
	
	std::cout << "Time to copy memory from device back to host: " << "\n";
    std::cout << "\t" << duration_device_host << " [us]" << "\n";
	
	double total_time = duration_cuda_malloc + duration_host_device + duration_compute + duration_device_host;
	
	std::cout << "Total time: " << "\n";
    std::cout << "\t" << total_time << " [us] (" << total_time/1e3 << " [ms]) (" << total_time/1e6 << " [s])" << "\n";

}