#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>

#include <cuda_runtime.h>
#include <helper_math.h>

#include "CLI/CLI.hpp"

#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <helper_cuda.h>       // helper function CUDA error checking and initialization
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

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

    std::vector<float3> positions;
    std::vector<float3> velocities;

    const int numParticles;
    const float mass;

};

__device__ float dist2(float3 A, float3 B)
{
    float3 C = A - B;
    return dot(C, C);
}

__global__ void n_body_calc(float3 *positions, float3 *velocities,
								 unsigned numParticles, float mass, float deltaTime)
{
    cg::grid_group g = cg::this_grid();

    unsigned i = g.thread_rank();

    if (i >= numParticles)
        return;

    float G = 6.6743e-11f;

	float3 cur_position = positions[i];
	float3 force = make_float3(0.0, 0.0, 0.0);

    for (unsigned j = 0; j < numParticles; ++j)
    {
        if (i == j)
            continue;

		float distance2 = dist2(cur_position, positions[j]);
        float3 direction = normalize(positions[j] - cur_position);

        force += G * mass * mass * direction / distance2;
    }

	float3 acceleration = force / mass;

    velocities[i] += acceleration * deltaTime;

    /// Sync
    cg::sync(g);

    positions[i] += velocities[i] * deltaTime;

    /// Sync again
    cg::sync(g);


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

	float3* d_positions = nullptr;
	float3* d_velocities = nullptr;

	unsigned size_memory = numParticles * sizeof(float3);

	/// Start measuring time
	auto start = std::chrono::steady_clock::now();

	cudaMalloc((void **)&d_positions, size_memory);
	cudaMalloc((void **)&d_velocities, size_memory);

	/// Wait for cudaMalloc to finish
	cudaDeviceSynchronize();

	/// Wait for malloc to finish
	auto device_malloc = std::chrono::steady_clock::now();

	cudaMemcpy(d_positions, &particle_set.positions[0], size_memory, cudaMemcpyHostToDevice);
	cudaMemcpy(d_velocities, &particle_set.velocities[0], size_memory, cudaMemcpyHostToDevice);

	/// Wait for memcpy to finish
	cudaDeviceSynchronize();

	auto memcpy_to_device = std::chrono::steady_clock::now();

	///-----------------------------------------------------------------------------------------
	cudaDeviceProp deviceProp;
    int devID = findCudaDevice(argc, (const char **)argv);
    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

    // This sample requires being run on a device that supports Cooperative Kernel Launch
    if (!deviceProp.cooperativeLaunch)
    {
        printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", devID);
        exit(EXIT_WAIVED);
    }

    int numSms = deviceProp.multiProcessorCount;

    int numBlocksPerSm = 0;
    checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, n_body_calc, workgroupSize, 0));

    dim3 dimGrid((numParticles + workgroupSize - 1) / workgroupSize, 1, 1);
    //dim3 dimGrid(numSms*numBlocksPerSm, 1, 1);
    dim3 dimBlock(workgroupSize, 1, 1);

    void* kernelArgs[] =
    {
        (void*)&d_positions,
        (void*)&d_velocities,
        (void*)&numParticles,
        (void*)&mass,
        (void*)&timeStep
    };

    ///-----------------------------------------------------------------------------------------


	for (unsigned iter = 0; iter < numIterations; ++iter)
	{
		//n_body_calc<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, numParticles, mass, timeStep);
        //checkCudaErrors(cudaLaunchCooperativeKernel((void *)n_body_calc, dimGrid, dimBlock, kernelArgs, 0, 0));
        checkCudaErrors(cudaLaunchKernel((void *)n_body_calc, dimGrid, dimBlock, kernelArgs, 0, 0));
		/*{
            cudaError_t cudaerr = cudaDeviceSynchronize();
            if (cudaerr != cudaSuccess)
                printf("kernel launch failed with error \"%s\".\n",
                       cudaGetErrorString(cudaerr));
        }*/

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
