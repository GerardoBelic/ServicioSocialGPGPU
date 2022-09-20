#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define GLM_ENABLE_EXPERIMENTAL
//#define GLM_FORCE_CUDA
#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"

// __device__ float dist2()
// {


    
// }

struct Particulas
{

    /**

    Una particula debe guardar su posicion y velocidad actual
    Se agrega el vector de aceleracion/fuerza para no tener que crear

    */
    Particulas(int _numParticulas, float _masa = 1.0f) : posiciones(_numParticulas), velocidades(_numParticulas), numParticulas(_numParticulas), masa(_masa)
    {
        std::random_device rd;
        std::mt19937 e2(rd());
        std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

        for (auto& posicion : posiciones)
        {
            posicion.x = dist(e2);
            posicion.y = dist(e2);
            posicion.z = dist(e2);
        }
    }

    std::vector<glm::vec3> posiciones;
    std::vector<glm::vec3> velocidades;

    const int numParticulas;
    const float masa;

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

int main()
{

	int numParticulas = 1024 * 2;
	float masa = 1.0e9f;
	float deltaTiempo = 1.0f;

	Particulas set_1(numParticulas, masa);
	
	glm::vec3* d_posiciones = nullptr;
	glm::vec3* d_velocidades = nullptr;
	
	int size_memory = numParticulas * sizeof(glm::vec3);
	
	cudaMalloc((void **)&d_posiciones, size_memory);
	cudaMalloc((void **)&d_velocidades, size_memory);
	
	cudaMemcpy(d_posiciones, &set_1.posiciones[0], size_memory, cudaMemcpyHostToDevice);
	cudaMemcpy(d_velocidades, &set_1.velocidades[0], size_memory, cudaMemcpyHostToDevice);

	dim3 threadsPerBlock(32);
	dim3 blocksPerGrid((set_1.numParticulas + 31) / 32);
	
	std::cout << set_1.posiciones[0].x << " " << set_1.posiciones[0].y << " " << set_1.posiciones[0].z << std::endl;
	
	for (int i = 0; i < 2000; ++i)
	{
		iteracionNCuerpos<<<blocksPerGrid, threadsPerBlock>>>(d_posiciones, d_velocidades, numParticulas, masa, deltaTiempo);
		integracionNCuerpos<<<blocksPerGrid, threadsPerBlock>>>(d_posiciones, d_velocidades, numParticulas, deltaTiempo);
		cudaMemcpy(&set_1.posiciones[0], d_posiciones, size_memory, cudaMemcpyDeviceToHost);
		std::cout << set_1.posiciones[0].x << " " << set_1.posiciones[0].y << " " << set_1.posiciones[0].z << std::endl;
	}
	
	//cudaMemcpy(&set_1.posiciones[0], d_posiciones, size_memory, cudaMemcpyDeviceToHost);
	
	


}