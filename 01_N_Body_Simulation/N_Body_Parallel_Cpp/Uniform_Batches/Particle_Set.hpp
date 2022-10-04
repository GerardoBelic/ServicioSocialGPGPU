#ifndef PARTICLE_SET
#define PARTICLE_SET

#include <vector>
#include <random>

#include <glm/glm.hpp>

struct Particle_Set
{
    /**
    A particle has a position and velocity
	All particles shares mass in this scenario
	Positions are randomly selected with uniform distribution
    */
    Particle_Set(unsigned _numParticles, float _mass) : positions(_numParticles), velocities(_numParticles), numParticles(_numParticles), mass(_mass)
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

#endif
