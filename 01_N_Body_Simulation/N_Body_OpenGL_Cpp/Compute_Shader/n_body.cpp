#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <regex>

#include <glm/glm.hpp>

#include "CLI/CLI.hpp"

struct Computation_Info
{

	Computation_Info(float _timeStep, unsigned _numIterations, unsigned _workGroupSize) :
					 timeStep(_timeStep), numIterations(_numIterations), workGroupSize(_workGroupSize) { }

	float timeStep;
	unsigned numIterations;

	unsigned workGroupSize;

};

bool check_shader_compile_status(GLuint obj) {
    GLint status;
    glGetShaderiv(obj, GL_COMPILE_STATUS, &status);
    if(status == GL_FALSE) {
        GLint length;
        glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &length);
        std::vector<char> log(length);
        glGetShaderInfoLog(obj, length, &length, &log[0]);
        std::cerr << &log[0];
        return false;
    }
    return true;
}

bool check_program_link_status(GLuint obj) {
    GLint status;
    glGetProgramiv(obj, GL_LINK_STATUS, &status);
    if(status == GL_FALSE) {
        GLint length;
        glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &length);
        std::vector<char> log(length);
        glGetProgramInfoLog(obj, length, &length, &log[0]);
        std::cerr << &log[0];
        return false;
    }
    return true;
}

// settings
const unsigned int SCR_WIDTH = 1280;
const unsigned int SCR_HEIGHT = 720;

struct Particle_Set
{

    /**

    Una particula debe guardar su posicion y velocidad actual
    Se agrega el vector de aceleracion/fuerza para no tener que crear

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

    const unsigned numParticles;
    const float mass;

};

auto parse_arguments(int argc, char **argv) -> std::tuple<Particle_Set, Computation_Info>
{
	/// Parse arguments to form the particle set and the computation info
	CLI::App app{"N-Body simulation"};

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

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "N Body Simulation", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    auto [particle_set, computation_info] = parse_arguments(argc, argv);

	/// Variables parsed from program arguments
	unsigned numParticles = particle_set.numParticles;
	float mass = particle_set.mass;
	float timeStep = computation_info.timeStep;
	unsigned numIterations = computation_info.numIterations;
	unsigned workgroupSize = computation_info.workGroupSize;

	/// generate vbos
    unsigned int positions_vbo, velocities_vbo;

    glGenBuffers(1, &positions_vbo);
    glGenBuffers(1, &velocities_vbo);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, velocities_vbo);
    glBufferStorage(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec3) * numParticles, &particle_set.velocities[0], GL_DYNAMIC_DRAW);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, positions_vbo);
    glBufferStorage(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec3) * numParticles, &particle_set.positions[0], GL_DYNAMIC_DRAW);

	unsigned int ssbos[] = {positions_vbo, velocities_vbo};
    glBindBuffersBase(GL_SHADER_STORAGE_BUFFER, 0, 2, ssbos);


    // we need these to properly pass the strings
    const char *source;
    int length;

    std::string n_body_vel_calculation = R"""(
        #version 460
        layout(local_size_x=32) in;

        layout(location = 0) uniform int numParticles;
        layout(location = 1) uniform float mass;
        layout(location = 2) uniform float dt;

        layout(std430, binding=0) buffer pblock { vec3 positions[]; };
        layout(std430, binding=1) buffer vblock { vec3 velocities[]; };

        float dist2(vec3 A, vec3 B)
        {
            vec3 C = A - B;
            return dot( C, C );
        }

        void main()
		{
            int i = int(gl_GlobalInvocationID);

            if (i >= numParticles)
                return;

            const float G = 6.6743e-11f;

            vec3 cur_position = positions[i];

            vec3 force = vec3(0.0);

            for (uint j = 0; j < numParticles; ++j)
            {
                if (i == j)
                    continue;

                vec3 neighbor_position = positions[j];

                float inv_distance2 = 1.0 / dist2(cur_position, neighbor_position);
                vec3 direction = normalize(neighbor_position - cur_position);

                force += G * mass * mass * inv_distance2 * direction;
            }

            vec3 acceleration = force / mass;

            velocities[i] += acceleration * dt;

        }   )""";

    /// Change number of threads per block (worksize)
    std::string worksize_kernel = std::string("local_size_x=") + std::to_string(workgroupSize);
    n_body_vel_calculation = std::regex_replace(n_body_vel_calculation, std::regex("local_size_x=32"), worksize_kernel);


    // program and shader handles
    unsigned int n_body_vel_program, n_body_vel_shader;

    // create and compiler vertex shader
    n_body_vel_shader = glCreateShader(GL_COMPUTE_SHADER);
    source = n_body_vel_calculation.c_str();
    length = n_body_vel_calculation.size();
    glShaderSource(n_body_vel_shader, 1, &source, &length);
    glCompileShader(n_body_vel_shader);
    if(!check_shader_compile_status(n_body_vel_shader)) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    // create program
    n_body_vel_program = glCreateProgram();

    // attach shaders
    glAttachShader(n_body_vel_program, n_body_vel_shader);

    // link the program and check for errors
    glLinkProgram(n_body_vel_program);
    check_program_link_status(n_body_vel_program);

    std::string n_body_pos_calculation = R"""(
        #version 460
        layout(local_size_x=32) in;

        layout(location = 0) uniform int numParticles;
        layout(location = 1) uniform float dt;

        layout(std430, binding=0) buffer pblock { vec3 positions[]; };
        layout(std430, binding=1) buffer vblock { vec3 velocities[]; };


        void main()
		{
            int i = int(gl_GlobalInvocationID);

            if (i >= numParticles)
                return;

            positions[i] += velocities[i] * dt;

        }   )""";

    /// Change number of threads per block (worksize)
    n_body_pos_calculation = std::regex_replace(n_body_pos_calculation, std::regex("local_size_x=32"), worksize_kernel);

    // program and shader handles
    unsigned int n_body_pos_program, n_body_pos_shader;

    // create and compiler vertex shader
    n_body_pos_shader = glCreateShader(GL_COMPUTE_SHADER);
    source = n_body_pos_calculation.c_str();
    length = n_body_pos_calculation.size();
    glShaderSource(n_body_pos_shader, 1, &source, &length);
    glCompileShader(n_body_pos_shader);
    if(!check_shader_compile_status(n_body_pos_shader)) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    // create program
    n_body_pos_program = glCreateProgram();

    // attach shaders
    glAttachShader(n_body_pos_program, n_body_pos_shader);

    // link the program and check for errors
    glLinkProgram(n_body_pos_program);
    check_program_link_status(n_body_pos_program);


    // Carga de variables uniformes
    glUseProgram(n_body_vel_program);
    glUniform1i(0, numParticles);
    glUniform1f(1, mass);
    glUniform1f(2, timeStep);

    glUseProgram(n_body_pos_program);
    glUniform1i(0, numParticles);
    glUniform1f(1, timeStep);

	/// Start measuring time
	unsigned queryID;
	glGenQueries(1, &queryID);
	glBeginQuery(GL_TIME_ELAPSED, queryID);

    // Despacho de programa
    for (unsigned iter = 0; iter < numIterations; ++iter)
    {
        glUseProgram(n_body_vel_program);
        glDispatchCompute((numParticles + workgroupSize - 1)/workgroupSize, 1, 1);

		/// Sync point
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        glUseProgram(n_body_pos_program);
        glDispatchCompute((numParticles + workgroupSize - 1)/workgroupSize, 1, 1);

        /// Sync point
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        //glBindBuffer(GL_SHADER_STORAGE_BUFFER, positions_vbo);
        //glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(glm::vec3) * numParticulas, &particle_set.positions[0]);
    }

	/// Finish measuring time
	glEndQuery(GL_TIME_ELAPSED);
	long long unsigned elapsedTime;
	glGetQueryObjectui64v(queryID, GL_QUERY_RESULT, &elapsedTime);


	std::cout << "Compute elapsed time: " << "\n";
    std::cout << "\t" << elapsedTime/1e3 << " [us] (" << elapsedTime/1e6 << " [ms]) (" << elapsedTime/1e9 << " [s])" << "\n";

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------

    glfwTerminate();
    return 0;
}
