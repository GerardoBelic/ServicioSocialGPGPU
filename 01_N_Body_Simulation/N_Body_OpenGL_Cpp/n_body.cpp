#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>

#include <glm/glm.hpp>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

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

// helper to check and display for shader linker error
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

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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

    // we need these to properly pass the strings
    const char *source;
    int length;

    std::string n_body_iteration = R"""(
        #version 460
        layout(local_size_x=32) in;

        layout(location = 0) uniform int numParticulas;
        layout(location = 1) uniform float masa;
        layout(location = 2) uniform float dt;

        layout(std430, binding=0) buffer pblock { vec3 posiciones[]; };
        layout(std430, binding=1) buffer vblock { vec3 velocidades[]; };

        float dist2(vec3 A, vec3 B)
        {
            vec3 C = A - B;
            return dot( C, C );
        }

        void main() {
            int i = int(gl_GlobalInvocationID);

            if (i >= numParticulas)
                return;

            float G = 6.6743e-11f;

            vec3 fuerza = vec3(0.0);

            for (int j = 0; j < numParticulas; ++j)
            {
                if (i == j)
                    continue;

                float distancia2 = dist2(posiciones[i], posiciones[j]);
                vec3 direccion = normalize(posiciones[j] - posiciones[i]);

                fuerza += G * masa * masa * direccion / distancia2;
            }

            vec3 aceleracion = fuerza / masa;

            velocidades[i] += aceleracion * dt;

        }   )""";

    // program and shader handles
    unsigned int acceleration_program, acceleration_shader;

    // create and compiler vertex shader
    acceleration_shader = glCreateShader(GL_COMPUTE_SHADER);
    source = n_body_iteration.c_str();
    length = n_body_iteration.size();
    glShaderSource(acceleration_shader, 1, &source, &length);
    glCompileShader(acceleration_shader);
    if(!check_shader_compile_status(acceleration_shader)) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    // create program
    acceleration_program = glCreateProgram();

    // attach shaders
    glAttachShader(acceleration_program, acceleration_shader);

    // link the program and check for errors
    glLinkProgram(acceleration_program);
    check_program_link_status(acceleration_program);

    std::string n_body_integration = R"""(
        #version 460
        layout(local_size_x=32) in;

        layout(location = 0) uniform int numParticulas;
        layout(location = 1) uniform float dt;

        layout(std430, binding=0) buffer pblock { vec3 posiciones[]; };
        layout(std430, binding=1) buffer vblock { vec3 velocidades[]; };


        void main() {
            int i = int(gl_GlobalInvocationID);

            if (i >= numParticulas)
                return;

            posiciones[i] += velocidades[i] * dt;

        }   )""";

    // program and shader handles
    unsigned int integrate_program, integrate_shader;

    // create and compiler vertex shader
    integrate_shader = glCreateShader(GL_COMPUTE_SHADER);
    source = n_body_integration.c_str();
    length = n_body_integration.size();
    glShaderSource(integrate_shader, 1, &source, &length);
    glCompileShader(integrate_shader);
    if(!check_shader_compile_status(integrate_shader)) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    // create program
    integrate_program = glCreateProgram();

    // attach shaders
    glAttachShader(integrate_program, integrate_shader);

    // link the program and check for errors
    glLinkProgram(integrate_program);
    check_program_link_status(integrate_program);

    ///Creacion de las particulas y copia a buffers
    int numParticulas = 1024 * 2;
	float masa = 1.0e9f;
	float deltaTiempo = 1.0f;

	Particulas set_1(numParticulas, masa);

	// generate positions_vbos and vaos
    unsigned int posiciones_vbo, velocidades_vbo;

    glGenBuffers(1, &posiciones_vbo);
    glGenBuffers(1, &velocidades_vbo);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, velocidades_vbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec3) * numParticulas, &set_1.velocidades[0], GL_STATIC_DRAW);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, posiciones_vbo);
    // fill with initial data
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec3) * numParticulas, &set_1.posiciones[0], GL_STATIC_DRAW);


	unsigned int ssbos[] = {posiciones_vbo, velocidades_vbo};
    glBindBuffersBase(GL_SHADER_STORAGE_BUFFER, 0, 2, ssbos);

    // Carga de variables uniformes
    glUseProgram(acceleration_program);
    glUniform1i(0, numParticulas);
    glUniform1f(1, masa);
    glUniform1f(2, deltaTiempo);

    glUseProgram(integrate_program);
    glUniform1i(0, numParticulas);
    glUniform1f(1, deltaTiempo);

    std::cout << set_1.posiciones[0].x << " " << set_1.posiciones[0].y << " " << set_1.posiciones[0].z << "\n";

    // Despacho de programa
    for (int i = 0; i < 2000; ++i)
    {
        glUseProgram(acceleration_program);

        glDispatchCompute((numParticulas + 31)/32, 1, 1);

        glUseProgram(integrate_program);

        glDispatchCompute((numParticulas + 31)/32, 1, 1);

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, posiciones_vbo);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(glm::vec3) * numParticulas, &set_1.posiciones[0]);

        std::cout << set_1.posiciones[0].x << " " << set_1.posiciones[0].y << " " << set_1.posiciones[0].z << "\n";
    }


    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------

    glfwTerminate();
    return 0;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}
