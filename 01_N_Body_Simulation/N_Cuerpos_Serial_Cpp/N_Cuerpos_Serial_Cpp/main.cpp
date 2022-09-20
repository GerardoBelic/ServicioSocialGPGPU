#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

struct Particulas
{

    /**

    Una particula debe guardar su posicion y velocidad actual
    Se agrega el vector de aceleracion/fuerza para no tener que crear

    */
    Particulas(int _numParticulas, float _masa = 1.0f) : posiciones(_numParticulas), velocidades(_numParticulas), numParticulas(_numParticulas), masa(_masa)
    {

        for (auto& posicion : posiciones)
        {
            std::random_device rd;
            std::mt19937 e2(rd());
            std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

            posicion.x = dist(e2);
            posicion.y = dist(e2);
        }

    }

    std::vector<glm::vec2> posiciones;
    std::vector<glm::vec2> velocidades;

    const int numParticulas;
    const float masa;

};

void iteracionNCuerpos(Particulas &particulas, const float deltaTiempo)
{

    std::vector<glm::vec2> &posiciones = particulas.posiciones;
    std::vector<glm::vec2> &velocidades = particulas.velocidades;
    const float masa = particulas.masa;
    const float numParticulas = particulas.numParticulas;

    ///Constante de gravitación universal
    const float G = 6.6743e-11f;

    /// Ciclo para calcular la nueva velocidad de cada particula calculando su vector de fuerza y aceleración
    for (int i = 0; i < numParticulas; ++i)
    {

        glm::vec2 fuerza(0.0f, 0.0f);

        /// Ciclo para sumar las fuerzas
        for (int j = 0; j < numParticulas; ++j)
        {

            /// Si i y j son el mismo cuerpo, no se efectua ninguna operacion
            if (i == j)
                continue;

            float distancia2 = glm::distance2(posiciones[i], posiciones[j]);
            glm::vec2 direccion = glm::normalize(posiciones[j] - posiciones[i]);

            fuerza += G * masa * masa * direccion / distancia2;
            
        }
        
        glm::vec2 aceleracion = fuerza / masa;

        velocidades[i] += aceleracion * deltaTiempo;

    }

    /// Ciclo para calcular las nuevas posiciones con las velocidades actualizadas
    for (int i = 0; i < numParticulas; ++i)
    {

        posiciones[i] += velocidades[i] * deltaTiempo;

    }

    std::cout << posiciones[0].x << ' ' << posiciones[0].y << '\n';

}


int main()
{

    Particulas set_1(1024 * 2, 1.0e9f);

    for (int i = 0; i < 20; ++i)
        iteracionNCuerpos(set_1, 1.0f);

}
