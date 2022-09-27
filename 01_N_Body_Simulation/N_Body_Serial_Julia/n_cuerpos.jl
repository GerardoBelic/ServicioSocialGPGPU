using StaticArrays
using LinearAlgebra
using Distributions
using Distances

mutable struct Particulas

    #Constructor que inicializa las posiciones aleatoriamente
    #y las velocidades en cero
    Particulas(numParticulas::Int32, masa::Float32) =
    begin
        
        posiciones = Vector{SVector{3, Float32}}(undef, numParticulas)
        velocidades = Vector{SVector{3, Float32}}(undef, numParticulas)

        for i in 1:numParticulas

            x = rand(Uniform(-1000.0f0, 1000.0f0))
            y = rand(Uniform(-1000.0f0, 1000.0f0))
            z = rand(Uniform(-1000.0f0, 1000.0f0))

            posiciones[i] = SA_F32[x, y, z]

            velocidades[i] = SA_F32[0, 0, 0]

        end

        new(posiciones, velocidades, numParticulas, masa)

    end

    posiciones::Vector{SVector{3, Float32}}
    velocidades::Vector{SVector{3, Float32}}

    numParticulas::Int32
    masa::Float32
    

end

function iteracion!(particulas::Particulas, ΔT::Float32)

    #Constante de gravitación universal
    G::Float32 = 6.6743f-11

    for i in 1:particulas.numParticulas

        fuerza::SVector{3, Float32} = SA_F32[0.0f0, 0.0f0, 0.0f0]

        for j in 1:particulas.numParticulas

            if i==j
                continue

            end

            distancia² = sqeuclidean(particulas.posiciones[i], particulas.posiciones[j])

            direccion = normalize(particulas.posiciones[j] - particulas.posiciones[i])

            fuerza += G * particulas.masa^2 * direccion / distancia²

        end

        aceleracion::SVector{3, Float32} = fuerza / particulas.masa

        particulas.velocidades[i] += aceleracion * ΔT

    end

    for i in 1:particulas.numParticulas

        particulas.posiciones[i] += particulas.velocidades[i] * ΔT

    end

end