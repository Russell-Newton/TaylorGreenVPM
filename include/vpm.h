#pragma once

#include <tuple>
#include <utility>
#include <vector>

namespace vpm {
    class Particle {
    public:

        Particle(std::tuple<double, double> Position, std::tuple<double, double> Velocity, double Vorticity) :
            Position{std::move(Position)}, Velocity{std::move(Velocity)}, Vorticity{Vorticity} {};

        Particle() : Particle({0, 0}, {0, 0}, 0) {};

        std::tuple<double, double> PeriodicDistanceVector(Particle Other, double DomainL);
        std::tuple<double, double> PeriodicDistanceVector(double x, double y, double DomainL);
        std::tuple<double, double> Position;
        std::tuple<double, double> Velocity;
        double Vorticity;
    };

    std::vector<std::tuple<double, double, double>> CalcDerivative(
            std::vector<Particle> Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity
    );

    std::tuple<double, double> CalcVelAtPoint(double x, double y, std::vector<Particle> Particles,
            double DomainL, double ParticleRadius);
}
