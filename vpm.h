#pragma once

#include <tuple>
#include <utility>
#include <vector>

namespace vpm {
    class Particle {
    public:

        Particle(std::tuple<double, double> Position, std::tuple<double, double> Velocity, double Vorticity) :
            Position{std::move(Position)}, Velocity{std::move(Velocity)}, Vorticity{Vorticity} {};

        //Particle() : Particle({0.0, 0.0}, {0.0, 0.0}, 0) {};
        Particle() : Position(std::make_tuple(0.0, 0.0)), Velocity(std::make_tuple(0.0, 0.0)), Vorticity(0) {};

        std::tuple<double, double> PeriodicDistanceVector(Particle Other, double DomainL);
        std::tuple<double, double> PeriodicDistanceVector(double x, double y, double DomainL);
        std::tuple<double, double> Position;
        std::tuple<double, double> Velocity;
        double Vorticity;
    };

    void CalcDerivative(
            Particle* Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity,
            int N,  
            std::tuple<double, double, double>* Out
    );

    std::tuple<double, double> CalcVelAtPoint(double x, double y, Particle* Particles,
            double DomainL, double ParticleRadius, int N);
}
