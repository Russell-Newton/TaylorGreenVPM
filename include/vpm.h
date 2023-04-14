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

        std::tuple<double, double> PeriodicDistanceVector(Particle Other, double DomainL) const;
        std::tuple<double, double> PeriodicDistanceVector(double x, double y, double DomainL) const;
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

inline std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL) const {
    auto [thisX, thisY] = this->Position;
    auto [otherX, otherY] = Other.Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0), DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0), DomainL) - DomainL / 2;

    return {diffX, diffY};
}


inline double Kernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * (1 - exp(-Rho * Rho / 2));
}


inline double ViscousKernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * exp(-Rho * Rho / 2);
}


inline std::tuple<double, double, double> Cross(std::tuple<double, double, double> a, std::tuple<double, double, double> b) {
    auto [ax, ay, az] = a;
    auto [bx, by, bz] = b;

    return {
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx
    };
}

inline std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL) const {
    auto [thisX, thisY] = this->Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return {diffX, diffY};
}
