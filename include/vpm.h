#pragma once

#include <tuple>
#include <utility>
#include <vector>
#include <iostream>

namespace vpm {
    struct Particle {
        Particle(double PositionX, double PositionY, double VelocityX, double VelocityY, double Vorticity) :
            PositionX{PositionX}, PositionY{PositionY}, VelocityX{VelocityX}, VelocityY{VelocityY}, Vorticity{Vorticity} {};

        Particle() : Particle(0, 0, 0, 0, 0) {};

        void PeriodicDistanceVector(Particle Other, double DomainL, double distanceVector[2]) const;
        void PeriodicDistanceVector(double x, double y, double DomainL, double distanceVector[2]) const;
        double PositionX;
        double PositionY;
        double VelocityX;
        double VelocityY;
        double Vorticity;

        friend std::ostream& operator<<(std::ostream& os, const Particle& part) {
            return os << part.PositionX << " " << part.PositionY << " " << part.VelocityX << " " << part.VelocityY << " " << part.Vorticity;
        }
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

inline void vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL, double distanceVector[2]) const {
    double thisX = this->PositionX;
    double thisY = this->PositionY;
    double otherX = Other.PositionX;
    double otherY = Other.PositionY;

    distanceVector[0] = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    distanceVector[1] = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
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

inline void vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL, double distanceVector[2]) const {
    double thisX = this->PositionX;
    double thisY = this->PositionY;

    distanceVector[0] = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    distanceVector[1] = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
}
