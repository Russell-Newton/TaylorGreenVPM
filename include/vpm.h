#pragma once

#include <tuple>
#include <utility>
#include <vector>
#include <iostream>
#include <accelmath.h>

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

    void Regrid(
            Particle* Particles,
            Particle* ParticlesRegrid,
            double DomainL,
            double ParticleRadius,
            int N
    );
}

#pragma acc routine seq
inline void vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL, double distanceVector[2]) const {
    double thisX = this->PositionX;
    double thisY = this->PositionY;
    double otherX = Other.PositionX;
    double otherY = Other.PositionY;

    distanceVector[0] = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    distanceVector[1] = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
}

#pragma acc routine seq
inline double Kernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * (1 - exp(-Rho * Rho / 2));
}

#pragma acc routine seq
inline double ViscousKernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * exp(-Rho * Rho / 2);
}

#pragma acc routine seq
inline void Cross(double a[3], double b[3], double c[3]) {

    c[0] =         a[1] * b[2] - a[2] * b[1];
    c[1] =     a[2] * b[0] - a[0] * b[2];
    c[2] =     a[0] * b[1] - a[1] * b[2];
}

#pragma acc routine seq
inline void vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL, double distanceVector[2]) const {
    double thisX = this->PositionX;
    double thisY = this->PositionY;

    distanceVector[0] = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    distanceVector[1] = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
}
