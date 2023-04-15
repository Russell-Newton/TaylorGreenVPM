#include <cmath>
#include "vpm.h"
#include <iostream>
#include <accelmath.h>

#pragma acc routine seq
std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL) {

    auto [thisX, thisY] = this->Position;
    auto [otherX, otherY] = Other.Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return std::make_tuple(diffX, diffY);
}

#pragma acc routine seq
std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL) {

    auto [thisX, thisY] = this->Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return std::make_tuple(diffX, diffY);
}

#pragma acc routine seq
double Kernel(double Rho) {

    return (1.0 / (2.0 * M_PI)) * (1 - exp(-Rho * Rho / 2));
}

#pragma acc routine seq
double ViscousKernel(double Rho) {

    return (1.0 / (2.0 * M_PI)) * exp(-Rho * Rho / 2);
}

#pragma acc routine seq
std::tuple<double, double, double> Cross(std::tuple<double, double, double> a, std::tuple<double, double, double> b) {

    auto [ax, ay, az] = a;
    auto [bx, by, bz] = b;

    return std::make_tuple(
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx
    );
}


std::tuple<double, double, double>* vpm::CalcDerivative(
        Particle* Particles,
        double DomainL,
        double ParticleRadius,
        double Viscosity) {

    int N = sizeof(Particles) / sizeof(Particles[0]) ;
    std::tuple<double, double, double> Out[N];


    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;
    double PeriodicDist, PeriodicDistanceX, PeriodicDistanceY;
    double Rho;
    vpm::Particle ThisParticle, OtherParticle;
    std::tuple<double, double, double> Vector1;
    double ResultX, ResultY;


#pragma acc parallel loop gang vector collapse(2) default(present) private( ThisParticle, OtherParticle, PeriodicDist, PeriodicDistanceX, PeriodicDistanceY, Vector1)
    for (size_t i = 0; i < N; i++) {
        
        for (size_t j = 0; j < N; j++) {
            

            if (i == j) continue;
            Out[i] = std::make_tuple(0.0, 0.0, 0.0);
            ThisParticle = Particles[i];
            OtherParticle = Particles[j];

            std::make_tuple(PeriodicDistanceX, PeriodicDistanceY) = ThisParticle.PeriodicDistanceVector(OtherParticle, DomainL);
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            Rho = PeriodicDist / ParticleRadius;

            Vector1 = std::make_tuple(
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
            );
            std::make_tuple(ResultX, ResultY, 0.0) = Cross(Vector1, std::make_tuple(0, 0, OtherParticle.Vorticity));
            std::get<0>(Out[i]) -= ResultX;
            std::get<1>(Out[i]) -= ResultY;
            std::get<2>(Out[i]) += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (OtherParticle.Vorticity - ThisParticle.Vorticity);
        }
    }

    return Out;
}

std::tuple<double, double> vpm::CalcVelAtPoint(double x, double y, std::vector<Particle> Particles,
            double DomainL, double ParticleRadius) {

    size_t N = Particles.size();
    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;

    std::tuple<double, double> Out(0., 0.);


    for (size_t j = 0; j < N; j++) {

        auto [PeriodicDistanceX, PeriodicDistanceY] = Particles[j].PeriodicDistanceVector(x, y, DomainL);
        PeriodicDistanceX = -PeriodicDistanceX;
        PeriodicDistanceY = -PeriodicDistanceY;
        double PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
        if (PeriodicDist < 1E-10)
            continue;
        double Rho = PeriodicDist / ParticleRadius;

        std::tuple<double, double, double> Vector1 = std::make_tuple(
                Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                0
        );
        auto [ResultX, ResultY, _] = Cross(Vector1, std::make_tuple(0, 0, Particles[j].Vorticity));
        std::get<0>(Out) -= ResultX;
        std::get<1>(Out) -= ResultY;           
    }

    return Out;
}
