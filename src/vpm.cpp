#include <cmath>
#include "vpm.h"
#include <iostream>


std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL) {
    auto [thisX, thisY] = this->Position;
    auto [otherX, otherY] = Other.Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return {diffX, diffY};
}

std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL) {
    auto [thisX, thisY] = this->Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return {diffX, diffY};
}


double Kernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * (1 - exp(-Rho * Rho / 2));
}


double ViscousKernel(double Rho) {
    return (1.0 / (2.0 * M_PI)) * exp(-Rho * Rho / 2);
}


std::tuple<double, double, double> Cross(std::tuple<double, double, double> a, std::tuple<double, double, double> b) {
    auto [ax, ay, az] = a;
    auto [bx, by, bz] = b;

    return {
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx
    };
}


std::vector<std::tuple<double, double, double>> vpm::CalcDerivative(
        std::vector<Particle> Particles,
        double DomainL,
        double ParticleRadius,
        double Viscosity) {
    size_t N = Particles.size();
    std::vector<std::tuple<double, double, double>> Out(N, std::tuple<double, double, double>(0, 0, 0));
    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;

    for (size_t i = 0; i < N; i++) {
        vpm::Particle ThisParticle = Particles[i];
        for (size_t j = 0; j < N; j++) {
            if (i == j) continue;
            vpm::Particle OtherParticle = Particles[j];

            auto [PeriodicDistanceX, PeriodicDistanceY] = ThisParticle.PeriodicDistanceVector(OtherParticle, DomainL);
            double PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            double Rho = PeriodicDist / ParticleRadius;

            std::tuple<double, double, double> Vector1 = {
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
            };
            auto [ResultX, ResultY, _] = Cross(Vector1, {0, 0, OtherParticle.Vorticity});
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

        std::tuple<double, double, double> Vector1 = {
                Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                0
        };
        auto [ResultX, ResultY, _] = Cross(Vector1, {0, 0, Particles[j].Vorticity});
        std::get<0>(Out) -= ResultX;
        std::get<1>(Out) -= ResultY;           
    }

    return Out;
}
