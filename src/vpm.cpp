#include <cmath>
#include "vpm.h"
#include <iostream>
#include <ctime>


std::vector<std::tuple<double, double, double>> vpm::CalcDerivative(
        std::vector<Particle> Particles,
        double DomainL,
        double ParticleRadius,
        double Viscosity,
        bool printTime) {
    size_t N = Particles.size();
    std::vector<std::tuple<double, double, double>> Out(N, std::tuple<double, double, double>(0, 0, 0));
    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;

    clock_t start = clock();

    for (size_t i = 0; i < N; i++) {
        vpm::Particle ThisParticle = Particles[i];
        for (size_t j = 0; j < N; j++) {
            if (i == j) continue;
            vpm::Particle OtherParticle = Particles[j];

            double periodicDistance[2];
            ThisParticle.PeriodicDistanceVector(OtherParticle, DomainL, periodicDistance);
            double periodicDist = sqrt(periodicDistance[0] * periodicDistance[0] + periodicDistance[1] * periodicDistance[1]);

            double Rho = periodicDist / ParticleRadius;

            std::tuple<double, double, double> Vector1 = {
                    Kernel(Rho) * periodicDistance[0] / (periodicDist * periodicDist),
                    Kernel(Rho) * periodicDistance[1] / (periodicDist * periodicDist),
                    0
            };
            auto [ResultX, ResultY, _] = Cross(Vector1, {0, 0, OtherParticle.Vorticity});
            std::get<0>(Out[i]) -= ResultX;
            std::get<1>(Out[i]) -= ResultY;
            std::get<2>(Out[i]) += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (OtherParticle.Vorticity - ThisParticle.Vorticity);
        }
    }

    clock_t end = clock();

    if (printTime) std::cout << "Time this step: " << static_cast<double>(end - start) / CLOCKS_PER_SEC << "s" << std::endl;

    return Out;
}

std::tuple<double, double> vpm::CalcVelAtPoint(double x, double y, std::vector<Particle> Particles,
            double DomainL, double ParticleRadius) {

    size_t N = Particles.size();
    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;

    std::tuple<double, double> Out(0., 0.);


    for (size_t j = 0; j < N; j++) {

        double periodicDistance[2];
        Particles[j].PeriodicDistanceVector(x, y, DomainL, periodicDistance);
        double PeriodicDistanceX = -periodicDistance[0];
        double PeriodicDistanceY = -periodicDistance[1];
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
