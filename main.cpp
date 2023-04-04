#include <cmath>
#include <iostream>
#include <ranges>
#include "vpm.h"

#define MYSTERY_FRACTION 0.1

int main() {
    // sim params
    double L = 2*M_PI; // meters
    size_t Resolution = 16; // initialize particles as 64*64 square grid
    auto ResolutionDouble = static_cast<double>(Resolution);
    double Viscosity = 1E-2; // kinematic Viscosity m^2/s
    double dt = 1E-1;
    size_t nt = 100;

    size_t N = Resolution * Resolution;
    double ParticleRad = L / ResolutionDouble * 2;
    double ParticleVol = M_PI * ParticleRad * ParticleRad;

    std::vector<vpm::Particle> Particles(N);

    // Initialize Particles
    for (size_t i = 0; i < N; i++) {
        double x = static_cast<double>(i % Resolution) * ParticleRad / 2.0 + ParticleRad / 4.0;
        double y = static_cast<double>(i - i % Resolution) / ResolutionDouble * (L / ResolutionDouble) + ParticleRad / 4.0;
        double u = std::cos(x) * std::sin(y);
        double v = -std::sin(x) * std::cos(y);
        double omega = -2 * std::cos(x) * std::cos(y);
        Particles[i] = {{x, y}, {u, v}, omega * ParticleVol * MYSTERY_FRACTION};
    }

    // plot
#ifdef PLOT
    plotField(Particles, L, 0.5);
#endif

    for (size_t t = 0; t < nt; t++) {
        auto Derivatives = vpm::CalcDerivative(Particles, L, ParticleRad, Viscosity);
        for (size_t i = 0; i < N; i++) {
            vpm::Particle& Particle = Particles[i];
            auto [dX, dY, dOmega] = Derivatives[i];
            std::get<0>(Particle.Position) += dX * dt;
            std::get<1>(Particle.Position) += dY * dt;
            std::get<0>(Particle.Velocity) = dX;
            std::get<1>(Particle.Velocity) = dY;
            Particle.Vorticity += dOmega * dt;

            if (std::get<0>(Particle.Position) < 0) {
                std::get<0>(Particle.Position) += L;
            }
            if (std::get<1>(Particle.Position) < 0) {
                std::get<1>(Particle.Position) += L;
            }
            if (std::get<0>(Particle.Position) > L) {
                std::get<0>(Particle.Position) -= L;
            }
            if (std::get<1>(Particle.Position) > L) {
                std::get<1>(Particle.Position) -= L;
            }
        }

        std::cout << t << std::endl;
#ifdef PLOT
        plotField(Particles, L, 0.5);
#endif
    }
}
