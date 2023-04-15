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


void vpm::CalcDerivative(
        Particle* Particles,
        double DomainL,
        double ParticleRadius,
        double Viscosity,
        int N, 
        std::tuple<double, double, double>* Out) {

    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;
    double PeriodicDist, PeriodicDistanceX, PeriodicDistanceY;
    double Rho;
    vpm::Particle ThisParticle, OtherParticle;
    std::tuple<double, double, double> Vector1;
    double ResultX, ResultY, Xsum, Ysum, Vorsum;

    

#pragma acc parallel loop gang default(present) private( ThisParticle)
    for (size_t i = 0; i < N; i++) {
        Xsum = 0.0;
        Ysum = 0.0;
        Vorsum = 0.0;   
        ThisParticle = Particles[i]; 
#pragma acc loop vector default(present) private(OtherParticle, Vector1, PeriodicDistanceX, PeriodicDistanceY) reduction(+:Xsum) reduction(+:Ysum) reduction(+:Vorsum)
        for (size_t j = 0; j < N; j++) {
            
            if (i == j) continue;
   
            
            OtherParticle = Particles[j];

            auto PerTuple = ThisParticle.PeriodicDistanceVector(OtherParticle, DomainL);
            PeriodicDistanceX = std::get<0>(PerTuple);
            PeriodicDistanceY = std::get<1>(PerTuple);

 
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            Rho = PeriodicDist / ParticleRadius;


            Vector1 = std::make_tuple(
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
            );

            auto ResTuple = Cross(Vector1, std::make_tuple(0, 0, OtherParticle.Vorticity));
            
            ResultX = std::get<0>(ResTuple);
            ResultY = std::get<1>(ResTuple);          
  

            Xsum += ResultX;
            Ysum += ResultY;
            Vorsum += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (OtherParticle.Vorticity - ThisParticle.Vorticity);
        }
        std::get<0>(Out[i]) = -Xsum;
        std::get<1>(Out[i]) = -Ysum;
        std::get<2>(Out[i]) = Vorsum;
    }


    return ;
    
}

std::tuple<double, double> vpm::CalcVelAtPoint(double x, double y, Particle* Particles,
            double DomainL, double ParticleRadius, int N) {

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
