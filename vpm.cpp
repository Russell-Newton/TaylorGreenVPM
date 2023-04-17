#include <cmath>
#include "vpm.h"
#include <iostream>
#include <accelmath.h>


std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(vpm::Particle Other, double DomainL) {

    auto [thisX, thisY] = this->Position;
    auto [otherX, otherY] = Other.Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return std::make_tuple(diffX, diffY);
}


std::tuple<double, double> vpm::Particle::PeriodicDistanceVector(double otherX, double otherY, double DomainL) {

    auto [thisX, thisY] = this->Position;

    double diffX = std::fmod((thisX - otherX + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;
    double diffY = std::fmod((thisY - otherY + DomainL / 2.0) + DomainL, DomainL) - DomainL / 2;

    return std::make_tuple(diffX, diffY);
}

#pragma acc routine seq
double Kernel_acc(double Rho) {

    return (1.0 / (2.0 * M_PI)) * (1 - exp(-Rho * Rho / 2));
}

#pragma acc routine seq
double ViscousKernel_acc(double Rho) {

    return (1.0 / (2.0 * M_PI)) * exp(-Rho * Rho / 2);
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

    return std::make_tuple(
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx
    );
}


void vpm::CalcDerivative(
        Particle* Particles,
        double L,
        double ParticleRadius,
        double Viscosity,
        int N, 
        double dt) {

 
    double ParticleRadius2, ParticleVol, Rho, PeriodicDist, PeriodicDistanceX, PeriodicDistanceY;
    double thisX, thisY, OtherX, OtherY ;
    double  dX, dY, dOmega;


#pragma acc parallel loop gang default(present) private(thisX, thisY, dX, dY, dOmega)
    for (size_t i = 0; i < N; i++) {
        dX = 0.0;
        dY= 0.0;
        dOmega = 0.0;   

        thisX = std::get<0>(Particles[i].Position);
        thisY = std::get<1>(Particles[i].Position);
        
#pragma acc loop gang vector default(present) private(OtherX, OtherY, PeriodicDistanceX, PeriodicDistanceY, PeriodicDist, Rho) reduction(+:dX) reduction(+:dY) reduction(+:dOmega)
        for (size_t j = 0; j < N; j++) {
            
            if (i == j) continue;
   
            ParticleRadius2 = ParticleRadius * ParticleRadius;
            ParticleVol = ParticleRadius2 * M_PI;

            OtherX = std::get<0>(Particles[j].Position);
            OtherY = std::get<1>(Particles[j].Position);

            PeriodicDistanceX = std::fmod((thisX - OtherX + L / 2.0) + L, L) - L / 2;
            PeriodicDistanceY = std::fmod((thisY - OtherY + L / 2.0) + L, L) - L / 2;
 
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            Rho = PeriodicDist / ParticleRadius;

            dX +=  (Kernel_acc(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist)) * Particles[j].Vorticity;
            dY +=  (-Kernel_acc(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist)) * Particles[j].Vorticity;


            dOmega += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel_acc(Rho) * ParticleVol * (Particles[j].Vorticity - Particles[i].Vorticity);
        }
        
        //printf("print %d %f %f %f \n", i, dX, dY, dOmega);

        std::get<0>(Particles[i].Position) += (-dX * dt);
        std::get<1>(Particles[i].Position) += (-dY * dt);
        std::get<0>(Particles[i].Velocity) = -dX;
        std::get<1>(Particles[i].Velocity) = -dY;
        Particles[i].Vorticity += dOmega * dt;

        
        if (std::get<0>(Particles[i].Position) < 0) {
            std::get<0>(Particles[i].Position) += L;
        }
        if (std::get<1>(Particles[i].Position) < 0) {
            std::get<1>(Particles[i].Position) += L;
        }
        if (std::get<0>(Particles[i].Position) > L) {
            std::get<0>(Particles[i].Position) -= L;
        }
        if (std::get<1>(Particles[i].Position) > L) {
            std::get<1>(Particles[i].Position) -= L;
        }
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
