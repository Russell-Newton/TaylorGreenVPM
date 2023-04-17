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
        double * ParticleX,
        double * ParticleY,
        double * ParticleVort,
        double * ParticleU,
        double * ParticleV,
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

        thisX = ParticleX[i];
        thisY = ParticleY[i];
        
#pragma acc loop gang vector default(present) private(OtherX, OtherY, PeriodicDistanceX, PeriodicDistanceY, PeriodicDist, Rho) reduction(+:dX) reduction(+:dY) reduction(+:dOmega)
        for (size_t j = 0; j < N; j++) {
            
            if (i == j) continue;
   
            ParticleRadius2 = ParticleRadius * ParticleRadius;
            ParticleVol = ParticleRadius2 * M_PI;

            OtherX = ParticleX[j];
            OtherY = ParticleY[j];

            PeriodicDistanceX = std::fmod((thisX - OtherX + L / 2.0) + L, L) - L / 2;
            PeriodicDistanceY = std::fmod((thisY - OtherY + L / 2.0) + L, L) - L / 2;
 
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            Rho = PeriodicDist / ParticleRadius;

            dX +=  (Kernel_acc(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist)) * ParticleVort[j];
            dY +=  (-Kernel_acc(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist)) * ParticleVort[j];


            dOmega += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel_acc(Rho) * ParticleVol * (ParticleVort[j] - ParticleVort[i]);
        }
        
        //printf("print %d %f %f %f \n", i, dX, dY, dOmega);

        ParticleX[i] += (-dX * dt);
        ParticleY[i] += (-dY * dt);
        ParticleU[i] = -dX;
        ParticleV[i] = -dY;
        ParticleVort[i] += dOmega * dt;

        
        if (ParticleX[i] < 0) {
            ParticleX[i] += L;
        }
        if (ParticleY[i] < 0) {
            ParticleY[i] += L;
        }
        if (ParticleX[i] > L) {
            ParticleX[i] -= L;
        }
        if (ParticleY[i] > L) {
            ParticleY[i] -= L;
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
