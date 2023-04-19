#include <cmath>
#include "vpm.h"
#include <iostream>
#include <accelmath.h>

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

            double PerDist[2];
            ThisParticle.PeriodicDistanceVector(OtherParticle, DomainL, PerDist);
            PeriodicDistanceX = PerDist[0];
            PeriodicDistanceY = PerDist[1];

 
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            Rho = PeriodicDist / ParticleRadius;


            double Vector1[3] = {
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
            };

            double Vector2[3] { 
                0,0,OtherParticle.Vorticity
            };

            double ResTuple[3];
            Cross(Vector1, Vector2, ResTuple);
            
            ResultX = ResTuple[0];
            ResultY = ResTuple[1];          
  

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

        double periodicDistance[2];
        Particles[j].PeriodicDistanceVector(x, y, DomainL, periodicDistance);
        double PeriodicDistanceX = -periodicDistance[0];
        double PeriodicDistanceY = -periodicDistance[1];
        double PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
        if (PeriodicDist < 1E-10)
            continue;
        double Rho = PeriodicDist / ParticleRadius;

        double Vector1[3] = {
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
        };

        double Vector2[3] { 
            0,0,Particles[j].Vorticity
        };

        double ResTuple[3];
        Cross(Vector1, Vector2, ResTuple);
        
        double ResultX = ResTuple[0];
        double ResultY = ResTuple[1];

        std::get<0>(Out) -= ResultX;
        std::get<1>(Out) -= ResultY;           
    }

    return Out;
}

void vpm::Regrid(
            Particle* Particles,
            Particle* ParticlesRegrid,
            double DomainL, 
            double ParticleRadius,
            int N) {

    size_t Resolution = (size_t) sqrt(N);
    double ResolutionDouble = sqrt(N);

    #pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < N; i++) {
        double x = static_cast<double>(i % Resolution) * ParticleRadius / 2.0 + ParticleRadius / 4.0;
        double y = static_cast<double>(i - i % Resolution) / ResolutionDouble * (DomainL / ResolutionDouble) + ParticleRadius / 4.0;

        ParticlesRegrid[i].PositionX = x;
        ParticlesRegrid[i].PositionY = y;

        double Vel[2];
        auto res = vpm::CalcVelAtPoint(x, y, Particles, DomainL, ParticleRadius, N);
        Vel[0] = std::get<0>(res);
        Vel[1] = std::get<1>(res);

        ParticlesRegrid[i].VelocityX = Vel[0];
        ParticlesRegrid[i].VelocityY = Vel[1];

        ParticlesRegrid[i].Vorticity = 0.0;

        for (int j = 0; j < N; j++) {
            double dist = sqrt( pow((Particles[j].PositionX - x), 2) + pow((Particles[j].PositionY - y), 2));
            double interpDist = dist*2/ParticleRadius;

            double interpFunc = interpDist < 1 ? (1-interpDist*interpDist)*(2 - interpDist)* 0.5 :
                                (interpDist < 2 ? (1 - interpDist)*(2 - interpDist)*(3 - interpDist)/6.0 :
                                                  0.0);

            ParticlesRegrid[i].Vorticity += interpFunc*Particles[j].Vorticity;
        }
    }

    #pragma acc parallel loop gang vector default(present)
    for (int i = 0; i < N; i++) {
        Particles[i] = ParticlesRegrid[i];
    }

}
