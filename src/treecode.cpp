#include <cmath>
#include <stack>
#include <cstdio>
#include "treecode.h"
#include <iostream>

std::vector<std::tuple<double, double, double>> vpm::CalcDerivativeTreeCode(
            std::vector<Particle> Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity,
            double OpeningAngle
    ) {
    size_t N = Particles.size();
    std::vector<std::tuple<double, double, double>> Out(N, std::tuple<double, double, double>(0, 0, 0));
    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;

    std::list<vpm::QuadTreeCodeNode> tree = std::list<vpm::QuadTreeCodeNode>();
    vpm::QuadTreeCodeNode root = vpm::QuadTreeCodeNode();
    tree.push_back(root);
    // root.build(Particles, DomainL, DomainL / 2, DomainL / 2, tree);
    tree.front().build(Particles, DomainL, DomainL / 2, DomainL / 2, tree);

    std::vector<vpm::QuadTreeCodeNode> treeAsVec = std::vector<vpm::QuadTreeCodeNode>(tree.begin(), tree.end());

    for (size_t i = 0; i < N; i++) {
        vpm::Particle thisParticle = Particles[i];
        // get contribution from tree
        std::stack<size_t> toCheck = std::stack<size_t>();
        toCheck.push(0);

        while (!toCheck.empty()) {
            QuadTreeCodeNode checking = treeAsVec[toCheck.top()];
            toCheck.pop();
            // std::cout << i << " " << checking << std::endl;

            // ignore empty nodes
            if (checking.isEmpty) continue;


            // explore further if needed
            double periodicDistance[2];
            thisParticle.PeriodicDistanceVector(checking.particle, DomainL, periodicDistance);
            double periodicDist = sqrt(periodicDistance[0] * periodicDistance[0] + periodicDistance[1] * periodicDistance[1]);
            if (checking.isLeaf && periodicDist <= 1e-10) continue;
            if (!checking.isLeaf && checking.size / periodicDist >= OpeningAngle) {
                toCheck.push(checking.q1ChildIdx);
                toCheck.push(checking.q2ChildIdx);
                toCheck.push(checking.q3ChildIdx);
                toCheck.push(checking.q4ChildIdx);
                continue;
            }

            double Rho = periodicDist / ParticleRadius;

            double Vector1[3] = {
                    Kernel(Rho) * periodicDistance[0] / (periodicDist * periodicDist),
                    Kernel(Rho) * periodicDistance[1] / (periodicDist * periodicDist),
                    0
            };
            double Vector2[3] = {0, 0, checking.particle.Vorticity};

            double ResTuple[3];
            Cross(Vector1, Vector2, ResTuple);

            double ResultX = ResTuple[0];
            double ResultY = ResTuple[1];          

            std::get<0>(Out[i]) -= ResultX;
            std::get<1>(Out[i]) -= ResultY;
            std::get<2>(Out[i]) += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (checking.particle.Vorticity - thisParticle.Vorticity);
        }
    }

    return Out;
}

void vpm::QuadTreeCodeNode::build(const std::vector<Particle>& particles, double _size, double _centerX, double _centerY, std::list<QuadTreeCodeNode>& tree) {
    size = _size;
    if (particles.empty()) {
        isLeaf = true;
        return;
    }

    isEmpty = false;
    if (particles.size() == 1) {
        particle = particles[0];
        isLeaf = true;
        return;
    }

    double centroidX = 0;
    double centroidY = 0;
    double centerVorticity = 0;
    for (auto _particle : particles) {
        centroidX += _particle.PositionX * _particle.Vorticity;
        centroidY += _particle.PositionY * _particle.Vorticity;
        centerVorticity += _particle.Vorticity;
    }

    if (centerVorticity == 0) {
        centroidX = _centerX;
        centroidY = _centerY;
    } else {
        centroidX /= centerVorticity;
        centroidY /= centerVorticity;
    }

    particle = Particle(centroidX, centroidY, 0, 0, centerVorticity);

    std::vector<Particle> particlesQ1 = std::vector<Particle>();
    std::vector<Particle> particlesQ2 = std::vector<Particle>();
    std::vector<Particle> particlesQ3 = std::vector<Particle>();
    std::vector<Particle> particlesQ4 = std::vector<Particle>();
    for (auto _particle : particles) {
        if (_particle.PositionX >= _centerX && _particle.PositionY >= _centerY) particlesQ1.push_back(_particle);
        else if (_particle.PositionX < _centerX && _particle.PositionY >= _centerY) particlesQ2.push_back(_particle);
        else if (_particle.PositionX < _centerX && _particle.PositionY < _centerY) particlesQ3.push_back(_particle);
        else if (_particle.PositionX >= _centerX && _particle.PositionY < _centerY) particlesQ4.push_back(_particle);
    }

    double subSize = _size / 2;
    double subCenterOffset = subSize / 2;

    q1ChildIdx = tree.size();
    tree.push_back({});
    tree.back().build(particlesQ1, subSize, _centerX + subCenterOffset, _centerY + subCenterOffset, tree);
    q2ChildIdx = tree.size();
    tree.push_back({});
    tree.back().build(particlesQ2, subSize, _centerX - subCenterOffset, _centerY + subCenterOffset, tree);
    q3ChildIdx = tree.size();
    tree.push_back({});
    tree.back().build(particlesQ3, subSize, _centerX - subCenterOffset, _centerY - subCenterOffset, tree);
    q4ChildIdx = tree.size();
    tree.push_back({});
    tree.back().build(particlesQ4, subSize, _centerX + subCenterOffset, _centerY - subCenterOffset, tree);
}

vpm::QuadTreeCodeNode::~QuadTreeCodeNode() {
}

void vpm::CalcDerivativeGPUTree(
        BoxNodeGPU * boxNodes,
        Particle* Particles,
        double DomainL,
        double ParticleRadius,
        double Viscosity,
        int N, 
        int numNodes,
        double OpeningAngle,
        std::tuple<double, double, double>* Out) {

    double ParticleRadius2 = ParticleRadius * ParticleRadius;
    double ParticleVol = ParticleRadius2 * M_PI;
    double PeriodicDist, PeriodicDistanceX, PeriodicDistanceY;
    double Rho;
    vpm::Particle ThisParticle, OtherParticle;
    std::tuple<double, double, double> Vector1;
    double ResultX, ResultY, Xsum, Ysum, Vorsum;

    BoxNodeGPU node;
    double nodeSize = DomainL / sqrt(numNodes);

    #pragma acc parallel loop gang vector default(present) private( ThisParticle)
    for (int i = 0; i < numNodes; i++) {
        boxNodes[i].numParticlesInThisBox = 0;
        boxNodes[i].particle.PositionX = 0.0;
        boxNodes[i].particle.PositionY = 0.0;
        boxNodes[i].particle.Vorticity = 0.0;
        // printf("%p", boxNodes[i].particlesInThisBox);
    }
    
    #pragma acc parallel loop gang vector default(present) private( ThisParticle)
    for (int i = 0; i < numNodes; i++) {
        double boxX = (i % (int)sqrt(numNodes)) * DomainL / sqrt(numNodes);
        double boxY = ((int)(i / sqrt(numNodes))) * DomainL / sqrt(numNodes);
        // printf("i: %d, box x: %f, box y: %f\n", i, boxX, boxY);
        

        double centroidX = 0.;
        double centroidY =0.;
        double centroidVorticity =0.;

        for (int j = 0; j < N; j++) {
            ThisParticle = Particles[j];
            // if (i == 0 && j == 3)
            //     printf("Particle 3 pos %f %f\n", ThisParticle.PositionX, ThisParticle.PositionY);
            if (ThisParticle.PositionX >= boxX && ThisParticle.PositionX < boxX+nodeSize &&
                ThisParticle.PositionY >= boxY && ThisParticle.PositionY < boxY+nodeSize) {

                boxNodes[i].particlesInThisBox[boxNodes[i].numParticlesInThisBox] = j;
                boxNodes[i].numParticlesInThisBox++;

                if (boxNodes[i].numParticlesInThisBox >= 395) {
                    printf("Box %d is full!\n", i);
                }

                centroidX += ThisParticle.PositionX * ThisParticle.Vorticity;
                centroidY += ThisParticle.PositionY * ThisParticle.Vorticity;
                centroidVorticity += ThisParticle.Vorticity;
                
                
            }            
        }
        
        boxNodes[i].particle.PositionX = centroidX / centroidVorticity;
        boxNodes[i].particle.PositionY = centroidY / centroidVorticity;
        boxNodes[i].particle.Vorticity = centroidVorticity;

        // if (i == 0) {
        //     printf("%f %f %f\n", boxNodes[i].particle.PositionX,boxNodes[i].particle.PositionY,boxNodes[i].particle.Vorticity);
        // }
    }
    
#pragma acc parallel loop gang default(present) private( ThisParticle, node,OtherParticle, Vector1, PeriodicDistanceX, PeriodicDistanceY, PeriodicDist)
    for (size_t i = 0; i < N; i++) {
        Xsum = 0.0;
        Ysum = 0.0;
        Vorsum = 0.0;   
        ThisParticle = Particles[i]; 

        #pragma acc loop vector default(present) private(node,OtherParticle, Vector1, PeriodicDistanceX, PeriodicDistanceY, PeriodicDist) reduction(+:Xsum) reduction(+:Ysum) reduction(+:Vorsum)
        for (size_t j = 0; j < numNodes; j++) {
            
            node = boxNodes[j];

            if (node.numParticlesInThisBox == 0) continue;

            double PerDist[2];
            ThisParticle.PeriodicDistanceVector(node.particle, DomainL, PerDist);
            PeriodicDistanceX = PerDist[0];
            PeriodicDistanceY = PerDist[1];
 
            PeriodicDist = sqrt(PeriodicDistanceX * PeriodicDistanceX + PeriodicDistanceY * PeriodicDistanceY);
            

            if (nodeSize/PeriodicDist >= OpeningAngle) {
                for (int k = 0; k < node.numParticlesInThisBox; k++) {
                    if (i == node.particlesInThisBox[k]) continue;
                    OtherParticle = Particles[node.particlesInThisBox[k]];
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
                continue;
            }

            // printf("Got here!\n");
            Rho = PeriodicDist / ParticleRadius;

            double Vector1[3] = {
                    Kernel(Rho) * PeriodicDistanceX / (PeriodicDist * PeriodicDist),
                    Kernel(Rho) * PeriodicDistanceY / (PeriodicDist * PeriodicDist),
                    0
            };

            double Vector2[3] { 
                0,0,node.particle.Vorticity
            };

            double ResTuple[3];
            Cross(Vector1, Vector2, ResTuple);
            
            ResultX = ResTuple[0];
            ResultY = ResTuple[1];          
  

            Xsum += ResultX;
            Ysum += ResultY;
            Vorsum += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (node.particle.Vorticity - ThisParticle.Vorticity);
        }
        std::get<0>(Out[i]) = -Xsum;
        std::get<1>(Out[i]) = -Ysum;
        std::get<2>(Out[i]) = Vorsum;
    }


    return ;
    
}