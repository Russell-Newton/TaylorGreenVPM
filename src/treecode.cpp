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

    vpm::QuadTreeCodeNode tree = vpm::QuadTreeCodeNode();
    tree.build(Particles, DomainL, DomainL / 2, DomainL / 2);

    for (size_t i = 0; i < N; i++) {
        vpm::Particle thisParticle = Particles[i];
        // get contribution from tree
        std::stack<QuadTreeCodeNode*> toCheck = std::stack<QuadTreeCodeNode*>();
        toCheck.push(&tree);

        while (!toCheck.empty()) {
            QuadTreeCodeNode* checking = toCheck.top();
            toCheck.pop();

            // ignore empty nodes
            if (checking->isEmpty) continue;


            // explore further if needed
            double periodicDistance[2];
            thisParticle.PeriodicDistanceVector(checking->particle, DomainL, periodicDistance);
            double periodicDist = sqrt(periodicDistance[0] * periodicDistance[0] + periodicDistance[1] * periodicDistance[1]);
            if (checking->isLeaf && periodicDist <= 1e-10) continue;
            if (!checking->isLeaf && checking->size / periodicDist >= OpeningAngle) {
                toCheck.push(checking->q1Child);
                toCheck.push(checking->q2Child);
                toCheck.push(checking->q3Child);
                toCheck.push(checking->q4Child);
                continue;
            }

            double Rho = periodicDist / ParticleRadius;

            std::tuple<double, double, double> Vector1 = {
                    Kernel(Rho) * periodicDistance[0] / (periodicDist * periodicDist),
                    Kernel(Rho) * periodicDistance[1] / (periodicDist * periodicDist),
                    0
            };
            auto [ResultX, ResultY, _] = Cross(Vector1, {0, 0, checking->particle.Vorticity});
            std::get<0>(Out[i]) -= ResultX;
            std::get<1>(Out[i]) -= ResultY;
            std::get<2>(Out[i]) += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (checking->particle.Vorticity - thisParticle.Vorticity);
        }
    }

    return Out;
}

void vpm::QuadTreeCodeNode::build(const std::vector<Particle>& particles, double _size, double _centerX, double _centerY) {
    this->size = _size;
    if (particles.empty()) {
        this->isLeaf = true;
        return;
    }

    this->isEmpty = false;
    if (particles.size() == 1) {
        this->particle = particles[0];
        this->isLeaf = true;
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

    this->particle = Particle(centroidX, centroidY, 0, 0, centerVorticity);

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

    q1Child = new vpm::QuadTreeCodeNode();
    q2Child = new vpm::QuadTreeCodeNode();
    q3Child = new vpm::QuadTreeCodeNode();
    q4Child = new vpm::QuadTreeCodeNode();

    double subSize = _size / 2;
    double subCenterOffset = subSize / 2;

    q1Child->build(particlesQ1, subSize, _centerX + subCenterOffset, _centerY + subCenterOffset);
    q2Child->build(particlesQ2, subSize, _centerX - subCenterOffset, _centerY + subCenterOffset);
    q3Child->build(particlesQ3, subSize, _centerX - subCenterOffset, _centerY - subCenterOffset);
    q4Child->build(particlesQ4, subSize, _centerX + subCenterOffset, _centerY - subCenterOffset);
}

vpm::QuadTreeCodeNode::~QuadTreeCodeNode() {
    delete q1Child;
    delete q2Child;
    delete q3Child;
    delete q4Child;
}
