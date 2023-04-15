#include <cmath>
#include <stack>
#include <cstdio>
#include <ctime>
#include "treecode.h"
#include <iostream>

std::vector<std::tuple<double, double, double>> vpm::CalcDerivativeTreeCode(
            std::vector<Particle> Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity,
            double OpeningAngle,
            bool printTime
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

    clock_t start = clock();

    for (size_t i = 0; i < N; i++) {
        vpm::Particle thisParticle = Particles[i];
        // get contribution from tree
        std::stack<size_t> toCheck = std::stack<size_t>();
        toCheck.push(0);

        while (!toCheck.empty()) {
            QuadTreeCodeNode checking = treeAsVec[toCheck.top()];
            toCheck.pop();

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

            std::tuple<double, double, double> Vector1 = {
                    Kernel(Rho) * periodicDistance[0] / (periodicDist * periodicDist),
                    Kernel(Rho) * periodicDistance[1] / (periodicDist * periodicDist),
                    0
            };
            auto [ResultX, ResultY, _] = Cross(Vector1, {0, 0, checking.particle.Vorticity});
            std::get<0>(Out[i]) -= ResultX;
            std::get<1>(Out[i]) -= ResultY;
            std::get<2>(Out[i]) += (2.0 * Viscosity / ParticleRadius2) * (1.0 / ParticleRadius2) * ViscousKernel(Rho) * ParticleVol * (checking.particle.Vorticity - thisParticle.Vorticity);
        }
    }

    clock_t end = clock();

    if (printTime) std::cout << "Time this step: " << static_cast<double>(end - start) / CLOCKS_PER_SEC << "s" << std::endl;

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
