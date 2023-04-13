#include <cmath>
#include <stack>
#include "treecode.h"

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
    tree.build(Particles, DomainL);

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
            auto [periodicDistanceX, periodicDistanceY] = thisParticle.PeriodicDistanceVector(checking->particle, DomainL);
            double periodicDist = sqrt(periodicDistanceX * periodicDistanceX + periodicDistanceY * periodicDistanceY);
            if (!checking->isLeaf && checking->size / periodicDist >= OpeningAngle) {
                toCheck.push(checking->q1Child);
                toCheck.push(checking->q2Child);
                toCheck.push(checking->q3Child);
                toCheck.push(checking->q4Child);
                continue;
            }

            double Rho = periodicDist / ParticleRadius;

            std::tuple<double, double, double> Vector1 = {
                    Kernel(Rho) * periodicDistanceX / (periodicDist * periodicDist),
                    Kernel(Rho) * periodicDistanceY / (periodicDist * periodicDist),
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

void vpm::QuadTreeCodeNode::build(const std::vector<Particle>& particles, double _size) {
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

    double centerX = 0;
    double centerY = 0;
    double centerVX = 0;
    double centerVY = 0;
    double centerVorticity = 0;
    for (auto _particle : particles) {
        centerX += std::get<0>(_particle.Position);
        centerY += std::get<1>(_particle.Position);
        centerVX += std::get<0>(_particle.Velocity);
        centerVY += std::get<1>(_particle.Velocity);
        centerVorticity += _particle.Vorticity;
    }

    centerX /= static_cast<double>(particles.size());
    centerY /= static_cast<double>(particles.size());
    centerVX /= static_cast<double>(particles.size());
    centerVY /= static_cast<double>(particles.size());
    centerVorticity /= static_cast<double>(particles.size());

    this->particle = Particle({centerX, centerY}, {centerVY, centerVX}, centerVorticity);

    std::vector<Particle> particlesQ1 = std::vector<Particle>();
    std::vector<Particle> particlesQ2 = std::vector<Particle>();
    std::vector<Particle> particlesQ3 = std::vector<Particle>();
    std::vector<Particle> particlesQ4 = std::vector<Particle>();
    for (auto _particle : particles) {
        double pX = std::get<0>(_particle.Position);
        double pY = std::get<1>(_particle.Position);
        if (pX >= centerX && pY >= centerY) particlesQ1.push_back(_particle);
        else if (pX < centerX && pY >= centerY) particlesQ2.push_back(_particle);
        else if (pX < centerX && pY < centerY) particlesQ3.push_back(_particle);
        else if (pX >= centerX && pY < centerY) particlesQ4.push_back(_particle);
    }

    q1Child = new vpm::QuadTreeCodeNode();
    q2Child = new vpm::QuadTreeCodeNode();
    q3Child = new vpm::QuadTreeCodeNode();
    q4Child = new vpm::QuadTreeCodeNode();

    q1Child->build(particlesQ1, _size / 2);
    q2Child->build(particlesQ2, _size / 2);
    q3Child->build(particlesQ3, _size / 2);
    q4Child->build(particlesQ4, _size / 2);
}

vpm::QuadTreeCodeNode::~QuadTreeCodeNode() {
    delete q1Child;
    delete q2Child;
    delete q3Child;
    delete q4Child;
}
