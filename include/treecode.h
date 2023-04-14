#pragma once

#include <utility>

#include "vpm.h"

namespace vpm {
    struct QuadTreeCodeNode {
        QuadTreeCodeNode *q1Child;
        QuadTreeCodeNode *q2Child;
        QuadTreeCodeNode *q3Child;
        QuadTreeCodeNode *q4Child;
        Particle particle;
        double size;
        bool isLeaf;
        bool isEmpty;

        QuadTreeCodeNode() : q1Child(nullptr), q2Child(nullptr), q3Child(nullptr), q4Child(nullptr), particle(), size(0), isLeaf(false), isEmpty(true) {};

        virtual ~QuadTreeCodeNode();

        void build(const std::vector<Particle>& particles, double _size, double _centerX, double _centerY);
    };

    std::vector<std::tuple<double, double, double>> CalcDerivativeTreeCode(
            std::vector<Particle> Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity,
            double OpeningAngle = 0.5
    );
}
