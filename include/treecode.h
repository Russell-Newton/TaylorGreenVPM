#pragma once

#include <utility>
#include <list>

#include "vpm.h"

namespace vpm {
    struct QuadTreeCodeNode {
        size_t q1ChildIdx;
        size_t q2ChildIdx;
        size_t q3ChildIdx;
        size_t q4ChildIdx;
        Particle particle;
        double size;
        bool isLeaf;
        bool isEmpty;

        QuadTreeCodeNode() : q1ChildIdx(-1), q2ChildIdx(-1), q3ChildIdx(-1), q4ChildIdx(-1), particle(), size(0), isLeaf(false), isEmpty(true) {};

        virtual ~QuadTreeCodeNode();

        void build(const std::vector<Particle>& particles, double _size, double _centerX, double _centerY, std::list<QuadTreeCodeNode>& tree);

        friend std::ostream& operator<<(std::ostream& os, const QuadTreeCodeNode& node) {
            return os << "Particle: " << node.particle << " Size: " << node.size << " isLeaf: " << node.isLeaf << " isEmpty: " << node.isEmpty << " q1ChildIdx: " << node.q1ChildIdx << " q2ChildIdx: " << node.q2ChildIdx << " q3ChildIdx: " << node.q3ChildIdx << " q4ChildIdx: " << node.q4ChildIdx;
        }
    };

    std::vector<std::tuple<double, double, double>> CalcDerivativeTreeCode(
            std::vector<Particle> Particles,
            double DomainL,
            double ParticleRadius,
            double Viscosity,
            double OpeningAngle = 0.5
    );
}
