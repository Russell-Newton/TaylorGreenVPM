#include <cmath>
#include <iostream>
#include <ranges>
#include "vpm.h"
#include <fstream>
#include "treecode.h"

#define MYSTERY_FRACTION 0.1
#define PLOT 1

void plotField(int timeStep, std::vector<vpm::Particle> Particles, double L, double ParticleRadius, size_t plotResolution) {

    std::ofstream outputVTK;
    outputVTK.open("vtk_out/vel_"+std::to_string(timeStep)+".vti");

    outputVTK << "<VTKFile type=\"ImageData\" version=\"2.2\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    outputVTK << " <ImageData WholeExtent=\"0 " << plotResolution << " 0 " << plotResolution << " 0 0\" Origin=\"0 0 0\" Spacing=\"" << L/plotResolution << " " << L/plotResolution << " " << L/plotResolution << "\">\n";
    outputVTK << "  <Piece Extent=\"0 " << plotResolution << " 0 " << plotResolution << " 0 0\">\n";
    outputVTK << "      <PointData>\n";
    outputVTK << "          <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"2\" format=\"ascii\">\n";

    for (int iY = 0; iY <= plotResolution; iY++)
        for (int iX = 0; iX <= plotResolution; iX++)
        {
            double xPoint = iX * L/plotResolution;
            double yPoint = iY * L/plotResolution;
            std::tuple<double, double> velAtPoint = vpm::CalcVelAtPoint(xPoint, yPoint, Particles, L, ParticleRadius);
            outputVTK <<  std::get<0>(velAtPoint) << " " << std::get<1>(velAtPoint) << std::endl;
        }

    outputVTK << "           </DataArray>\n";
    outputVTK << "      </PointData>\n";
    outputVTK << "      <CellData>\n";
    outputVTK << "      </CellData>\n";
    outputVTK << "  </Piece>\n";
    outputVTK << " </ImageData>\n";
    outputVTK << "</VTKFile>\n";

    std::ofstream outputVTP;
    outputVTP.open("vtk_out/particles_"+std::to_string(timeStep)+".vtp");

    outputVTP << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    outputVTP << " <PolyData>\n";
    outputVTP << "  <Piece NumberOfPoints=\"" << Particles.size() << "\" NumberOfVerts=\"" << Particles.size() << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"; 
    outputVTP << "      <PointData>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"2\" format=\"ascii\">\n";

    for (vpm::Particle& part : Particles) {
        outputVTP << part.VelocityX << " " << part.VelocityY << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "      </PointData>\n";
    outputVTP << "      <CellData>\n";
    outputVTP << "      </CellData>\n";
    outputVTP << "      <Points>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";

    for (vpm::Particle& part : Particles) {
        outputVTP << part.PositionX << " " << part.PositionY << " 0.0 " << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Points>\n";

    outputVTP << "      <Verts>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";

    for (int i = 0; i < Particles.size(); i++) {
        outputVTP << i << " " << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << Particles.size() << "\" RangeMax=\"" << Particles.size() << "\">\n";

    for (int i = 0; i < Particles.size(); i++) {
        outputVTP << i+1 << " " << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Verts>\n";

    outputVTP << "      <Lines>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Lines>\n";

    outputVTP << "      <Strips>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Strips>\n";

    outputVTP << "      <Polys>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n";
    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Polys>\n";
    outputVTP << "  </Piece>\n";
    outputVTP << " </PolyData>\n";
    outputVTP << "</VTKFile>\n";
}

int main() {
    // sim params
    double L = 2*M_PI; // meters
    size_t Resolution = 64; // initialize particles as 64*64 square grid
    auto ResolutionDouble = static_cast<double>(Resolution);
    double Viscosity = 1E-2; // kinematic Viscosity m^2/s
    double dt = 1E-1;
    size_t nt = 100;

    size_t N = Resolution * Resolution;
    double ParticleRad = L / ResolutionDouble * 2;
    double ParticleVol = M_PI * ParticleRad * ParticleRad;

    std::vector<vpm::Particle> Particles(N);

    // plot params
    size_t PlotResolution = 16;

    // Initialize Particles
    for (size_t i = 0; i < N; i++) {
        double x = static_cast<double>(i % Resolution) * ParticleRad / 2.0 + ParticleRad / 4.0;
        double y = static_cast<double>(i - i % Resolution) / ResolutionDouble * (L / ResolutionDouble) + ParticleRad / 4.0;
        double u = std::cos(x) * std::sin(y);
        double v = -std::sin(x) * std::cos(y);
        double omega = -2 * std::cos(x) * std::cos(y);
        Particles[i] = vpm::Particle(x, y, u, v, omega * ParticleVol * MYSTERY_FRACTION);
    }

    // plot
#ifdef PLOT
    plotField(0, Particles, L, ParticleRad, PlotResolution);
#endif

    for (size_t t = 1; t < nt; t++) {
        auto Derivatives = vpm::CalcDerivativeTreeCode(Particles, L, ParticleRad, Viscosity, 0.2);
        for (size_t i = 0; i < N; i++) {
            vpm::Particle& Particle = Particles[i];
            auto [dX, dY, dOmega] = Derivatives[i];
            Particle.PositionX += dX * dt;
            Particle.PositionY += dY * dt;
            Particle.VelocityX = dX;
            Particle.VelocityY = dY;
            Particle.Vorticity += dOmega * dt;

            if (Particle.PositionX < 0) {
                Particle.PositionX += L;
            }
            if (Particle.PositionY < 0) {
                Particle.PositionY += L;
            }
            if (Particle.PositionX > L) {
                Particle.PositionX -= L;
            }
            if (Particle.PositionY > L) {
                Particle.PositionY -= L;
            }
        }

        std::cout << t << std::endl;
#ifdef PLOT
        plotField(t, Particles, L, 0.5, PlotResolution);
#endif
    }
}
