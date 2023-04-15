#include <cmath>
#include <iostream>
#include "vpm.h"
#include <fstream>
#include <cmath>

#define MYSTERY_FRACTION 0.1
#define PLOT 0


void plotField(int timeStep, vpm::Particle* Particles, double L, double ParticleRadius, int N,  size_t plotResolution) {

/*    
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
            std::tuple<double, double> velAtPoint = vpm::CalcVelAtPoint(xPoint, yPoint, Particles, L, ParticleRadius, N);
            outputVTK <<  std::get<0>(velAtPoint) << " " << std::get<1>(velAtPoint) << std::endl;
        }

    outputVTK << "           </DataArray>\n";
    outputVTK << "      </PointData>\n";
    outputVTK << "      <CellData>\n";
    outputVTK << "      </CellData>\n";
    outputVTK << "  </Piece>\n";
    outputVTK << " </ImageData>\n";
    outputVTK << "</VTKFile>\n";

    outputVTK.close();


    std::ofstream outputVTP;
    outputVTP.open("vtk_out/particles_"+std::to_string(timeStep)+".vtp");

    outputVTP << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    outputVTP << " <PolyData>\n";
    outputVTP << "  <Piece NumberOfPoints=\"" << N << "\" NumberOfVerts=\"" << N << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"; 
    outputVTP << "      <PointData>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"2\" format=\"ascii\">\n";

    for (int i = 0; i < N; i++) {
        vpm::Particle part = Particles[i];
        outputVTP << std::get<0>(part.Velocity) << " " << std::get<1>(part.Velocity) << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "      </PointData>\n";
    outputVTP << "      <CellData>\n";
    outputVTP << "      </CellData>\n";
    outputVTP << "      <Points>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";

    for (int i = 0; i < N; i++) {
        vpm::Particle part = Particles[i];
        outputVTP << std::get<0>(part.Position) << " " << std::get<1>(part.Position) << " 0.0 " << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "      </Points>\n";

    outputVTP << "      <Verts>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";

    for (int i = 0; i < N; i++) {
        outputVTP << i << " " << std::endl;
    }

    outputVTP << "           </DataArray>\n";
    outputVTP << "          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << N << "\" RangeMax=\"" << N << "\">\n";

    for (int i = 0; i < N; i++) {
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
*/
    
}


int main() {
    // sim params
    double L = 2*M_PI; // meters
    size_t Resolution = 16; // initialize particles as 64*64 square grid
    auto ResolutionDouble = static_cast<double>(Resolution);
    double Viscosity = 1E-2; // kinematic Viscosity m^2/s
    double dt = 1E-1;
    size_t nt = 100;
    double dX, dY, dOmega;
    int N = Resolution * Resolution;

    double ParticleRad = L / ResolutionDouble * 2;
    double ParticleVol = M_PI * ParticleRad * ParticleRad;


    #pragma acc enter data copyin(L, Viscosity, ParticleRad, dt, N)

    
    
    vpm::Particle Particles[N];
    #pragma acc enter data create(Particles[0:N])

    std::tuple<double, double, double> Derivatives[N];
    #pragma acc enter data create(Derivatives[0:N])

    // plot params
   //size_t PlotResolution = 16;

    // Initialize Particles
    for (size_t i = 0; i < N; i++) {
        double x = static_cast<double>(i % Resolution) * ParticleRad / 2.0 + ParticleRad / 4.0;
        double y = static_cast<double>(i - i % Resolution) / ResolutionDouble * (L / ResolutionDouble) + ParticleRad / 4.0;
        double u = std::cos(x) * std::sin(y);
        double v = -std::sin(x) * std::cos(y);
        double omega = -2 * std::cos(x) * std::cos(y);
        Particles[i] = {std::make_tuple(x, y), std::make_tuple(u, v), omega * ParticleVol * MYSTERY_FRACTION};
        Derivatives[i] = std::make_tuple(0.0, 0.0, 0.0);
    }
    #pragma acc update device(Particles[0:N])
    #pragma acc update device(Derivatives[0:N])

    // plot

#ifdef PLOT
//    plotField(0, Particles, L, ParticleRad, N, PlotResolution);
#endif


    for (size_t t = 1; t < nt; t++) {
        
        vpm::CalcDerivative(Particles, L, ParticleRad, Viscosity, N, Derivatives);

        #pragma acc parallel loop gang vector default(present) private(dX, dY, dOmega)
        for (size_t i = 0; i < N; i++) {

            std::make_tuple(dX, dY, dOmega) = Derivatives[i];
            std::get<0>(Particles[i].Position) += dX * dt;
            std::get<1>(Particles[i].Position) += dY * dt;
            std::get<0>(Particles[i].Velocity) = dX;
            std::get<1>(Particles[i].Velocity) = dY;
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

        std::cout << t << std::endl;

        #pragma acc update host(Particles[0:N])


#ifdef PLOT
//        plotField(t, Particles, L, 0.5, N,  PlotResolution);
#endif
    }

    #pragma acc exit data delete(Derivatives[0:N])
    #pragma acc exit data delete(Particles[0:N])
    #pragma acc exit data delete(L, Viscosity, dt)

}
