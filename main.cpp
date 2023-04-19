#include <cmath>
#include <iostream>
#include "vpm.h"
#include <fstream>
#include "treecode.h"
#include <cmath>
#include <math.h>

#define MYSTERY_FRACTION 0.1
#define PLOT 1


void plotField(int timeStep, vpm::Particle* Particles, double L, double ParticleRadius, int N,  size_t plotResolution) {

    
    FILE* fp;
    
    std::string s1 = "vtk_out/vel_"+ std::to_string(timeStep)+ ".vti";
    fp = fopen(s1.c_str(), "w"); 

    fprintf(fp,"<VTKFile type=\"ImageData\" version=\"2.2\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp," <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\" %f %f %f\">\n", plotResolution, plotResolution, L/plotResolution, L/plotResolution, L/plotResolution);

    fprintf(fp, "  <Piece Extent=\"0 %d 0 %d 0 0\">\n", plotResolution, plotResolution);
    fprintf(fp,"      <PointData>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"2\" format=\"ascii\">\n");

    for (int iY = 0; iY <= plotResolution; iY++)
        for (int iX = 0; iX <= plotResolution; iX++)
        {
            double xPoint = iX * L/plotResolution;
            double yPoint = iY * L/plotResolution;
            std::tuple<double, double> velAtPoint = vpm::CalcVelAtPoint(xPoint, yPoint, Particles, L, ParticleRadius, N);
            fprintf(fp, "%f %f\n", std::get<0>(velAtPoint), std::get<1>(velAtPoint));
        }

    fprintf(fp, "           </DataArray>\n");
    fprintf(fp,  "      </PointData>\n");
    fprintf(fp, "      <CellData>\n");
    fprintf(fp, "      </CellData>\n");
    fprintf(fp, "  </Piece>\n");
    fprintf(fp, " </ImageData>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);

    std::string s2 = "vtk_out/particles_"+std::to_string(timeStep)+".vtp";
    fp = fopen(s2.c_str(), "w");

    fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp," <PolyData>\n");
    fprintf(fp,"  <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n", N, N); 
    fprintf(fp,"      <PointData>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"2\" format=\"ascii\">\n");

    for (int i = 0; i < N; i++) {
        vpm::Particle part = Particles[i];
         fprintf(fp,"%f %f\n", part.VelocityX , part.VelocityY);
   }

    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </PointData>\n");
    fprintf(fp,"      <CellData>\n");
    fprintf(fp,"      </CellData>\n");
    fprintf(fp,"      <Points>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n");

    for (int i = 0; i < N; i++) {
        vpm::Particle part = Particles[i];
        fprintf(fp,"%f %f 0.0\n", part.PositionX , part.PositionY);
   }

    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </Points>\n");

    fprintf(fp,"      <Verts>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n");

    for (int i = 0; i < N; i++) {
        fprintf(fp,"%d\n", i);
    }

    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", N, N);

    for (int i = 0; i < N; i++) {
        fprintf(fp,"%d\n", i+1);
    }

    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </Verts>\n");

    fprintf(fp,"      <Lines>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </Lines>\n");

    fprintf(fp,"      <Strips>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </Strips>\n");

    fprintf(fp,"      <Polys>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"          <DataArray type=\"Float64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\">\n");
    fprintf(fp,"           </DataArray>\n");
    fprintf(fp,"      </Polys>\n");
    fprintf(fp,"  </Piece>\n");
    fprintf(fp," </PolyData>\n");
    fprintf(fp,"</VTKFile>\n");

    fclose(fp);
    

    
}


int main() {
    // sim params
    double L = 2*M_PI; // meters
    size_t Resolution = 192; // initialize particles as 64*64 square grid
    auto ResolutionDouble = static_cast<double>(Resolution);
    double Viscosity = 1E-3; // kinematic Viscosity m^2/s
    double dt = 1E-1;
    size_t nt = 100;
    double dX, dY, dOmega;
    int N = Resolution * Resolution;

    double ParticleRad = L / ResolutionDouble * 2.0;
    double ParticleVol = 3.14 * ParticleRad * ParticleRad;

    int GPUGridLevel = 16;
    int numNodes = GPUGridLevel*GPUGridLevel;
    double gridSize = L/GPUGridLevel;

    #pragma acc enter data copyin(L, Viscosity, ParticleRad, dt, N, numNodes, gridSize)

    vpm::BoxNodeGPU BoxNodes[numNodes];    
    #pragma acc enter data create(BoxNodes[0:numNodes])
    
    vpm::Particle Particles[N];
    #pragma acc enter data create(Particles[0:N])

    vpm::Particle ParticlesRegrid[N];
    #pragma acc enter data create(ParticlesRegrid[0:N])

    std::tuple<double, double, double> Derivatives[N];
    #pragma acc enter data create(Derivatives[0:N])

    // plot params
   size_t PlotResolution = 16;

    // Initialize Particles
    for (size_t i = 0; i < N; i++) {
        double x = static_cast<double>(i % Resolution) * ParticleRad / 2.0 + ParticleRad / 4.0;
        double y = static_cast<double>(i - i % Resolution) / ResolutionDouble * (L / ResolutionDouble) + ParticleRad / 4.0;
        double u = std::cos(x) * std::sin(y);
        double v = -std::sin(x) * std::cos(y);
        double omega = -2 * std::cos(x) * std::cos(y);
        
        Particles[i] = {x, y, u, v, omega * ParticleVol * MYSTERY_FRACTION};
        Derivatives[i] = std::make_tuple(0.0, 0.0, 0.0);
    }
    #pragma acc update device(Particles[0:N])
    #pragma acc update device(Derivatives[0:N])

    
    // plot

#ifdef PLOT
    plotField(0, Particles, L, ParticleRad, N, PlotResolution);
#endif

    clock_t start = clock();

    for (size_t t = 1; t < nt; t++) {

        // vpm::CalcDerivative(Particles, L, ParticleRad, Viscosity, N, Derivatives);
        CalcDerivativeGPUTree(BoxNodes, Particles, L, ParticleRad, Viscosity, N, numNodes, 0.2, Derivatives);

        #pragma acc parallel loop gang vector default(present) private(dX, dY, dOmega)
        for (int i = 0; i < N; i++) {

            dX = std::get<0>(Derivatives[i]);
            dY = std::get<1>(Derivatives[i]);
            dOmega = std::get<2>(Derivatives[i]);

            Particles[i].PositionX += dX * dt;
            Particles[i].PositionY += dY * dt;
            Particles[i].VelocityX = dX;
            Particles[i].VelocityY = dY;
            Particles[i].Vorticity += dOmega * dt;

            if (Particles[i].PositionX < 0) {
                Particles[i].PositionX += L;
            }
            if (Particles[i].PositionY < 0) {
                Particles[i].PositionY += L;
            }
            if (Particles[i].PositionX > L) {
                Particles[i].PositionX -= L;
            }
            if (Particles[i].PositionY > L) {
                Particles[i].PositionY -= L;
            }
        }

        std::cout << t << std::endl;
        
//         #pragma acc update host(Particles[0:N])


// #ifdef PLOT
//        plotField(t, Particles, L, 0.5, N,  PlotResolution);
// #endif
    }

    clock_t end = clock();
    printf("Time: %f\n", double(end-start)/CLOCKS_PER_SEC);


#ifdef PLOT
        #pragma acc update host(Particles[0:N])
       plotField(99, Particles, L, 0.5, N,  PlotResolution);
#endif

    #pragma acc exit data delete(Derivatives[0:N])
    #pragma acc exit data delete(Particles[0:N])
    #pragma acc exit data delete(BoxNodes[0:numNodes])
    #pragma acc exit data delete(L, Viscosity, dt)

}
