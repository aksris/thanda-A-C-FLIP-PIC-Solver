//
//  fluidSolver.hpp
//  Thanda

#ifndef fluidSolver_hpp
#define fluidSolver_hpp
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include "macgriddata.h"
#include "../scene/scene.hpp"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unordered_map>

using namespace glm;

//for marking cells
enum geomtype {AIR = 0, FLUID = 1, SOLID = 2};

class Particle{

public:
    Particle();
    ivec3 gridIdx;
    glm::vec3 pos, speed;
    unsigned char r,g,b,a; // Color
    float size;
    float cameradistance; // *Squared* distance to the camera. if dead : -1.0f

    bool operator<(const Particle& that) const {
        // Sort in reverse order : far particles drawn first.
        return this->cameradistance > that.cameradistance;
    }
};

class MACGrid{
public:
    MACGrid(const ivec3& resolution, const vec3& containerBounds, float cellSize);
    ~MACGrid();
    void initialize();
    MACGrid& operator=(const MACGrid& val);

    /*The grid data structure needs to be on the heap for maximum memory usage in calculating pressure (for MICCG(0)) */

    //Velocity grid
    MACGridDataX* vel_U;
    MACGridDataY* vel_V;
    MACGridDataZ* vel_W;

    //A separate data structure for storing FLIP velocities, could be optimized withing vel_U, vel_V and vel_W grids
    MACGridDataX* flip_vel_U;
    MACGridDataY* flip_vel_V;
    MACGridDataZ* flip_vel_W;

    MACGridDataX* save_kernel_wt_U;
    MACGridDataY* save_kernel_wt_V;
    MACGridDataZ* save_kernel_wt_W;
    //pressure grid; size is same as velocity grid
    MACGridData* P;

};

class FluidSolver{
public:
    FluidSolver(const Scene &iScene);
    ~FluidSolver();
    MACGrid* grid;

    float delta;
    int num_cells;
    vec3 containerBounds;
    ivec3 resolution;
    float cellSize;

    Scene scene;

    int MaxParticles;
    std::vector<Particle> ParticlesContainer;

    std::vector<Particle> particle_save;
    std::vector<Particle> particle_save_pic;//to save particle velocity; cannot reuse particle_save

    std::map<int, ivec3> reposition_map;

    /*
     Constructs a data structure with the number of cells based on the scene variables:
     resolution and bounds of the simulation container.
    */
    void constructMACGrid();

    /*
      Splatting particle velocity to grid by interpolating surrounding grid cells.
      9 cells around every i,j,k are set for trilinear interpolation to get velocities back to particles
    */
    void storeParticleVelocityToGrid();
    /*
      Storing current at time step n velocities on the grid, to then calculate FLIP and PIC vel
    */
    void storeCurrentGridVelocities();
    /*
      Subtract pressure gradient for setting up Pressure calculation matrix
    */
    void SubtractPressureGradient();

    /*
      Particle step function; that calls all functions to solve Navier-Stokes equation
    */
    void step(const float &dt);
    /*zero out velocity grid for every iteration*/
    void clearGrid();

    /*
      Marker grid for extrapolation; start by marking every cell that is not already fluid;
      i.e., mark all near fluid cells as fluid
    */
    void initializeMarkerGrid();
    /*
      Set coefficients for A matrix in the pressure solve; Actual pressure solve done by Eigen
    */
    void buildMatrixA(std::vector<Eigen::Triplet<double> > &coefficients, long n);
    /*
      Caclulate divergence for each direction: U, V, W
    */
    void buildDivergences(Eigen::VectorXd& u);
    void fillPressureGrid(Eigen::VectorXd x);

    /*add gravity forces to cell velocity*/
    void CalculateGravityToCell(float delta);
    /*Helper function for matrix A coefficient*/
    void insertCoefficient(int id, int i, int j, int k, double w, std::vector<Eigen::Triplet<double> > &coeffs);

    void calculateNewGridVelocities();
    /*reset boundary cells velocities to zero as they do not interact in pressure solve*/
    void setBoundaryVelocitiesToZero();
    /*
      Different velocity interpolation based on PIC and FLIP
    */
    void FlipSolve();
    void PicSolve();

    void ProjectPressure();

    /*advect using RK2*/
    vec3 integratePos(const vec3& pos, const vec3& speed, const float& time_step, bool RK2);

    /*trilinear interpolation*/
    void ExtrapolateVelocity();

    void genParticles();

};
#endif /* fluidSolver_hpp */
