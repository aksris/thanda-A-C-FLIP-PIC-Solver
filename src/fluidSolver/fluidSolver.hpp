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
#include <iostream>
using namespace glm;

class Particle{

public:
    Particle();
    int gridPos;
    glm::vec3 pos, speed;
    unsigned char r,g,b,a; // Color
    float size, angle, mass, density;
    float life; // Remaining life of the particle. if <0 : dead and unused.
    float cameradistance; // *Squared* distance to the camera. if dead : -1.0f

    bool operator<(const Particle& that) const {
        // Sort in reverse order : far particles drawn first.
        return this->cameradistance > that.cameradistance;
    }
};

class FluidSolver{
public:
    FluidSolver();

    std::vector<float> forces; // at the center of grid
    std::vector<float> vel_u; //u,v,w components of velocity
    std::vector<float> vel_v;
    std::vector<float> vel_w;

    std::vector<std::vector<Particle*>> hash_grid; //improve here for faster particle to grid access
    std::vector<float> vel_u_current;
    std::vector<float> vel_v_current;
    std::vector<float> vel_w_current;

    std::vector<float> vel_u_change;
    std::vector<float> vel_v_change;
    std::vector<float> vel_w_change;



    int i_size;
    int j_size;

    int LastUsedParticle; int MaxParticles;
    std::vector<Particle> ParticlesContainer;

    std::vector<Particle> particle_save; //to save particle velocity

    //cell index returner
    int findGridIndex(int x, int y, int z);

    void constructMACGrid(glm::vec3 containerBounds);
    void initMACGrid(Particle &p);
    void storeParticleVelocityToGrid();
    void storeCurrentGridVelocities();

    void projectPressure();
    void calculateNewGridVelocities();
    void setBoundaryVelocitiesToZero(const glm::vec3 containerBounds);
    void flipSolve();

    void calculateDensity(Particle &p);
    void calculateGravityForces(Particle& p);
    void ExtrapolateVelocity();
    int findUnusedParticles();
    void sortParticles();
    void particlesInit();
    void genParticles(float particle_separation, float boundx, float boundy, float boundz);

};
#endif /* fluidSolver_hpp */
