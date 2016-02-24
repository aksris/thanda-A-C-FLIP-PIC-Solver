#ifndef PARTICLES
#define PARTICLES
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
    glm::vec3 pos, speed;
    unsigned char r,g,b,a; // Color
    float size, angle, weight;
    float life; // Remaining life of the particle. if <0 : dead and unused.
    float cameradistance; // *Squared* distance to the camera. if dead : -1.0f

    bool operator<(const Particle& that) const {
        // Sort in reverse order : far particles drawn first.
        return this->cameradistance > that.cameradistance;
    }
};

class ParticleSystem{
public:
    ParticleSystem();

    int LastUsedParticle; int MaxParticles;
    std::vector<Particle> ParticlesContainer;

    int findUnusedParticles();
    void sortParticles();
    void particlesInit();
    void genParticles(float particle_separation, float boundx, float boundy, float boundz);
//    void simulateParticles(GLfloat *g_particule_position_size_data, GLubyte *g_particule_color_data, vec3 CameraPosition, double delta);
};
#endif // PARTICLES

