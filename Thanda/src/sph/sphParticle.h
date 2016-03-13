#ifndef SPHPARTICLES
#define SPHPARTICLES
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

#include "../particles.h"

using namespace glm;

class SPHParticle : public Particle { 

public:
    SPHParticle();
    float density, pressure;
};

#endif // SPHPARTICLES