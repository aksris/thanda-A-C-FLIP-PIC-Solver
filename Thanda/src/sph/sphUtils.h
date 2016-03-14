#include "glm\glm.hpp"
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

// defines for SPH
#define SPH_PARTICLE_SPACING 0.05f // h
#define SPH_STIFFNESS 1000.0f // k
#define SPH_VISCOSITY 25.0f // mu
#define SPH_MASS 0.125f // particle mass
#define SPH_REST_DENSITY 1000.0f
#define SPH_TIMESTEP 0.04f//0.001f
#define SPH_H 0.100001f // neighbor search radius (aka grid cell size)

float clamp(float val, float min, float max);
float kernelPoly(float h, float dist);
glm::vec3 kernelGradSpiky(float h, glm::vec3 distVector, float dist);
float kernelLaplacianViscous(float h, float dist);