#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include "glm\glm.hpp"
#include "../particles.h"
#include "../scene/container.h"
#include "sphUtils.h"

class SPHSolver {
public:
	Container *container;
	ParticleSystem *particles;
	int numParticles;

	SPHSolver(Container *container, ParticleSystem *particles);
	void step();
private:
	void stepSingle(); // run a neighbor search and do computations per neighbor
	float computeDensityContrib(Particle *p, Particle *neighbor);
};

#endif