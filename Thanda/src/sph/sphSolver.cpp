#include "sphSolver.h"

SPHSolver::SPHSolver(Container *container, ParticleSystem *particles) {
	this->container = container;
	this->particles = particles;
	numParticles = particles->ParticlesContainer.size();
}

void SPHSolver::step() {
	// testing collision detection. TODO: remove
	glm::vec3 nor;
	glm::vec3 pos;

	float dt = SPH_TIMESTEP;

	for (int i = 0; i < numParticles; i++) {
		Particle &p = particles->ParticlesContainer[i]; // shortcut
		p.speed += glm::vec3(0.0f, -9.81f, 0.0f) * dt * 0.5f;
		p.pos += p.speed * dt;

		// testing collision detection. TODO: remove
		if (container->collisionDetect(&p, dt, pos, nor)) {
			p.r = abs(nor.r) * 220;
			p.g = abs(nor.g) * 220;
			p.b = abs(nor.b) * 220;
		}
		else {
			p.r = 0;
			p.g = 0;
			p.b = 220;
		}
	}
}

void SPHSolver::stepSingle(Particle *p) {
	// reset particle's pressure and density
	p->pressure = 0.0f;
	p->density = 0.0f;

	// do a neighbor search
	std::vector<Particle*> neighbors;
	naiveNeighborSearch(p, neighbors);
}

void SPHSolver::naiveNeighborSearch(Particle *p, std::vector<Particle*> &neighbors) {
	for (int i = 0; i < numParticles; i++) {
		float dist = glm::length(p->pos - particles->ParticlesContainer[i].pos);
		if (dist < SPH_H) {
			neighbors.push_back(&particles->ParticlesContainer[i]);
		}
	}
}