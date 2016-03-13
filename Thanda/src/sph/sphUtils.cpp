#include "sphUtils.h"

SPHUtils::SPHUtils() {
	particleSpacing = 0.05; 
	stiffness = 1000;  
	viscosity = 25;
	mass = 0.125; 
	restDensity = 1000;
	timeStep = 0.001;
	h = 0.100001;
}


float SPHUtils::clamp(float val, float min, float max) {
	if (val < min) return min;
	else if (val > max) return max;
	else return val;
}

// Applies the poly smoothing kernel to the given value h using the
// given distance dist. Used for density calculations.
float SPHUtils::kernelPoly(float h, float dist) {
	if (dist > h) return 0;
	float a = pow((h * h - dist * dist), 3);
	return 315 * a / (64 * M_PI * pow(h, 9));
}

// Applies the gradient of the spiky smoothing kernel to the given 
// value h using the given distance vector distvector and distance
// dist. Used for pressure force calculations.
glm::vec3 SPHUtils::kernelGradSpiky(float h, glm::vec3 distVector, float dist) {
	if (dist == 0 || dist > h) return glm::vec3(0, 0, 0);
	float a = -45 * pow((h - dist), 2);
	float b = dist * M_PI * pow(h, 6);
	return distVector * (a / b);
}

// Applies the laplacian of the viscous smoothing kernel to the given 
// value h using the given distance dist. Used for viscosity force
// calculations.
float SPHUtils::kernelLaplacianViscous(float h, float dist) {
	if (dist > h) return 0;
	return (45 * (h - dist) / (M_PI * pow(h, 6)));
}
