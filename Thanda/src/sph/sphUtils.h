#include "glm\glm.hpp"
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

class SPHUtils {
	public:
		float particleSpacing; // h 
		float stiffness;  // k
		float viscosity; // mu
		float mass; // particle mass
		float restDensity;
		float timeStep;
		float h;
		static float clamp(float val, float min, float max);
		static float kernelPoly(float h, float dist);
		static glm::vec3 kernelGradSpiky(float h, glm::vec3 distVector, float dist);
		static float kernelLaplacianViscous(float h, float dist);

		SPHUtils();
};