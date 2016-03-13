#include "glm\glm.hpp"
#include "GL\glew.h"
#include <vector>
#include "../particles.h"

// class containing a container, defaulting to a cube.
// do NOT make containers until you have a context up and running!
class Container {
public:
	GLuint vertexBuffer;
	GLuint elementBuffer;

	glm::mat4 modelMatrix;

	std::vector<GLuint> gl_indices;
	std::vector<glm::vec3> gl_positions;

	bool collisionDetect(Particle *p);
};