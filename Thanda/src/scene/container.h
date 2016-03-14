#ifndef CONTAINER_H
#define CONTAINER_H

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

	Container();
	~Container();
	virtual bool collisionDetect(Particle *p, float dt, glm::vec3 &coll_Pos, glm::vec3 &coll_Nor);
};

// helpers for cube intersection
bool inUnitCube(glm::vec3 point);
float rayPlaneISX(glm::vec3 pos, glm::vec3 dir, glm::vec3 planePos, glm::vec3 norm);
bool nearlyEqual(float a, float b, float epsilon);
void checkSlab(glm::vec3 pos, glm::vec3 dir, float &nearMax, glm::vec3 &nearNorm, float &farMin, glm::vec3 &farNorm, glm::vec3 norm, glm::vec3 pos1, glm::vec3 pos2);

#endif