#ifndef SCENE_H
#define SCENE_H

#include <glm/glm.hpp>
class Scene
{
public:
    Scene();

    glm::vec3 containerBounds; //scale of the container with <0, 0, 0> as center
    glm::vec3 particleBounds; //cube right now; <0, 0, 0>

    float particle_separation;
};

#endif // SCENE_H
