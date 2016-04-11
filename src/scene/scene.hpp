//
//  scene.hpp
//  Thanda
#ifndef SCENE_H
#define SCENE_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <GL/glew.h>
#include <json/json.h>
#include <string.h>
#include <glm/glm.hpp>
#include <tinyobj/tiny_obj_loader.h>

class Scene
{
public:
    Scene();

    glm::vec3 containerBounds; //scale of the container with <0, 0, 0> as center
    glm::vec3 particleBounds; //cube right now; <0, 0, 0>

    std::vector<glm::vec3> positions;

    glm::ivec3 resolution;

    float particle_separation;

    void parseScene(const char* filename, Scene& scene);
    GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path);
    void LoadOBJ(const char* filename, Scene& scene);
    GLuint loadDDS(const char * imagepath);
    GLuint loadBMP_custom(const char * imagepath);
};

#endif // SCENE_H
