//
//  viewer.hpp
//  Thanda

#ifndef viewer_hpp
#define viewer_hpp

#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "../scene/scene.hpp"
#include <stdio.h>
#include <string>
#include <vector>
#include "../fluidSolver/fluidSolver.hpp"
#include "../camera/camera.hpp"
#include <fstream>
#include <algorithm>

using namespace std;
#include "../geom/cube.h"
#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <openvdb/openvdb.h>
#include <openvdb_points/openvdb.h>
#include <openvdb_points/tools/PointDataGrid.h>
#include <openvdb_points/tools/PointConversion.h>
#include <openvdb_points/tools/PointCount.h>
#endif


class Viewer {
public:
    Viewer(int width, int height, const Scene& s);
    ~Viewer();
    void initializeGL();
    void initializeShader();
//    void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);

    void drawCube();
    void drawParticles(int ParticlesCount);
    void allocateParticleBuffers(GLuint billboard_vertex_buffer, GLuint particles_position_buffer, GLuint particles_color_buffer);
    void display();

    int w, h;
    Scene scene;
    Cube cube;
    Camera camera;
    FluidSolver *fluid;
    GLuint CameraRight_worldspace_ID, CameraUp_worldspace_ID,  ViewProjMatrixID, programID, TextureID, programIDGeometry, geomMatrixID;
    GLuint elementbuffer, cubevertexbuffer;
    GLubyte* g_particule_color_data;

    GLFWwindow* window;

    glm::mat4 ProjectionMatrix;
    glm::mat4 ViewMatrix;
    // We will need the camera's position in order to sort the particles
    // w.r.t the camera's distance.
    // There should be a getCameraPosition() function in common/controls.cpp,
    // but this works too.
    glm::vec3 CameraPosition;

#ifdef __APPLE__
    // Create some point positions
    std::vector<openvdb::Vec3f> positions;
#endif


};
#endif /* viewer_hpp */
