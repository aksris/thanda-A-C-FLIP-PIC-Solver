#ifndef CUBE_H
#define CUBE_H
#include "geom.hpp"

class Cube: public Geometry
{
public:
    Cube();
    std::vector<GLfloat> g_cube_vertex_buffer_data;
    std::vector<GLuint> cub_idx;
    void push_vert_data();
    void push_indices_data();
};

#endif // CUBE_H
