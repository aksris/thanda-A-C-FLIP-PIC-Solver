//
//  geom.hpp
//  Thanda

#ifndef geom_hpp
#define geom_hpp
#include "GL/glew.h"
#include <vector>
class Geometry{
public:
    Geometry();
    std::vector<GLfloat> g_cube_vertex_buffer_data;
};

#endif /* geom_hpp */
