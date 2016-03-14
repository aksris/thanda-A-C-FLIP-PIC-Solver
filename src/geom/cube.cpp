#include "cube.h"

Cube::Cube()
{
    this->push_vert_data();
    this->push_indices_data();
}

void Cube::push_vert_data(){
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
}

void Cube::push_indices_data(){
    cub_idx.push_back(0);
    cub_idx.push_back(1);
    cub_idx.push_back(1);
    cub_idx.push_back(2);
    cub_idx.push_back(2);
    cub_idx.push_back(3);
    cub_idx.push_back(3);
    cub_idx.push_back(0);

    cub_idx.push_back(0);
    cub_idx.push_back(4);
    cub_idx.push_back(1);
    cub_idx.push_back(5);
    cub_idx.push_back(2);
    cub_idx.push_back(6);
    cub_idx.push_back(3);
    cub_idx.push_back(7);

    cub_idx.push_back(4);
    cub_idx.push_back(5);
    cub_idx.push_back(5);
    cub_idx.push_back(6);
    cub_idx.push_back(6);
    cub_idx.push_back(7);
    cub_idx.push_back(7);
    cub_idx.push_back(4);
}
