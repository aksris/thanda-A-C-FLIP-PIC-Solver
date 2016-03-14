//
//  fluidSolver.cpp
//  Thanda


#include "fluidSolver.hpp"
Particle::Particle(){
    pos = glm::vec3(0.f, 0.f, 0.f);
    speed = glm::vec3(0.f, 0.f, 0.f);
    r = 0;
    g = 0;
    b = 220;
    a = 230;
    size = 0.1f;
    angle = 45.f;
    mass = 1.f;
    life = 1.f;
    cameradistance = 10.f;
    density = 1.f;
}

FluidSolver::FluidSolver(){
    LastUsedParticle = 0;
    MaxParticles = 200000;
}

int FluidSolver::findUnusedParticles(){

    for(int i=LastUsedParticle; i<MaxParticles; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    for(int i=0; i<LastUsedParticle; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    return 0; // All particles are taken, override the first one
}

void FluidSolver::sortParticles(){
    std::sort(&ParticlesContainer[0], &ParticlesContainer[MaxParticles]);
}

void FluidSolver::constructMACGrid(glm::vec3 containerBounds){
    int num_cells = 2 * (containerBounds.x * containerBounds.y + containerBounds.y * containerBounds.z + containerBounds.x * containerBounds.z);

    forces.resize(num_cells);
    vel_u.resize(num_cells + containerBounds.y * containerBounds.z);
    vel_v.resize(num_cells + containerBounds.x * containerBounds.z);
    vel_w.resize(num_cells + containerBounds.y * containerBounds.x);
    hash_grid.resize(num_cells);
    i_size = containerBounds.x;
    j_size = containerBounds.y;
}

void FluidSolver::initMACGrid(Particle &p){
    int index = findGridIndex(p.pos.x, p.pos.y, p.pos.z);
    p.gridPos = index;
    hash_grid.at(index).push_back(&p);
    std::cout << "index: " << index << std::endl;
}

void FluidSolver::particlesInit(){
    for(int i=0; i<MaxParticles; i++){
//        ParticlesContainer[i].life = -1.0f;
        ParticlesContainer[i].cameradistance = -1.0f;
        ParticlesContainer[i].size = 0.1f;
    }
}

void FluidSolver::genParticles(float particle_separation, float boundx, float boundy, float boundz){
    Particle p;
    for(float i = 0; i < boundx; i+= particle_separation){
        for(float j = 0; j < boundy ; j+= particle_separation){
            for(float k = 0; k <boundz ; k+= particle_separation){
                p.pos = glm::vec3(i, j, k);
                ParticlesContainer.push_back(p);
            }
        }
    }
}

int FluidSolver::findGridIndex(int x, int y, int z){
    return x + y * i_size + z * j_size * i_size;
}

//std::pair<Cell, int> FluidSolver::fillMACGrid(int cellsize, Particle p){
//        int pGridx = floor(p.pos.x / cellsize);
//        int pGridy = floor(p.pos.y / cellsize);
//        int pGridz = floor(p.pos.z / cellsize);

//        int cellindex = findGridIndex(pGridx, pGridy, pGridz);
//        p.gridPos = cellindex;
////        grid[cellindex].particles.push_back(&p);
//        Cell cell;
//        cell.particles.push_back(&p);
//        std::pair<Cell, int> mac_unit;
//        mac_unit.first = cell;
//        mac_unit.second = cellindex;
//        return mac_unit;
//}

void FluidSolver::calculateDensity(Particle &p){
    // for all neighboring particles: mass_i * kernel_function(p_i - x, h) / max_density

}

void FluidSolver::flipSolve(){
    for (int i = 0; i < vel_u.size(); ++i){
        vel_u_change[i] = vel_u[i] - vel_u_current[i];
    }
    for (int i = 0; i < vel_v.size(); ++i){
        vel_v_change[i] = vel_v[i] - vel_v_current[i];
    }
    for (int i = 0; i < vel_w.size(); ++i){
        vel_w_change[i] = vel_w[i] - vel_w_current[i];
    }

    particle_save = ParticlesContainer;



}

void FluidSolver::calculateGravityForces(Particle &p){
    p.speed += glm::vec3(0.f, -9.81f, 0.f);
}

void FluidSolver::storeParticleVelocityToGrid(){
    //for all the grid indices, calculate vel_u, vel_v, vel_w
    float u = 0.f, v = 0.f, w = 0.f, fract_part, dec_part;
    int index = -1;
    for(auto particle_list : hash_grid){
        while(particle_list.size() > 0){
            for(Particle* p : particle_list){
               //u component of velocity
               fract_part = modf(p->pos.x, &dec_part);
               u += fract_part * p->speed.x;
               //v component of velocity
               fract_part = modf(p->pos.y, &dec_part);
               v += fract_part * p->speed.y;
               //w component of velocity
               fract_part = modf(p->pos.z, &dec_part);
               w += fract_part * p->speed.z;
               index = p->gridPos;
            }
            vel_u[index] = u / particle_list.size();
            vel_v[index] = v / particle_list.size();
            vel_w[index] = w / particle_list.size();
        }
    }
}

void FluidSolver::calculateNewGridVelocities(){

}

void FluidSolver::ExtrapolateVelocity(){

}

void FluidSolver::setBoundaryVelocitiesToZero(const glm::vec3 containerBounds){
//    for(int x = ; x < containerBounds.x ; ++x){
//        vel_u[findGridIndex(x, )]
//    }
    for (int x = 0; x < containerBounds.x; ++x){
        vel_u[findGridIndex(x, 0, 0)] = 0.f;
        vel_u[findGridIndex(x, containerBounds.y - 1, 0)] = 0.f;
        vel_u[findGridIndex(x, containerBounds.y - 1, containerBounds.z - 1)] = 0.f;
        vel_u[findGridIndex(x, 0, containerBounds.z - 1)] = 0.f;
    }
    for (int y = 0; y < containerBounds.y; ++y){
        vel_v[findGridIndex(0, y, 0)] = 0.f;
        vel_v[findGridIndex(containerBounds.x - 1, y, 0)] = 0.f;
        vel_u[findGridIndex(containerBounds.x - 1, y, containerBounds.z - 1)] = 0.f;
        vel_u[findGridIndex(0, y, containerBounds.z - 1)] = 0.f;
    }
    for (int z = 0; z < containerBounds.z; ++z){
        vel_w[findGridIndex(0, 0, z)] = 0.f;
        vel_w[findGridIndex(containerBounds.x - 1, 0, z)] = 0.f;
        vel_w[findGridIndex(containerBounds.x - 1, containerBounds.y - 1, z)] = 0.f;
        vel_w[findGridIndex(0, containerBounds.y - 1, z)] = 0.f;
    }
}

void FluidSolver::projectPressure(){
    //check divergence
    float div_u, div_v, div_w;
    for (int i = 0; i < vel_u.size() - 1; ++i){
        div_u = vel_u[i + 1] - vel_u[i];
    }
    for (int i = 0; i < vel_v.size() - 1; ++i){
        div_v = vel_v[i + 1] - vel_v[i];
    }
    for (int i = 0; i < vel_w.size() - 1; ++i){
        div_w = vel_w[i + 1] - vel_w[i];
    }
    float divergence = div_u + div_v + div_w;
}

void FluidSolver::storeCurrentGridVelocities(){
    vel_u_current = vel_u;
    vel_v_current = vel_v;
    vel_w_current = vel_w;
}
