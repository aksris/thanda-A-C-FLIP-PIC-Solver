//
//  fluidSolver.cpp
//  Thanda
//#define DEBUG
#define VISCOSITY 0.f
#define EPSILON 0.00001f
#define SEED 8

#include "fluidSolver.hpp"
#include "../scene/scene.hpp"

Particle::Particle(){
    pos = glm::vec3(0.f, 0.f, 0.f);
    speed = glm::vec3(0.f, 0.f, 0.f);
    r = 0;
    g = 0;
    b = 220;
    a = 230;
    size = 0.1f;
    angle = 45.f;
    mass = 10.f;
    life = 1.f;
    cameradistance = 10.f;
    density = 1.f;
}

MACGrid::MACGrid(const ivec3& resolution, const vec3& containerBounds,float cellSize){
    vel_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
    vel_V = new MACGridDataY(ivec3(resolution.x, resolution.y+1, resolution.z), containerBounds, cellSize);
    vel_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z+1), containerBounds, cellSize);

    flip_vel_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
    flip_vel_V = new MACGridDataY(ivec3(resolution.x, resolution.y+1, resolution.z), containerBounds, cellSize);
    flip_vel_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z+1), containerBounds, cellSize);

    save_kernel_wt_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
    save_kernel_wt_V = new MACGridDataY(ivec3(resolution.x, resolution.y+1, resolution.z), containerBounds, cellSize);
    save_kernel_wt_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z+1), containerBounds, cellSize);

    P = new MACGridData(ivec3(resolution.x, resolution.y, resolution.z), containerBounds, cellSize);
}

MACGrid::~MACGrid(){
    delete vel_U;
    delete vel_V;
    delete vel_W;

    delete save_kernel_wt_U;
    delete save_kernel_wt_V;
    delete save_kernel_wt_W;

    delete P;
}

void MACGrid::initialize(){
    vel_U->MACGridDataInitialize();
    vel_V->MACGridDataInitialize();
    vel_W->MACGridDataInitialize();
    P->MACGridDataInitialize();
    save_kernel_wt_U->MACGridDataInitialize();
    save_kernel_wt_V->MACGridDataInitialize();
    save_kernel_wt_W->MACGridDataInitialize();
}

MACGrid& MACGrid::operator =(const MACGrid& val){
    if (&val == this)
    {
        return *this;
    }
    vel_U->data = val.vel_U->data;
    vel_U->mData = val.vel_U->mData;

    vel_V->data = val.vel_V->data;
    vel_V->mData = val.vel_V->mData;

    vel_W->data = val.vel_W->data;
    vel_W->mData = val.vel_W->mData;

    P->data = val.P->data;
    P->mData = val.P->mData;

    save_kernel_wt_U = val.save_kernel_wt_U;
    save_kernel_wt_V = val.save_kernel_wt_V;
    save_kernel_wt_W = val.save_kernel_wt_W;
    return *this;
}

FluidSolver::FluidSolver(const ivec3& resolution, const vec3& containerBounds){
    LastUsedParticle = 0;
    MaxParticles = 200000;
    delta = 0.1f;

    this->resolution = resolution;
    this->containerBounds = containerBounds;
    this->cellSize = containerBounds.x/resolution.x;

    grid = new MACGrid(this->resolution, this->containerBounds, this->cellSize);
//    tmp = new MACGrid(this->resolution, this->containerBounds, this->cellSize);
}

FluidSolver::~FluidSolver(){
    delete grid;
//    delete tmp;
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

void FluidSolver::constructMACGrid(const Scene& scene){
    resolution = scene.resolution;
    containerBounds = scene.containerBounds;

    num_cells = (resolution.x * resolution.y * resolution.z);
//    grid->initialize();
}

void FluidSolver::initMACGrid(Particle &p){

}

void FluidSolver::calculateGravityForces(Particle &p, float delta){
    p.speed += glm::vec3(0.f, -9.81f , 0.f) * delta;
    p.pos += p.speed * delta;
}

void FluidSolver::particlesInit(){
    for(int i=0; i<MaxParticles; i++){
        ParticlesContainer[i].cameradistance = -1.0f;
        ParticlesContainer[i].size = 0.1f;
    }
}
float rand(int LO, int HI){
    return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
}

void FluidSolver::genParticles(float particle_separation, float boundx, float boundy, float boundz){
    Particle p;
    int iter;
    for(int i = 1; i < (int)boundx; i++){
        for(int j = 1; j < (int)boundy; j++){
            for(int k = 1; k < (int)boundz; k++){
                iter=0;
                while(iter < SEED){
//                    P->pos = vec3(1.8,1.8,1.6);
                    p.pos = vec3(rand(i,i+1), rand(j,j+1), rand(k, k+1)) * this->cellSize;
                    p.speed = vec3(0.f, 0.f, 0.f);
                    p.gridIdx = vec3(i,j,k);
                    ParticlesContainer.push_back(p);
                    iter++;
                }
            }
        }
    }
}

void FluidSolver::calculateDensity(Particle &p){
    // for all neighboring particles: mass_i * kernel_function(p_i - x, h) / max_density

}

void FluidSolver::FlipSolve(){
    int x = grid->P->resolution.x;
    int y = grid->P->resolution.y;
    int z = grid->P->resolution.z;
    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                grid->flip_vel_U->setCell(i, j, k, (*grid->vel_U)(i, j, k) - (*grid->flip_vel_U)(i, j, k));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                grid->flip_vel_V->setCell(i, j, k, (*grid->vel_V)(i, j, k) - (*grid->flip_vel_V)(i, j, k)) ;
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z + 1; ++k){
                grid->flip_vel_W->setCell(i, j, k, (*grid->vel_W)(i, j, k) - (*grid->flip_vel_W)(i, j, k));
            }
        }
    }

    particle_save = ParticlesContainer;
    //for every particle, set the change in velocity + current particle velocity
    //interpolate
    for(int i = 0; i < particle_save.size(); i++){
        particle_save.at(i).speed.x = grid->flip_vel_U->interpolate(particle_save.at(i).pos);
        particle_save.at(i).speed.y = grid->flip_vel_V->interpolate(particle_save.at(i).pos);
        particle_save.at(i).speed.z = grid->flip_vel_W->interpolate(particle_save.at(i).pos);
    }

}

void FluidSolver::PicSolve(){
    particle_save_pic = ParticlesContainer;
    //for every particle, set the new grid velocity
    //interpolate
    for(int i = 0; i < particle_save_pic.size(); i++){
        particle_save_pic.at(i).speed.x = grid->vel_U->interpolate(particle_save_pic.at(i).pos);
        particle_save_pic.at(i).speed.y = grid->vel_V->interpolate(particle_save_pic.at(i).pos);
        particle_save_pic.at(i).speed.z = grid->vel_W->interpolate(particle_save_pic.at(i).pos);
    }

}

void FluidSolver::naiveNeighborSearch(Particle *p, std::vector<Particle> &neighbors) {
    for (int i = 0; i < ParticlesContainer.size(); i++) {
        float dist = glm::length(p->pos - ParticlesContainer[i].pos);
        if (dist < 1.4f) {
            neighbors.push_back(ParticlesContainer[i]);
        }
    }
}

float Sqrlength(const glm::vec3& p0, const glm::vec3& p1){
    float a = p0.x - p1.x;
    float b = p0.y - p1.y;
    float c = p0.z - p1.z;
    return a*a + b*b + c*c;
}

float Smooth(const float& r2, const float& h) {
    return glm::max(1.0f-r2/(h*h), 0.0f);
}

float Sharpen(const float& r2, const float& h) {
    return glm::max(h*h/glm::max(r2,(float)1.0e-5) - 1.0f, 0.0f);
}

float h(float r){
    if(r >= 0.f && r <= 1.f)
        return 1.f - r;
    else if(-1.f <= r && r < 0.f)
        return 1.f + r;
    else
        return 0.f;
}

float StiffKernel(const vec3& r, const float& cell_width){
//    return h(r.x/cell_width)*h(r.y/cell_width)*h(r.z/cell_width);
    float r2 = length2(r);
    return max(1.f - r2 / (cell_width*cell_width), EPSILON);
}

void FluidSolver::CalculateGravityToCell(float delta){
    vec3 speed;
    for (int i = 0; i < resolution.x; i++){
        for (int j = 0; j < resolution.y; j++){
            for (int k = 0; k < resolution.x; k++){
                speed = glm::vec3(0.f, -9.8f , 0.f) * delta;
                grid->vel_V->setCellAdd(i,j,k, speed.y);
            }
        }
    }
}

void FluidSolver::storeParticleVelocityToGrid(){
    //for all the grid indices, calculate vel_u, vel_v, vel_w
    vec3 index, pos, r; int x,y,z;
    float h = this->cellSize, weight = 0.f;
    float fract_partx,fract_party,fract_partz;
    
    for(int i = 0; i < ParticlesContainer.size(); ++i){
        
        const Particle& par = ParticlesContainer.at(i);
        
        vec3 w_pos = par.pos;
        index = par.gridIdx;
        x = index.x, y = index.y, z = index.z;
        weight = 0.f;
        
        grid->P->setCellMark(x, y, z, FLUID );

        //################## X direction #############################################
        pos = grid->vel_U->worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_U->setCellAdd(x, y, z, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_U->setCellAdd(x+1, y, z, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_U->setCellAdd(x+1, y+1, z, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_U->setCellAdd(x, y+1, z, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid->vel_U->setCellAdd(x, y+1, z+1, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_U->setCellAdd(x, y, z+1, par.speed.x * weight );
        grid->save_kernel_wt_U->setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_U->setCellAdd(x+1, y, z+1, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid->vel_U->setCellAdd(x+1, y+1, z+1, par.speed.x * weight);
        grid->save_kernel_wt_U->setCellAdd(x+1, y+1, z+1, weight);

        //######################## y direction splatting #####################################
        pos = grid->vel_V->worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        float speedY = par.speed.y;
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_V->setCellAdd(x, y, z, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_V->setCellAdd(x+1, y, z, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_V->setCellAdd(x+1, y+1, z, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_V->setCellAdd(x, y+1, z, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid->vel_V->setCellAdd(x, y+1, z+1, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_V->setCellAdd(x, y, z+1, speedY * weight );
        grid->save_kernel_wt_V->setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_V->setCellAdd(x+1, y, z+1, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid->vel_V->setCellAdd(x+1, y+1, z+1, speedY * weight);
        grid->save_kernel_wt_V->setCellAdd(x+1, y+1, z+1, weight);

        //############### z direction #######################################################
        pos = grid->vel_W->worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_W->setCellAdd(x, y, z, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid->vel_W->setCellAdd(x+1, y, z, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_W->setCellAdd(x+1, y+1, z, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid->vel_W->setCellAdd(x, y+1, z, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid->vel_W->setCellAdd(x, y+1, z+1, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_W->setCellAdd(x, y, z+1, par.speed.z * weight );
        grid->save_kernel_wt_W->setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid->vel_W->setCellAdd(x+1, y, z+1, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid->vel_W->setCellAdd(x+1, y+1, z+1, par.speed.z * weight);
        grid->save_kernel_wt_W->setCellAdd(x+1, y+1, z+1, weight);

    }

    x = grid->P->resolution.x;
    y = grid->P->resolution.y;
    z = grid->P->resolution.z;
    for(int i = 0; i < x+1; ++i){
        for(int j = 0; j < y ; ++j){
            for(int k = 0; k < z; ++k){
                grid->vel_U->setCell(i, j, k, (*grid->vel_U)(i,j,k) / max((*grid->save_kernel_wt_U)(i,j,k), EPSILON));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                grid->vel_V->setCell(i, j, k, (*grid->vel_V)(i,j,k) / max((*grid->save_kernel_wt_V)(i,j,k), EPSILON));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y ; ++j){
            for(int k = 0; k < z+1; ++k){
                grid->vel_W->setCell(i, j, k, (*grid->vel_W)(i,j,k) / max((*grid->save_kernel_wt_W)(i,j,k), EPSILON));
            }
        }
    }
}

//void FluidSolver::RepositionParticles(){
//    for(auto i : reposition_map){
//        Particle& p = ParticlesContainer.at(i.first);
//#ifdef DEBUG
//        if(p.gridIdx.x < 0 || p.gridIdx.x > resolution.x - 1 ||
//                p.gridIdx.y < 0 || p.gridIdx.y > resolution.y - 1 ||
//                p.gridIdx.z < 0 || p.gridIdx.z > resolution.z - 1){
//            int a = 1;
//        }
//#endif
//
//        // change particle info
//        p.pos = vec3(i.second) * cellSize + vec3(rand(0,1));
//        p.gridIdx = vec3(i.second);
//        
//        // mark new grid cell where particle was repositioned to
//        this->grid->P->setCellMark(i.second.x, i.second.y, i.second.z, FLUID);
//    }
//    
//    reposition_map.clear();
//}

void FluidSolver::initializeMarkerGrid(){
    int x = resolution.x;
    int y = resolution.y;
    int z = resolution.z;

    //neighborhood of 0, based on distance

    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(i == 0 || i == x - 1 || j == 0 || j == y - 1 || k == 0 || k == z - 1)
                {
                    grid->P->setCellMark(i, j, k, SOLID );
                }
                else{
                    grid->P->setCellMark(i,j,k, AIR );
                }
            }
        }
    }

    // mark particles to be repositioned and the appropriate neighbor cell to reposition into
//    for(int p = 0; p< ParticlesContainer.size(); ++p){
//        int i = ParticlesContainer.at(p).gridIdx.x, j = ParticlesContainer.at(p).gridIdx.y, k = ParticlesContainer.at(p).gridIdx.z;
//        vec3 list[6] = { vec3(i-1,j,k), vec3(i+1,j,k),
//                         vec3(i,j-1,k), vec3(i,j+1,k),
//                         vec3(i,j,k-1), vec3(i,j,k+1)
//        };
//        if(grid->P->getCellMark(i,j,k) == SOLID){
//            for(int n = 0; n < 6; ++n){
//                if(list[n][0] > 0 && list[n][0]< x-1 && list[n][1] > 0 &&
//                        list[n][1] < y-1 && list[n][2] > 0 && list[n][2] < z-1){
//                    if(grid->P->getCellMark(list[n]) != SOLID){
//                        reposition_map[p] = ivec3(list[n][0], list[n][1], list[n][2]);
//                        break;
//                    }
//                }
//            }
//        }
//        else{
//            grid->P->setCellMark(i,j,k, FLUID );
//        }
//    }
}


void FluidSolver::insertCoefficient(int id, int i, int j, int k, double w, std::vector<Eigen::Triplet<double>> &coeffs, int n){
    n = 5;
    int id1 = k * resolution.x * resolution.y + j * resolution.x + i;
    coeffs.push_back(Eigen::Triplet<double>(id, id1, w));
}

void FluidSolver::buildMatrixA(std::vector<Eigen::Triplet<double>>& coefficients, long n){
    float density = 1000.f;
    for(int i = 0; i < resolution.x; ++i){
        for(int j = 0; j < resolution.y; ++j){
            for(int k = 0; k < resolution.z; ++k){
                int id = k * resolution.x * resolution.y + j * resolution.x + i; //id for matrix
                float scale = delta / (density * cellSize * cellSize);
                float Adiag = 0.f;
                if(grid->P->getCellMark(i,j,k) == FLUID){
                    //x
                    if(i > 0 && grid->P->getCellMark(i-1,j,k) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i-1,j,k, -scale, coefficients, n);
                    }
                    if(i < n && grid->P->getCellMark(i+1,j,k) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i+1,j,k, -scale, coefficients, n);
                    }
                    if(i > 0 && grid->P->getCellMark(i-1,j,k) == AIR){
                        Adiag += scale;
                    }
                    if(i < n && grid->P->getCellMark(i+1,j,k) == AIR){
                        Adiag += scale;
                    }

                    //y
                    if(j > 0 && grid->P->getCellMark(i,j-1,k) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i,j-1,k, -scale, coefficients, n);
                    }
                    if(j < n && grid->P->getCellMark(i,j+1,k) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i,j+1,k, -scale, coefficients, n);
                    }
                    if(j > 0 && grid->P->getCellMark(i,j-1,k) == AIR){
                        Adiag += scale;
                    }
                    if(j < n && grid->P->getCellMark(i,j+1,k) == AIR){
                        Adiag += scale;
                    }

                    //z
                    if(k > 0 && grid->P->getCellMark(i,j,k-1) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i,j,k-1, -scale, coefficients, n);
                    }
                    if(k < n && grid->P->getCellMark(i,j,k+1) == FLUID){
                        Adiag += scale;
                        insertCoefficient(id, i,j,k+1, -scale, coefficients, n);
                    }
                    if(k > 0 && grid->P->getCellMark(i,j,k-1) == AIR){
                        Adiag += scale;
                    }
                    if(k < n && grid->P->getCellMark(i,j,k+1) == AIR){
                        Adiag += scale;
                    }
                }
                insertCoefficient(id, i,j,k, Adiag, coefficients, n);
            }
        }
    }
}

void FluidSolver::buildDivergences(Eigen::VectorXd& rhs){
    float scale = 1/this->cellSize;
    double divergence = 0.f;
    for(int i = 0; i < resolution.x; ++i){
        for(int j = 0; j < resolution.y; ++j){
            for(int k = 0; k < resolution.z; ++k){
                if(grid->P->getCellMark(i,j,k) == FLUID){
                    int id = k * resolution.x * resolution.y + j * resolution.x + i;
                    divergence =scale * (
                                        (*grid->vel_U)(i+1, j, k) -
                                        (*grid->vel_U)(i, j, k)
                                         +
                                        (*grid->vel_V)(i, j+1, k) -
                                        (*grid->vel_V)(i, j, k)
                                         +
                                        (*grid->vel_W)(i, j, k+1) -
                                        (*grid->vel_W)(i, j, k)
                                        );
#ifdef DEBUG
                    if(id == 31){
                        std::cout << (*grid->vel_U)(i+1, j, k) << ", " <<
                                     (*grid->vel_U)(i, j, k)  << ", " <<
                                     (*grid->vel_V)(i, j+1, k) << ", " <<
                                     (*grid->vel_V)(i, j, k)    << ", " <<
                                     (*grid->vel_W)(i, j, k+1) << ", " <<
                                     (*grid->vel_W)(i, j, k) << std::endl;
                    }
#endif
                    rhs[id] = -divergence;
                    
                    if (grid->P->getCellMark(i-1, j, k)==SOLID) rhs[id] -= scale * (*grid->vel_U)(i, j, k);
                    if (grid->P->getCellMark(i+1, j, k)==SOLID) rhs[id] += scale * (*grid->vel_U)(i+1, j, k);
                    
                    if (grid->P->getCellMark(i, j-1, k)==SOLID) rhs[id] -= scale * (*grid->vel_V)(i, j, k);
                    if (grid->P->getCellMark(i, j+1, k)==SOLID) rhs[id] += scale * (*grid->vel_V)(i, j+1, k);

                    if (grid->P->getCellMark(i, j, k-1)==SOLID) rhs[id] -= scale * (*grid->vel_W)(i, j, k);
                    if (grid->P->getCellMark(i, j, k+1)==SOLID) rhs[id] += scale * (*grid->vel_W)(i, j, k+1);
                }
            }
        }
    }
}

void FluidSolver::fillPressureGrid(Eigen::VectorXd x){
    int iter = 0;
    for(int i = 0; i < resolution.x; ++i){
        for(int j = 0; j < resolution.y; ++j){
            for(int k = 0; k < resolution.z; ++k){
                if (std::isnan(x[iter]))
                    std::cout<<"ERROR with PCG"<<std::endl;
//                grid->P->setCell(i,j,k, max((float)x[iter++], 0.f));
//                grid->P->setCell(i,j,k, (float)x[iter++]);
                grid->P->setCell(i,j,k, max((float)x[iter++], 0.f));
//#ifdef DEBUG
                if(grid->P->getCellMark(i,j,k) == FLUID){
                    float x = (*grid->P)(i,j,k);
                    std::cout << "Pressure at "<< i << " "<< j << " " << k << ": "<<x<< std::endl;
                }
//#endif
            }
        }
    }
}

void FluidSolver::ProjectPressure(){
    int x = resolution.x, y = resolution.y, z = resolution.z;
    long m = x * y * z;


    Eigen::VectorXd p(m);
    p.setZero();

    Eigen::VectorXd rhs(m);
    rhs.setZero();
    buildDivergences(rhs);
    
    std::vector<Eigen::Triplet<double>> coefficients;
    buildMatrixA(coefficients, m);
    Eigen::SparseMatrix<double> A(m,m);
    A.setZero();
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    
    
    // solve
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;  // performs a Cholesky factorization of A
//    cg.compute(A);
//    p = cg.solve(rhs);

    Eigen::IncompleteCholesky<double, Eigen::Lower|Eigen::Upper> pcg(A); //calls pcg.compute internally
    p = pcg.solve(rhs);
    
    // fill eigen vector into grid.P
    fillPressureGrid(p);

//    std::cout << u(31) << std::endl;

    SubtractPressureGradient();
//    buildDivergences(rhs);
//    std::cout << u(31) << std::endl;
//    int a =1;
}

void FluidSolver::SubtractPressureGradient(){
    float density = 1000.f;
    float dx = cellSize;
    float scale = delta / (density * dx);

    int x = resolution.x;
    int y = resolution.y;
    int z = resolution.z;

    //loop over i,j,k
    //update vel_U
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(grid->P->getCellMark(i-1,j,k) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
                    if(grid->P->getCellMark(i-1,j,k) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid->vel_U->setCell(i,j,k,0.f);
                    }
                    else{
                        grid->vel_U->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i-1,j,k))));
                    }
                }
                if(grid->P->getCellMark(i,j-1,k) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
                    if(grid->P->getCellMark(i,j-1,k) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid->vel_V->setCell(i,j,k,0.f);
                    }
                    else{
                        grid->vel_V->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i,j-1,k))));
                    }
                }
                if(grid->P->getCellMark(i,j,k-1) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
                    if(grid->P->getCellMark(i,j,k-1) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid->vel_W->setCell(i,j,k,0.f);
                    }
                    else{
                        grid->vel_W->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i,j,k-1))));
                    }
                }
            }
        }
    }
    
//    for(int i = 0; i < x+1; ++i){
//        for(int j = 0; j < y; ++j){
//            for(int k = 0; k < z; ++k){
//                if(grid->P->getCellMark(i-1,j,k) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
//                    if(grid->P->getCellMark(i-1,j,k) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
//                        //set cell value as usolid at i,j,k
//                        grid->vel_U->setCell(i,j,k,0.f);
//                    }
//                    else{
//                        grid->vel_U->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i-1,j,k))));
//                    }
//                }
//            }
//        }
//    }
//    //update vel_V
//    for(int i = 0; i < x; ++i){
//        for(int j = 0; j < y+1; ++j){
//            for(int k = 0; k < z; ++k){
//                if(grid->P->getCellMark(i,j-1,k) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
//                    if(grid->P->getCellMark(i,j-1,k) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
//                        //set cell value as usolid at i,j,k
//                        grid->vel_V->setCell(i,j,k,0.f);
//                    }
//                    else{
//                        grid->vel_V->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i,j-1,k))));
//                    }
//                }
//            }
//        }
//    }
//    //update vel_W
//    for(int i = 0; i < x; ++i){
//        for(int j = 0; j < y; ++j){
//            for(int k = 0; k < z+1; ++k){
//                if(grid->P->getCellMark(i,j,k-1) == FLUID || grid->P->getCellMark(i,j,k) == FLUID){
//                    if(grid->P->getCellMark(i,j,k-1) == SOLID || grid->P->getCellMark(i,j,k) == SOLID){
//                        //set cell value as usolid at i,j,k
//                        grid->vel_W->setCell(i,j,k,0.f);
//                    }
//                    else{
//                        grid->vel_W->setCellAdd(i,j,k, -(scale * ((*grid->P)(i,j,k) - (*grid->P)(i,j,k-1))));
//                    }
//                }
//            }
//        }
//    }
}

vec3 EulerStep (const vec3 pos, const vec3 speed, float time_step){
    return pos + speed * time_step;
}

vec3 FluidSolver::integratePos(const vec3& pos, const vec3& speed, const float& time_step, bool RK2){
    vec3 new_pos(0.f);
    if(RK2){
        //RK2 integration
        vec3 k1 = speed * (pos) * time_step / 2.f;
        vec3 k2 = speed * (pos + k1) * time_step;
//        new_pos = pos + k2;
        new_pos = pos + speed * time_step;
    }
    else{
        //RK4
        vec3 k1 = speed * (pos) * time_step / 2.f;
        vec3 k2 = speed * (pos + k1) * time_step;
        vec3 k3 = speed * (pos + k2) * time_step;
        vec3 k4 = speed * (pos + k3);
        new_pos = pos + (0.1666666f) * (k1 + 2.f * k2 + 2.f * k3 + k4);
    }
    
    new_pos.x = glm::clamp(new_pos.x, cellSize*1.f, cellSize*(resolution.x-1));
    new_pos.y = glm::clamp(new_pos.y, cellSize*1.f, cellSize*(resolution.y-1));
    new_pos.z = glm::clamp(new_pos.z, cellSize*1.f, cellSize*(resolution.z-1));
    
    return new_pos;
}

void FluidSolver::ExtrapolateVelocity(){
    int x = grid->P->resolution.x;
    int y = grid->P->resolution.y;
    int z = grid->P->resolution.z;

    MACGridData extrapolate_grid_U(resolution, containerBounds, cellSize);
    MACGridData extrapolate_grid_V(resolution, containerBounds, cellSize);
    MACGridData extrapolate_grid_W(resolution, containerBounds, cellSize);
    extrapolate_grid_U.MACGridDataInitialize();
    extrapolate_grid_V.MACGridDataInitialize();
    extrapolate_grid_W.MACGridDataInitialize();

    //making three grid structures for all axes, so that we store the (fluid or near fluid) and (solid/air or near solid/(solid or air)) mark
    //this is to extrapolate in to the right cells so that tri lerp does not lerp with an empty cell

    //setting a temporary marker grid to get extrapolated velocities

    //x direction
    for(int i = 0; i < x+1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                //marking fluid or near fluid cells as FLUID
                if((i < x && grid->P->getCellMark(i,j,k) == FLUID) || (i > 0 && grid->P->getCellMark(i-1,j,k) == FLUID)){
                    extrapolate_grid_U.setCellMark(i,j,k, FLUID );
                }
                //marking solid for extrapolation
                else if((i < x && (grid->P->getCellMark(i,j,k) == SOLID)) &&
                        (i > 0 && (grid->P->getCellMark(i-1,j,k) == SOLID )) //check if me solid and near solid
                        || (i < x && (grid->P->getCellMark(i,j,k) == AIR)) &&
//                        (i > 0 && (grid->P->getCellMark(i-1,j,k) == SOLID || grid->P->getCellMark(i-1,j,k) == AIR)) // not bridson
                        (i > 0 && (grid->P->getCellMark(i-1,j,k) == AIR)) // bridson
                        ) //me air, and near (solid or air)
                {
                    extrapolate_grid_U.setCellMark(i,j,k,SOLID );
                }
            }
        }
    }
    //y direction
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y+1; ++j){
            for(int k = 0; k < z; ++k){
                //marking fluid or near fluid cells as FLUID
                if((j < y && grid->P->getCellMark(i,j,k) == FLUID) || (j > 0 && grid->P->getCellMark(i,j-1,k) == FLUID)){
                    extrapolate_grid_V.setCellMark(i,j,k, FLUID);
                }
                //marking solid for extrapolation
                else if((j < y && (grid->P->getCellMark(i,j,k) == SOLID)) &&
                        (j > 0 && (grid->P->getCellMark(i,j-1,k) == SOLID )) //check if me solid and near solid
                        || (j < y && (grid->P->getCellMark(i,j,k) == AIR)) &&
//                        (j > 0 && (grid->P->getCellMark(i,j-1,k) == SOLID || grid->P->getCellMark(i,j-1,k) == AIR))) //me air, and near (solid or air)
                    (j > 0 && (grid->P->getCellMark(i,j-1,k) == AIR))) //me air, and near (solid or air)

                {
                    extrapolate_grid_V.setCellMark(i,j,k,SOLID );
                }
            }
        }
    }
    //z direction
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z+1; ++k){
                //marking fluid or near fluid cells as FLUID
                if((k < z && grid->P->getCellMark(i,j,k) == FLUID) || (k > 0 && grid->P->getCellMark(i,j,k-1) == FLUID)){
                    extrapolate_grid_W.setCellMark(i,j,k, FLUID );
                }
                //marking solid for extrapolation
                else if((k < z && (grid->P->getCellMark(i,j,k) == SOLID)) &&
                        (k > 0 && (grid->P->getCellMark(i,j,k-1) == SOLID )) //check if me solid and near solid
                        || (k < z && (grid->P->getCellMark(i,j,k) == AIR)) &&
//                        (k > 0 && (grid->P->getCellMark(i,j,k-1) == SOLID || grid->P->getCellMark(i,j,k-1) == AIR))) //me air, and near (solid or air)
                    (k > 0 && (grid->P->getCellMark(i,j,k-1) == AIR))) //me air, and near (solid or air)
                {
                    extrapolate_grid_W.setCellMark(i,j,k,SOLID);
                }
            }
        }
    }

    //neighborhood of 0, based on distance
    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z + 1; ++k){
                //x faces
                if((j < y && k < z) && extrapolate_grid_U.getCellMark(i, j, k) == SOLID){
                    unsigned int wsum = 0;
                    float sum = 0.0f;
                    glm::vec3 q[6] = { glm::vec3(i-1,j,k), glm::vec3(i+1,j,k),
                                       glm::vec3(i,j-1,k), glm::vec3(i,j+1,k),
                                       glm::vec3(i,j,k-1), glm::vec3(i,j,k+1) };
                    for(unsigned int qk = 0; qk < 6; ++qk){
                        if(q[qk][0] >= 0 && q[qk][0]< x + 1 && q[qk][1] >= 0 &&
                                q[qk][1] < y && q[qk][2] >= 0 && q[qk][2] < z ) {

                            if(extrapolate_grid_U.getCellMark(q[qk][0], q[qk][1], q[qk][2]) == FLUID){
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid->vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;
                                sum += (*grid->vel_U)(q[qk][0],q[qk][1],q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        grid->vel_U->setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid->vel_V(i,j+1,k) << std::endl;
                }
                //y faces
                if((i < x && k < z) && extrapolate_grid_V.getCellMark(i, j, k) == SOLID){
                    unsigned int wsum = 0;
                    float sum = 0.0f;
                    glm::vec3 q[6] = { glm::vec3(i-1,j,k), glm::vec3(i+1,j,k),
                                       glm::vec3(i,j-1,k), glm::vec3(i,j+1,k),
                                       glm::vec3(i,j,k-1), glm::vec3(i,j,k+1) };
                    for(unsigned int qk = 0; qk < 6; ++qk){
                        if(q[qk][0] >= 0 && q[qk][0]< x && q[qk][1] >= 0 &&
                                q[qk][1] < y+1 && q[qk][2] >= 0 && q[qk][2] < z ) {

                            if(extrapolate_grid_V.getCellMark(q[qk][0], q[qk][1], q[qk][2]) == FLUID){
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid->vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;

                                sum += (*grid->vel_V)(q[qk][0],q[qk][1],
                                        q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " : "<< sum/wsum<<std::endl;
                        grid->vel_V->setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid->vel_V(i,j+1,k) << std::endl;
                }
                //z faces
                if((j < y && i < x) && extrapolate_grid_W.getCellMark(i, j, k) == SOLID){
                    unsigned int wsum = 0;
                    float sum = 0.0f;
                    glm::vec3 q[6] = { glm::vec3(i-1,j,k), glm::vec3(i+1,j,k),
                                       glm::vec3(i,j-1,k), glm::vec3(i,j+1,k),
                                       glm::vec3(i,j,k-1), glm::vec3(i,j,k+1) };
                    for(unsigned int qk = 0; qk < 6; ++qk){
                        if(q[qk][0] >= 0 && q[qk][0]< x && q[qk][1] >= 0 &&
                                q[qk][1] < y && q[qk][2] >= 0 && q[qk][2] < z+1 ) {

                            if(extrapolate_grid_W.getCellMark(q[qk][0], q[qk][1], q[qk][2]) == FLUID){
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid->vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;

                                sum += (*grid->vel_W)(q[qk][0],q[qk][1],
                                        q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        grid->vel_W->setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid->vel_V(i,j+1,k) << std::endl;
                }
            }

        }
    }
    
    // on the boundary there's only 1 FLUID neighbor (+1 in the dimension) - more efficient to set just that neighbor
}


void FluidSolver::setBoundaryVelocitiesToZero(){

    int x = grid->P->resolution.x;
    int y = grid->P->resolution.y;
    int z = grid->P->resolution.z;
    //neighborhood of 0, based on distance

    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(i == 0 || i == x){
                    grid->vel_U->setCell(i, j, k, 0.f);
                }
                if (i > 0 && i < x && ((grid->P->getCellMark(i,j,k) == SOLID && grid->P->getCellMark(i-1,j,k) == FLUID) ||
                                       (grid->P->getCellMark(i,j,k) == FLUID && grid->P->getCellMark(i-1,j,k) == SOLID))){
                    grid->vel_U->setCell(i, j, k, 0.f);
                }
            }
        }
    }


    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                if(j == 0 || j == y){
                    grid->vel_V->setCell(i, j, k, 0.f);
                }
                if (j > 0 && j < y && ((grid->P->getCellMark(i,j,k) == SOLID && grid->P->getCellMark(i,j-1,k) == FLUID) ||
                                       (grid->P->getCellMark(i,j,k) == FLUID && grid->P->getCellMark(i,j-1,k) == SOLID))){
                    grid->vel_V->setCell(i, j, k, 0.f);
                }
            }
        }
    }

    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z + 1; ++k){
                if(k == 0 || k == z){
                    grid->vel_W->setCell(i, j, k, 0.f);
                }
                if (k > 0 && k < z && ((grid->P->getCellMark(i,j,k) == SOLID && grid->P->getCellMark(i,j,k-1) == FLUID) ||
                                       (grid->P->getCellMark(i,j,k) == FLUID && grid->P->getCellMark(i,j,k-1) == SOLID))){
                    grid->vel_W->setCell(i, j, k, 0.f);
                }
            }
        }
    }
}

void FluidSolver::storeCurrentGridVelocities(){
//    tmp=grid;
    *this->grid->flip_vel_U = *this->grid->vel_U;
    *this->grid->flip_vel_V = *this->grid->vel_V;
    *this->grid->flip_vel_W = *this->grid->vel_W;
    
}

void FluidSolver::clearGrid(){
    std::fill(grid->save_kernel_wt_U->data.begin(), grid->save_kernel_wt_U->data.end(), 0.f);
    std::fill(grid->save_kernel_wt_V->data.begin(), grid->save_kernel_wt_V->data.end(), 0.f);
    std::fill(grid->save_kernel_wt_W->data.begin(), grid->save_kernel_wt_W->data.end(), 0.f);
    std::fill(grid->vel_U->data.begin(), grid->vel_U->data.end(), 0.f);
    std::fill(grid->vel_V->data.begin(), grid->vel_V->data.end(), 0.f);
    std::fill(grid->vel_W->data.begin(), grid->vel_W->data.end(), 0.f);
    std::fill(grid->P->data.begin(), grid->P->data.end(), 0.f);
}

void FluidSolver::step(){

    this->initializeMarkerGrid();
//    this->RepositionParticles();

//#define DEBUG
    
    // Step 3 - Store Particle Velocity at current time step to MACGrid
    int i = this->ParticlesContainer.at(0).gridIdx[0];
    int j = this->ParticlesContainer.at(0).gridIdx[1];
    int k = this->ParticlesContainer.at(0).gridIdx[2];
    
    this->storeParticleVelocityToGrid();
#ifdef DEBUG
    std::cout << "[thanda] storeParticleVelocityToGrid Iteration" << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j, k) << " at idx " << i << " "<< j << " "<< k << " : "<< (*grid->vel_V)(i,j,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< (*grid->vel_V)(i,j+1,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i+1, j, k) << " at idx " << i+1 << " "<< j << " "<< k << " : "<< (*grid->vel_V)(i+1,j,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i+1, j+1, k) << " at idx " << i+1 << " "<< j+1 << " "<< k << " : "<< (*grid->vel_V)(i+1,j+1,k) << std::endl; // is this valid

    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j, k+1) << " at idx " << i << " "<< j << " "<< k+1 << " : "<< (*grid->vel_V)(i,j,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i, j+1, k+1) << " at idx " << i << " "<< j+1 << " "<< k+1 << " : "<< (*grid->vel_V)(i,j+1,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i+1, j, k+1) << " at idx " << i+1 << " "<< j << " "<< k+1 << " : "<< (*grid->vel_V)(i+1,j,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid->P->getCellMark(i+1, j+1, k+1) << " at idx " << i+1 << " "<< j+1 << " "<< k+1 << " : "<< (*grid->vel_V)(i+1,j+1,k+1) << std::endl; // is this valid

#endif
    
    
    // Step 4 - Store a temporary copy of Grid Velocities for FLIP
    this->storeCurrentGridVelocities();

    // Step 5 - Add Body Forces like Gravity to MACGrid
    this->CalculateGravityToCell(delta);
    
    this->ExtrapolateVelocity();
    this->setBoundaryVelocitiesToZero();
    
    //pressure solve
    this->ProjectPressure();
    
    this->ExtrapolateVelocity();

    // Step  - Calculate new flip & pic velocities for each particle
//    this->FlipSolve();
    this->PicSolve();

    // Step - Lerp(FLIPVelocity, PICVelocity, 0.95)
    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
        this->ParticlesContainer.at(i).speed = (1.f - VISCOSITY) * this->particle_save_pic.at(i).speed;
//        + VISCOSITY * this->particle_save.at(i).speed;
        this->ParticlesContainer.at(i).pos = this->integratePos(this->ParticlesContainer.at(i).pos,
                                                                this->ParticlesContainer.at(i).speed, delta, true );
    }

    

//     Step - Collision Response
//    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
//        
//        Particle& particle = this->ParticlesContainer.at(i); // careful - by reference;
//        
//        float dampingFactor = 1.f;
//        
////            if(particle.pos.x  < EPSILON || particle.pos.x  > upperBounds.x - EPSILON)
////            {
////                currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(-dampingFactor,1,1));
////                currParticle.setPredictedPosition(currParticle.getPosition() + timeStep * currParticle.getVelocity());
////            }
//            
//        if (particle.pos.y < EPSILON){
//            if(particle.speed.y < EPSILON){
//                particle.speed *= (vec3(1.f, -dampingFactor, 1.f));
//            }
//            particle.pos.y = EPSILON ;
//            particle.pos = this->integratePos(particle.pos, particle.speed, delta, true);
//        }
//
//        else if(particle.pos.y  > containerBounds.y - EPSILON){// does not work
//            if(particle.speed.y > EPSILON){
//                particle.speed *= (vec3(1.f, -dampingFactor, 1.f));
//            }
//            particle.pos.y = containerBounds.y - EPSILON ;
//            particle.pos = this->integratePos(particle.pos, particle.speed, delta, true);
//        }
//    
//            
////            if(particle.pos.z  < lowerBounds.z + EPSILON || particle.pos.z  > upperBounds.z - EPSILON)
////            {
////                currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(1,1,-dampingFactor));
////                currParticle.setPredictedPosition(currParticle.getPosition() + timeStep * currParticle.getVelocity());
////            }
//            
//            
//        particle.gridIdx = ivec3(min(max(particle.pos.x, 0.f), containerBounds.x),
//                                 min(max(particle.pos.y, 0.f), containerBounds.y),
//                                 min(max(particle.pos.z, 0.f), containerBounds.z));
////        grid->P->setCellMark(particle.gridIdx.x, particle.gridIdx.y, particle.gridIdx.z, FLUID );
//    }
    this->clearGrid();
}
