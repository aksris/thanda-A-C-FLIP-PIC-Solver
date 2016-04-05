//
//  fluidSolver.cpp
//  Thanda
#define DEBUG
#define VISCOSITY 0.f
#define EPSILON 0.00001f

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
    mass = 10.f;
    life = 1.f;
    cameradistance = 10.f;
    density = 1.f;
}

MACGrid::MACGrid(){
}

void MACGrid::initialize(){
    vel_U.MACGridDataInitialize();
    vel_V.MACGridDataInitialize();
    vel_W.MACGridDataInitialize();
    P.MACGridDataInitialize();
    save_kernel_wt_U.MACGridDataInitialize();
    save_kernel_wt_V.MACGridDataInitialize();
    save_kernel_wt_W.MACGridDataInitialize();
}

MACGrid& MACGrid::operator =(const MACGrid& val){
    if (&val == this)
    {
        return *this;
    }
    vel_U.data = val.vel_U.data;
    vel_U.mData = val.vel_U.mData;

    vel_V.data = val.vel_V.data;
    vel_V.mData = val.vel_V.mData;

    vel_W.data = val.vel_W.data;
    vel_W.mData = val.vel_W.mData;

    P.data = val.P.data;
    P.mData = val.P.mData;

    save_kernel_wt_U = val.save_kernel_wt_U;
    save_kernel_wt_V = val.save_kernel_wt_V;
    save_kernel_wt_W = val.save_kernel_wt_W;
    return *this;
}

FluidSolver::FluidSolver(){
    LastUsedParticle = 0;
    MaxParticles = 200000;
    delta = (1.f/24);
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

void FluidSolver::constructMACGrid(glm::vec3 containerBoundsOri){
    containerBounds = containerBoundsOri;

    num_cells = (containerBounds.x * containerBounds.y * containerBounds.z);
    grid.initialize();
    i_size = containerBounds.x;
    j_size = containerBounds.y;
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
                while(iter < 8){
//                    p.pos = vec3(1.8,1.8,1.6);
                    p.pos = vec3(rand(i,i+1), rand(j,j+1), rand(k, k+1));
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
    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;
    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                tmp.vel_U.setCell(i, j, k, grid.vel_U(i, j, k) - tmp.vel_U(i, j, k));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                tmp.vel_V.setCell(i, j, k, grid.vel_V(i, j, k) - tmp.vel_V(i, j, k)) ;
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z + 1; ++k){
                tmp.vel_W.setCell(i, j, k, grid.vel_W(i, j, k) - tmp.vel_W(i, j, k));
            }
        }
    }

    particle_save = ParticlesContainer;
    //for every particle, set the change in velocity + current particle velocity
    //interpolate
    for(int i = 0; i < particle_save.size(); i++){
        particle_save.at(i).speed.x = tmp.vel_U.interpolate(particle_save.at(i).pos);
        particle_save.at(i).speed.y = tmp.vel_V.interpolate(particle_save.at(i).pos);
        particle_save.at(i).speed.z = tmp.vel_W.interpolate(particle_save.at(i).pos);
    }

}

void FluidSolver::PicSolve(){
    particle_save_pic = ParticlesContainer;
    //for every particle, set the new grid velocity
    //interpolate
    for(int i = 0; i < particle_save_pic.size(); i++){
        particle_save_pic.at(i).speed.x = grid.vel_U.interpolate(particle_save_pic.at(i).pos);
        particle_save_pic.at(i).speed.y = grid.vel_V.interpolate(particle_save_pic.at(i).pos);
        particle_save_pic.at(i).speed.z = grid.vel_W.interpolate(particle_save_pic.at(i).pos);
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
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
            for (int k = 0; k < 5; k++){
                speed = glm::vec3(0.f, -9.8f , 0.f) * delta;
                grid.vel_V.setCellAdd(i,j,k, speed.y);
            }
        }
    }
}

void FluidSolver::storeParticleVelocityToGrid(){
    //for all the grid indices, calculate vel_u, vel_v, vel_w
    vec3 index, pos, r; int x,y,z;
    float h = 1.f, weight = 0.f;
    float fract_partx,fract_party,fract_partz;
    for(int i = 0; i < ParticlesContainer.size(); ++i){
        vec3 w_pos = ParticlesContainer.at(i).pos;
        index = ParticlesContainer.at(i).gridIdx;
        x = index.x, y = index.y, z = index.z;

        grid.P.setCellMark(x, y, z, FLUID, true);

        //################## X direction #############################################
        pos = grid.vel_U.worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_U.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_U.setCellAdd(x+1, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_U.setCellAdd(x+1, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_U.setCellAdd(x, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid.vel_U.setCellAdd(x, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_U.setCellAdd(x, y, z+1, ParticlesContainer.at(i).speed.y * weight );
        grid.save_kernel_wt_U.setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_U.setCellAdd(x+1, y, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid.vel_U.setCellAdd(x+1, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_U.setCellAdd(x+1, y+1, z+1, weight);

        //######################## y direction splatting #####################################
        pos = grid.vel_V.worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_V.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_V.setCellAdd(x+1, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_V.setCellAdd(x+1, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_V.setCellAdd(x, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid.vel_V.setCellAdd(x, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_V.setCellAdd(x, y, z+1, ParticlesContainer.at(i).speed.y * weight );
        grid.save_kernel_wt_V.setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_V.setCellAdd(x+1, y, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid.vel_V.setCellAdd(x+1, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_V.setCellAdd(x+1, y+1, z+1, weight);

        //############### z direction #######################################################
        pos = grid.vel_W.worldToLocal(w_pos);
        x = pos.x, y = pos.y, z = pos.z;
        fract_partx = (pos[0] - x);
        fract_party = (pos[1] - y);
        fract_partz = (pos[2] - z);

        //splatting to all neighbors that are/will be involved in trilerp
        weight = (1-fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_W.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x, y, z, weight);

        weight = (fract_partx)*(1-fract_party)*(1-fract_partz);
        grid.vel_W.setCellAdd(x+1, y, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x+1, y, z, weight);

        weight = (fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_W.setCellAdd(x+1, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x+1, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(1-fract_partz);
        grid.vel_W.setCellAdd(x, y+1, z, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x, y+1, z, weight);

        weight = (1-fract_partx)*(fract_party)*(fract_partz);
        grid.vel_W.setCellAdd(x, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x, y+1, z+1, weight);

        weight = (1-fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_W.setCellAdd(x, y, z+1, ParticlesContainer.at(i).speed.y * weight );
        grid.save_kernel_wt_W.setCellAdd(x, y, z+1, weight);

        weight = (fract_partx)*(1-fract_party)*(fract_partz);
        grid.vel_W.setCellAdd(x+1, y, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x+1, y, z+1, weight);

        weight = (fract_partx)*(fract_party)*(fract_partz);
        grid.vel_W.setCellAdd(x+1, y+1, z+1, ParticlesContainer.at(i).speed.y * weight);
        grid.save_kernel_wt_W.setCellAdd(x+1, y+1, z+1, weight);

    }

    x = grid.P.containerBounds.x;
    y = grid.P.containerBounds.y;
    z = grid.P.containerBounds.z;
    for(int i = 0; i < x+1; ++i){
        for(int j = 0; j < y ; ++j){
            for(int k = 0; k < z; ++k){
                grid.vel_U.setCell(i, j, k, grid.vel_U(i,j,k) / max(grid.save_kernel_wt_U(i,j,k), EPSILON));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                grid.vel_V.setCell(i, j, k, grid.vel_V(i,j,k) / max(grid.save_kernel_wt_V(i,j,k), EPSILON));
            }
        }
    }
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y ; ++j){
            for(int k = 0; k < z+1; ++k){
                grid.vel_W.setCell(i, j, k, grid.vel_W(i,j,k) / max(grid.save_kernel_wt_W(i,j,k), EPSILON));
            }
        }
    }

    initializeMarkerGrid();
}

void FluidSolver::initializeMarkerGrid(){
    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;

    //neighborhood of 0, based on distance

    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(i == 0 || i == x - 1 || j == 0 || j == y - 1 || k == 0 || k == z - 1)
                {
                    grid.P.setCellMark(i, j, k, SOLID, true);
                }
            }
        }
    }
}


void FluidSolver::insertCoefficient(int id, int i, int j, int k, double w, std::vector<Eigen::Triplet<double>> &coeffs, int n){
    n = 5;
    int id1 = k * n * n + j * n + i;
    coeffs.push_back(Eigen::Triplet<double>(id, id1, w));
}

void FluidSolver::buildMatrixA(std::vector<Eigen::Triplet<double>>& coefficients, long n){
    n = 5;
    for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 5; ++j){
            for(int k = 0; k < 5; ++k){
                int id = k * n * n + j * n + i; //id for matrix
                int scale = 1;
                if(grid.P.getCellMark(i,j,k) == FLUID){
                    //x
                    if(i > 0 && grid.P.getCellMark(i-1,j,k) == FLUID){
                        scale++;
                        insertCoefficient(id, i-1,j,k, -1, coefficients, n);
                    }
                    if(i < n && grid.P.getCellMark(i+1,j,k) == FLUID){
                        scale++;
                        insertCoefficient(id, i+1,j,k, -1, coefficients, n);
                    }
                    else if(i < n && grid.P.getCellMark(i+1,j,k) == AIR){
                        scale++;
                    }

                    //y
                    if(j > 0 && grid.P.getCellMark(i,j-1,k) == FLUID){
                        scale++;
                        insertCoefficient(id, i,j-1,k, -1, coefficients, n);
                    }
                    if(j < n && grid.P.getCellMark(i,j+1,k) == FLUID){
                        scale++;
                        insertCoefficient(id, i,j+1,k, -1, coefficients, n);
                    }
                    else if(j < n && grid.P.getCellMark(i,j+1,k) == AIR){
                        scale++;
                    }

                    //z
                    if(k > 0 && grid.P.getCellMark(i,j,k-1) == FLUID){
                        scale++;
                        insertCoefficient(id, i,j,k-1, -1, coefficients, n);
                    }
                    if(k < n && grid.P.getCellMark(i,j,k+1) == FLUID){
                        scale++;
                        insertCoefficient(id, i,j,k+1, -1, coefficients, n);
                    }
                    else if(k < n && grid.P.getCellMark(i,j,k+1) == AIR){
                        scale++;
                    }
                }
                insertCoefficient(id, i,j,k, scale, coefficients, n);
            }
        }
    }
}

void FluidSolver::buildDivergences(Eigen::VectorXd& u, int n){
    float h = 1.f;
    int iter = 0;
    for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 5; ++j){
            for(int k = 0; k < 5; ++k){
                double divergence = (grid.vel_U(i+1, j, k) -
                                    grid.vel_U(i, j, k) +
                                    grid.vel_V(i, j+1, k) -
                                    grid.vel_V(i, j, k) +
                                    grid.vel_W(i, j, k+1) -
                                    grid.vel_W(i, j, k)
                                    ) / h;

                u[iter++] = -divergence;
            }
        }
    }
}

void FluidSolver::fillPressureGrid(Eigen::VectorXd x, int n){
    int iter = 0;
    n = 5;
    for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 5; ++j){
            for(int k = 0; k < 5; ++k){
                grid.P.setCell(i,j,k, x[iter++]);
#ifdef DEBUG
                if(grid.P.getCellMark(i,j,k) == FLUID){
                    float x = grid.P(i,j,k);
                    std::cout << "Pressure at "<< i << " "<< j << " " << k << ": "<<x<< std::endl;
                }
#endif
            }
        }
    }
}

void FluidSolver::ProjectPressure(){
    int x = containerBounds.x, y = containerBounds.y, z = containerBounds.z;
    long m = x * y * z;

    std::vector<Eigen::Triplet<double>> coefficients;
    Eigen::VectorXd u(m);

    buildDivergences(u, m);
    buildMatrixA(coefficients, m);
    Eigen::SparseMatrix<double> A(m,m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    //solve
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> chol;  // performs a Cholesky factorization of A
    chol.compute(A);
    Eigen::VectorXd p = chol.solve(u);
    fillPressureGrid(p, m);

    SubtractPressureGradient();
}

void FluidSolver::SubtractPressureGradient(){
    float density = 1000.f;
    float dx = 1.f;
    float scale = delta / (density * dx);

    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;

    //loop over i,j,k
    //update vel_U
    for(int i = 0; i < x+1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(grid.P.getCellMark(i-1,j,k) == FLUID || grid.P.getCellMark(i,j,k) == FLUID){
                    if(grid.P.getCellMark(i-1,j,k) == SOLID || grid.P.getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid.vel_U.setCell(i,j,k,0.f);
                    }
                    else{
                        grid.vel_U.setCellAdd(i,j,k, -(scale * (grid.P(i,j,k) - grid.P(i-1,j,k))));
                    }
                }
                else{
                    //set as unknown
//                    grid.vel_U.setCell(i,j,k, UNKNOWN);
                }
            }
        }
    }
    //update vel_V
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y+1; ++j){
            for(int k = 0; k < z; ++k){
                if(grid.P.getCellMark(i,j-1,k) == FLUID || grid.P.getCellMark(i,j,k) == FLUID){
                    if(grid.P.getCellMark(i,j-1,k) == SOLID || grid.P.getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid.vel_V.setCell(i,j,k,0.f);
                    }
                    else{
                        grid.vel_V.setCellAdd(i,j,k, -(scale * (grid.P(i,j,k) - grid.P(i,j-1,k))));
                    }
                }
                else{
                    //set as unknown
//                    grid.vel_V.setCell(i,j,k, UNKNOWN);
                }
            }
        }
    }
    //update vel_W
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z+1; ++k){
                if(grid.P.getCellMark(i,j,k-1) == FLUID || grid.P.getCellMark(i,j,k) == FLUID){
                    if(grid.P.getCellMark(i,j,k-1) == SOLID || grid.P.getCellMark(i,j,k) == SOLID){
                        //set cell value as usolid at i,j,k
                        grid.vel_W.setCell(i,j,k,0.f);
                    }
                    else{
                        grid.vel_W.setCellAdd(i,j,k, -(scale * (grid.P(i,j,k) - grid.P(i,j,k-1))));
                    }
                }
                else{
                    //set as unknown
//                    grid.vel_W.setCell(i,j,k, UNKNOWN);
                }
            }
        }
    }
}

vec3 EulerStep (const vec3 pos, const vec3 speed, float time_step){
    return pos + speed * time_step;
}

vec3 FluidSolver::integratePos(const vec3 pos, const vec3 speed, float time_step, bool RK2){
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
    return new_pos;
}

void FluidSolver::ExtrapolateVelocity(){
    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;

    MACGridData extrapolate_grid_U;
    MACGridData extrapolate_grid_V;
    MACGridData extrapolate_grid_W;
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
                if((i < x && grid.P.getCellMark(i,j,k) == FLUID) || (i > 0 && grid.P.getCellMark(i-1,j,k) == FLUID)){
                    extrapolate_grid_U.setCellMark(i,j,k, FLUID, true);
                }
                //marking solid for extrapolation
                else if((i < x && (grid.P.getCellMark(i,j,k) == SOLID)) &&
                        (i > 0 && (grid.P.getCellMark(i-1,j,k) == SOLID )) //check if me solid and near solid
                        || (i < x && (grid.P.getCellMark(i,j,k) == AIR)) &&
                        (i > 0 && (grid.P.getCellMark(i-1,j,k) == SOLID || grid.P.getCellMark(i-1,j,k) == AIR))) //me air, and near (solid or air)
                {
                    extrapolate_grid_U.setCellMark(i,j,k,SOLID, true);
                }
            }
        }
    }
    //y direction
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y+1; ++j){
            for(int k = 0; k < z; ++k){
                //marking fluid or near fluid cells as FLUID
                if((j < y && grid.P.getCellMark(i,j,k) == FLUID) || (j > 0 && grid.P.getCellMark(i,j-1,k) == FLUID)){
                    extrapolate_grid_V.setCellMark(i,j,k, FLUID, true);
                }
                //marking solid for extrapolation
                else if((j < y && (grid.P.getCellMark(i,j,k) == SOLID)) &&
                        (j > 0 && (grid.P.getCellMark(i,j-1,k) == SOLID )) //check if me solid and near solid
                        || (j < y && (grid.P.getCellMark(i,j,k) == AIR)) &&
                        (j > 0 && (grid.P.getCellMark(i,j-1,k) == SOLID || grid.P.getCellMark(i,j-1,k) == AIR))) //me air, and near (solid or air)
                {
                    extrapolate_grid_V.setCellMark(i,j,k,SOLID, true);
                }
            }
        }
    }
    //z direction
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z+1; ++k){
                //marking fluid or near fluid cells as FLUID
                if((k < z && grid.P.getCellMark(i,j,k) == FLUID) || (k > 0 && grid.P.getCellMark(i,j,k-1) == FLUID)){
                    extrapolate_grid_W.setCellMark(i,j,k, FLUID, true);
                }
                //marking solid for extrapolation
                else if((k < z && (grid.P.getCellMark(i,j,k) == SOLID)) &&
                        (k > 0 && (grid.P.getCellMark(i,j,k-1) == SOLID )) //check if me solid and near solid
                        || (k < z && (grid.P.getCellMark(i,j,k) == AIR)) &&
                        (k > 0 && (grid.P.getCellMark(i,j,k-1) == SOLID || grid.P.getCellMark(i,j,k-1) == AIR))) //me air, and near (solid or air)
                {
                    extrapolate_grid_W.setCellMark(i,j,k,SOLID, true);
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
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid.vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;
                                sum += grid.vel_U(q[qk][0],q[qk][1],
                                        q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        grid.vel_U.setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid.vel_V(i,j+1,k) << std::endl;
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
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid.vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;

                                sum += grid.vel_V(q[qk][0],q[qk][1],
                                        q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " : "<< sum/wsum<<std::endl;
                        grid.vel_V.setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid.vel_V(i,j+1,k) << std::endl;
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
                                //std::cout << "[thanda] ExtrapolateVelocity for "<< i << " "<< j << " "<< k << " from " << q[qk][0] << " "<< q[qk][1] << " "<< q[qk][2] << " : " << this->grid.vel_V(q[qk][0], q[qk][1], q[qk][2]) << std::endl;

                                wsum ++;

                                sum += grid.vel_W(q[qk][0],q[qk][1],
                                        q[qk][2]);
                            }
                        }
                    }
                    if(wsum){
                        grid.vel_W.setCell(i,j,k,sum/wsum);
                    }
                    //std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid.vel_V(i,j+1,k) << std::endl;
                }
            }

        }
    }
}


void FluidSolver::setBoundaryVelocitiesToZero(const glm::vec3 containerBounds){

    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;
    //neighborhood of 0, based on distance

    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z; ++k){
                if(i == 0 || i == x){
                    grid.vel_U.setCell(i, j, k, 0.f);
                }
                if (i > 0 && i < x && ((grid.P.getCellMark(i,j,k) == SOLID && grid.P.getCellIndex(i-1,j,k) == FLUID) ||
                                       (grid.P.getCellMark(i,j,k) == FLUID && grid.P.getCellIndex(i-1,j,k) == SOLID))){
                    grid.vel_U.setCell(i, j, k, 0.f);
                }
            }
        }
    }


    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                if(j == 0 || j == y){
                    grid.vel_V.setCell(i, j, k, 0.f);
                }
                if (j > 0 && j < y && ((grid.P.getCellMark(i,j,k) == SOLID && grid.P.getCellIndex(i,j-1,k) == FLUID) ||
                                       (grid.P.getCellMark(i,j,k) == FLUID && grid.P.getCellIndex(i,j,k) == SOLID))){
                    grid.vel_V.setCell(i, j, k, 0.f);
                }
            }
        }
    }

    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z + 1; ++k){
                if(k == 0 || k == z){
                    grid.vel_W.setCell(i, j, k, 0.f);
                }
                if (k > 0 && k < z && ((grid.P.getCellMark(i,j,k) == SOLID && grid.P.getCellIndex(i,j,k-1) == FLUID) ||
                                       (grid.P.getCellMark(i,j,k) == FLUID && grid.P.getCellIndex(i,j,k-1) == SOLID))){
                    grid.vel_W.setCell(i, j, k, 0.f);
                }
            }
        }
    }
}

void FluidSolver::storeCurrentGridVelocities(){
    tmp = grid;
}

void FluidSolver::clearGrid(){
    std::fill(grid.save_kernel_wt_U.data.begin(), grid.save_kernel_wt_U.data.end(), 0.f);
    std::fill(grid.save_kernel_wt_V.data.begin(), grid.save_kernel_wt_V.data.end(), 0.f);
    std::fill(grid.save_kernel_wt_W.data.begin(), grid.save_kernel_wt_W.data.end(), 0.f);
    std::fill(grid.vel_U.data.begin(), grid.vel_U.data.end(), 0.f);
    std::fill(grid.vel_V.data.begin(), grid.vel_V.data.end(), 0.f);
    std::fill(grid.vel_W.data.begin(), grid.vel_W.data.end(), 0.f);
    std::fill(grid.P.data.begin(), grid.P.data.end(), 0.f);
    std::fill(grid.P.mData.begin(), grid.P.mData.end(), 0.f);
}

void FluidSolver::step(){

//#define DEBUG
    
    // Step 3 - Store Particle Velocity at current time step to MACGrid
    int i = this->ParticlesContainer.at(0).gridIdx[0];
    int j = this->ParticlesContainer.at(0).gridIdx[1];
    int k = this->ParticlesContainer.at(0).gridIdx[2];
    
    this->storeParticleVelocityToGrid();
#ifdef DEBUG
        std::cout << "[thanda] storeParticleVelocityToGrid Iteration" << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j, k) << " at idx " << i << " "<< j << " "<< k << " : "<< this->grid.vel_V(i,j,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j+1, k) << " at idx " << i << " "<< j+1 << " "<< k << " : "<< this->grid.vel_V(i,j+1,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i+1, j, k) << " at idx " << i+1 << " "<< j << " "<< k << " : "<< this->grid.vel_V(i+1,j,k) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i+1, j+1, k) << " at idx " << i+1 << " "<< j+1 << " "<< k << " : "<< this->grid.vel_V(i+1,j+1,k) << std::endl; // is this valid

    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j, k+1) << " at idx " << i << " "<< j << " "<< k+1 << " : "<< this->grid.vel_V(i,j,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i, j+1, k+1) << " at idx " << i << " "<< j+1 << " "<< k+1 << " : "<< this->grid.vel_V(i,j+1,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i+1, j, k+1) << " at idx " << i+1 << " "<< j << " "<< k+1 << " : "<< this->grid.vel_V(i+1,j,k+1) << std::endl;
    std::cout << "[thanda] Cell Type :"<< this->grid.P.getCellMark(i+1, j+1, k+1) << " at idx " << i+1 << " "<< j+1 << " "<< k+1 << " : "<< this->grid.vel_V(i+1,j+1,k+1) << std::endl; // is this valid

#endif
    
    
    // Step 4 - Add Body Forces like Gravity to MACGrid
    this->CalculateGravityToCell(delta);
    
    // Step 5 - Store a temporary copy of Grid Velocities for FLIP
    this->storeCurrentGridVelocities();
    
    this->setBoundaryVelocitiesToZero(containerBounds);
    
    //pressure solve
    this->ProjectPressure();
    
    this->setBoundaryVelocitiesToZero(containerBounds);
    this->ExtrapolateVelocity();

    // Step  - Calculate new flip & pic velocities for each particle
//    this->FlipSolve();
    this->PicSolve();

    // Step - Lerp(FLIPVelocity, PICVelocity, 0.95)
    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
        this->ParticlesContainer.at(i).speed = (1.f - VISCOSITY) * this->particle_save_pic.at(i).speed;
//        + VISCOSITY * this->particle_save.at(i).speed;
        this->ParticlesContainer.at(i).pos = this->integratePos(this->ParticlesContainer.at(i).pos,
                                                                this->ParticlesContainer.at(i).speed, delta, true);
    }

    

//     Step - Collision Response
    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
        
        Particle& particle = this->ParticlesContainer.at(i); // careful - by reference;
        
        float dampingFactor = 0.5f;
        
//            if(particle.pos.x  < EPSILON || particle.pos.x  > upperBounds.x - EPSILON)
//            {
//                currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(-dampingFactor,1,1));
//                currParticle.setPredictedPosition(currParticle.getPosition() + timeStep * currParticle.getVelocity());
//            }
            
        if (particle.pos.y < EPSILON){
            if(particle.speed.y < EPSILON){
                particle.speed *= (vec3(1.f, -dampingFactor, 1.f));
            }
            particle.pos.y = EPSILON ;
        }

        else if(particle.pos.y  > containerBounds.y - EPSILON){// does not work
            if(particle.speed.y > EPSILON){
                particle.speed *= (vec3(1.f, -dampingFactor, 1.f));
            }
            particle.pos.y = containerBounds.y - EPSILON ;
        }
    
            
//            if(particle.pos.z  < lowerBounds.z + EPSILON || particle.pos.z  > upperBounds.z - EPSILON)
//            {
//                currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(1,1,-dampingFactor));
//                currParticle.setPredictedPosition(currParticle.getPosition() + timeStep * currParticle.getVelocity());
//            }
            
            
        particle.pos = this->integratePos(particle.pos, particle.speed, delta, true);
        particle.gridIdx = glm::ivec3(particle.pos);
    }
    this->clearGrid();
}
