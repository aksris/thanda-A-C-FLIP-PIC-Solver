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
    save_kernel_wt_V.MACGridDataInitialize();
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

    save_kernel_wt_V = val.save_kernel_wt_V;
    return *this;
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
    for(int i = 1; i < (int)boundx+1; i++){
        for(int j = 1; j < (int)boundy+1; j++){
            for(int k = 1; k < (int)boundz+1; k++){
                int iter = 0;
                while(iter < 8){
                    p.pos = vec3(rand(i,i+1), rand(j,j+1), rand(k, k+1));
                    p.speed = vec3(0.f, -8.f, 0.f);
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
//        particle_save_pic.at(i).speed.x = grid.vel_U.interpolate(particle_save_pic.at(i).pos);
        particle_save_pic.at(i).speed.y = grid.vel_V.interpolate(particle_save_pic.at(i).pos);
//        particle_save_pic.at(i).speed.z = grid.vel_W.interpolate(particle_save_pic.at(i).pos);
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
                speed = glm::vec3(0.f, -10.f , 0.f) * delta;
                grid.vel_V.setCellAdd(i,j,k, speed.y);
            }
        }
    }
}

void FluidSolver::storeParticleVelocityToGrid(){
    //for all the grid indices, calculate vel_u, vel_v, vel_w
    float h = 1.f;
    for(int i = 0; i < ParticlesContainer.size(); ++i){
        vec3 index = ParticlesContainer.at(i).gridIdx;
        int x = index.x, y = index.y, z = index.z;

        grid.P.setCellMark(x, y, z, FLUID, true);

        vec3 r = ParticlesContainer.at(i).pos - (index + vec3(0.f, 0.5f, 0.5f));
        grid.vel_U.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.x * StiffKernel(r, h) * 0.125f);


        //y direction splatting
        r = ParticlesContainer.at(i).pos - (index + vec3(0.5f, 0.f, 0.5f));

        grid.vel_V.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x, y, z, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x+1, y, z) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x+1, y, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x+1, y, z, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x-1, y, z) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x-1, y, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x-1, y, z, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x, y+1, z) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x, y+1, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x, y+1, z, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x, y-1, z) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x, y-1, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x, y-1, z, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x, y, z+1) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x, y, z+1, ParticlesContainer.at(i).speed.y * StiffKernel(r, h) );
        grid.save_kernel_wt_V.setCellAdd(x, y, z+1, StiffKernel(r, h));

        r = ParticlesContainer.at(i).pos - (vec3(x, y, z-1) + vec3(0.5f, 0.f, 0.5f));
        grid.vel_V.setCellAdd(x, y, z-1, ParticlesContainer.at(i).speed.y * StiffKernel(r, h));
        grid.save_kernel_wt_V.setCellAdd(x, y, z-1, StiffKernel(r, h));


        r = ParticlesContainer.at(i).pos - (index + vec3(0.5f, 0.5f, 0.f));
        grid.vel_W.setCellAdd(x, y, z, ParticlesContainer.at(i).speed.y * StiffKernel(r, h) * 0.125f);

    }

    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                grid.vel_V.setCell(i, j, k, grid.vel_V(i,j,k) / max(grid.save_kernel_wt_V(i,j,k), EPSILON));
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
                if(grid.P.getCellMark(i, j, k) != FLUID){
                    grid.P.setCellMark(i, j, k, AIR, true);
                    if(i == 0 || i == x - 1 || j == 0 || j == y - 1 || k == 0 || k == z - 1)
                    {
                        grid.P.setCellMark(i, j, k, SOLID, true);
                    }
                }
            }
        }
    }
}

void FluidSolver::calculateNewGridVelocities(){
    float h = 1.4f;
    int x = grid.P.containerBounds.x;
    int y = grid.P.containerBounds.y;
    int z = grid.P.containerBounds.z;
    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y; ++j){
            for( int k = 0; k < z; ++k){
                if(i>0 && i<x){
                    float pf = grid.P(i,j,k);
                    float pb = grid.P(i-1,j,k);
                    pf = grid.P(i,j,k);
                    pb = grid.P(i-1,j,k);
                    float xval = grid.vel_U(i,j,k);
                    xval -= (pf-pb)/h;
                    grid.vel_U.setCell(i,j,k,xval);
                }
            }
        }
    }

    for( int i = 0; i < x ; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                if(j>0 && j<y){
                    float pf = grid.P(i,j,k);
                    float pb = grid.P(i,j-1,k);
                    pf = grid.P(i,j,k);
                    pb = grid.P(i,j-1,k);
                    float xval = grid.vel_V(i,j,k);
                    xval -= (pf-pb)/h;
                    grid.vel_V.setCell(i,j,k,xval);
                }
            }
        }
    }

    for( int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            for(int k = 0; k < z + 1; ++k){
                if(k>0 && k<z){
                    float pf = grid.P(i,j,k);
                    float pb = grid.P(i,j,k-1);
                    pf = grid.P(i,j,k);
                    pb = grid.P(i,j,k-1);
                    float xval = grid.vel_W(i,j,k);
                    xval -= (pf-pb)/h;
                    grid.vel_W.setCell(i,j,k,xval);
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

                u[iter++] = divergence;
            }
        }
    }
}

void FluidSolver::fillPressureGrid(Eigen::VectorXd x, int n){
    int iter = 0;
    n = 5;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                grid.P.setCell(i,j,k, x[iter++]);
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
    //neighborhood of 0, based on distance
    for(int i = 0; i < x + 1; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z + 1; ++k){
                for(unsigned int n = 0; n < 3; ++n){
                    if(n!=0 && i>x-1){ continue; };
                    if(n!=1 && j>y-1){ continue; };
                    if(n!=2 && k>z-1){ continue; };
                    if(grid.P.getCellMark(i, j, k) == SOLID || grid.P.getCellMark(i, j, k) == AIR ){
                        unsigned int wsum = 0;
                        float sum = 0.0f;
                        glm::vec3 q[6] = { glm::vec3(i-1,j,k), glm::vec3(i+1,j,k),
                                           glm::vec3(i,j-1,k), glm::vec3(i,j+1,k),
                                           glm::vec3(i,j,k-1), glm::vec3(i,j,k+1) };
                        for(unsigned int qk = 0; qk < 6; ++qk){
                            if(q[qk][0] >= 0 && q[qk][0]< x +(n==0) && q[qk][1] >= 0 &&
                                    q[qk][1] < y+(n==1) && q[qk][2] >= 0 && q[qk][2] < z+(n==2) ) {
                                if(grid.P.getCellMark(q[qk][0], q[qk][1], q[qk][2]) == FLUID){
                                    wsum ++;
                                    if(n == 0){
                                        sum += grid.vel_U(q[qk][0],q[qk][1],
                                                q[qk][2]);
                                    }else if(n == 1){
                                        if(grid.vel_V(q[qk][0],q[qk][1],
                                                      q[qk][2]) > -EPSILON)
                                        {
                                            int a = 1;
                                            sum += grid.vel_V(q[qk][0],q[qk][1],
                                                    q[qk][2]);
                                        }
                                    }else if(n == 2){
                                        sum += grid.vel_W(q[qk][0],q[qk][1],
                                                q[qk][2]);
                                    }
                                }
                            }
                        }
                        if(wsum){
                            if(n==0){
                                grid.vel_U.setCell(i,j,k,sum/wsum);
                            }else if(n==1){
                                grid.vel_V.setCell(i,j,k,sum/wsum);
                            }else if(n==2){
                                grid.vel_W.setCell(i,j,k,sum/wsum);
                            }
                        }
                    }
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
            }
        }
    }


    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y + 1; ++j){
            for(int k = 0; k < z; ++k){
                if(j == 0 || j == y){
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
            }
        }
    }
}

void FluidSolver::storeCurrentGridVelocities(){
    tmp = grid;
}

void FluidSolver::clearGrid(){
    std::fill(grid.vel_U.data.begin(), grid.vel_U.data.end(), 0.f);
    std::fill(grid.vel_V.data.begin(), grid.vel_V.data.end(), 0.f);
    std::fill(grid.vel_W.data.begin(), grid.vel_W.data.end(), 0.f);
    std::fill(grid.P.data.begin(), grid.P.data.end(), 0.f);
    std::fill(grid.P.mData.begin(), grid.P.mData.end(), 0.f);
}

void FluidSolver::step(){
    // Step 3 - Store Particle Velocity at current time step to MACGrid
    std::cout << "kuch text " << std::endl;
    this->storeParticleVelocityToGrid();
    std::cout << this->grid.vel_V.data.at(31) << std::endl;
    
    // Step 4 - Add Body Forces like Gravity to MACGrid
//    this->CalculateGravityToCell(delta);
    
    // Step 5 - Store a temporary copy of Grid Velocities for FLIP
    this->storeCurrentGridVelocities();
    
//    this->setBoundaryVelocitiesToZero(vec3(5.f));
    //pressure solve
//    this->ProjectPressure();
    
//    this->calculateNewGridVelocities();
//    this->setBoundaryVelocitiesToZero(vec3(5.f));
    this->ExtrapolateVelocity();
    std::cout << this->grid.vel_V.data.at(31) << std::endl;

    // Step  - Calculate new flip & pic velocities for each particle
//    this->FlipSolve();
    this->PicSolve();
    std::cout << this->grid.vel_V.data.at(31) << std::endl;

    // Step - Lerp(FLIPVelocity, PICVelocity, 0.95)
//    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
//        this->ParticlesContainer.at(i).speed = (1.f - VISCOSITY) * this->particle_save_pic.at(i).speed
//        /*+ VISCOSITY * this->particle_save.at(i).speed*/;
//        this->ParticlesContainer.at(i).pos = this->integratePos(this->ParticlesContainer.at(i).pos,
//                                                                this->ParticlesContainer.at(i).speed, delta, true);
//        if(this->ParticlesContainer.at(i).gridIdx == vec3(1.f)){
//            std::cout << this->ParticlesContainer.at(i).speed.y << std::endl;
//        }
//    }

    // Step - Collision Response
//    for(int i = 0; i < this->ParticlesContainer.size(); ++i){
//        if (this->ParticlesContainer.at(i).pos.y < EPSILON){
//            if(this->ParticlesContainer.at(i).speed.y < EPSILON){
//                this->ParticlesContainer.at(i).speed *= (vec3(1.f, -1.f, 1.f));
//            }
//            this->ParticlesContainer.at(i).pos.y = EPSILON ;
//            this->ParticlesContainer.at(i).pos = this->integratePos(this->ParticlesContainer.at(i).pos,
//                                                                    this->ParticlesContainer.at(i).speed, delta, true);
//        }
//    }
    this->clearGrid();
}
