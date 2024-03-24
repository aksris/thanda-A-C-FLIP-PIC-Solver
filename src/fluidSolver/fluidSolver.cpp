//
//  fluidSolver.cpp
//  Thanda

/* Debugging and profiling are under same macro: DEBUG */
//#define DEBUG
#define EPSILON 0.00001f

#include <algorithm>
#include <execution>

#include "fluidSolver.hpp"
#include "../scene/scene.hpp"


/*
 Random number generator between the range LO to HI
*/
float rand(int LO, int HI) {
	return LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
}

Particle::Particle() {
	pos = glm::vec3(0.f, 0.f, 0.f);
	speed = glm::vec3(0.f, 0.f, 0.f);
	r = 0;
	g = 0;
	b = 220; //blue colored partciles to look more fluid-y
	a = 230;
	size = 0.05f;
	cameradistance = 0.f;
}

MACGrid::MACGrid(const ivec3& resolution, const vec3& containerBounds, float cellSize) {
	vel_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
	vel_V = new MACGridDataY(ivec3(resolution.x, resolution.y + 1, resolution.z), containerBounds, cellSize);
	vel_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z + 1), containerBounds, cellSize);

	flip_vel_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
	flip_vel_V = new MACGridDataY(ivec3(resolution.x, resolution.y + 1, resolution.z), containerBounds, cellSize);
	flip_vel_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z + 1), containerBounds, cellSize);

	save_kernel_wt_U = new MACGridDataX(ivec3(resolution.x + 1, resolution.y, resolution.z), containerBounds, cellSize);
	save_kernel_wt_V = new MACGridDataY(ivec3(resolution.x, resolution.y + 1, resolution.z), containerBounds, cellSize);
	save_kernel_wt_W = new MACGridDataZ(ivec3(resolution.x, resolution.y, resolution.z + 1), containerBounds, cellSize);

	P = new MACGridData(ivec3(resolution.x, resolution.y, resolution.z), containerBounds, cellSize);
}

MACGrid::~MACGrid() {
	delete vel_U;
	delete vel_V;
	delete vel_W;

	delete flip_vel_U;
	delete flip_vel_V;
	delete flip_vel_W;

	delete save_kernel_wt_U;
	delete save_kernel_wt_V;
	delete save_kernel_wt_W;

	delete P;
}

void MACGrid::initialize() {
	vel_U->MACGridDataInitialize();
	vel_V->MACGridDataInitialize();
	vel_W->MACGridDataInitialize();
	P->MACGridDataInitialize();
	save_kernel_wt_U->MACGridDataInitialize();
	save_kernel_wt_V->MACGridDataInitialize();
	save_kernel_wt_W->MACGridDataInitialize();
}

MACGrid& MACGrid::operator =(const MACGrid& val) {
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

FluidSolver::FluidSolver(const Scene& iScene) {
	MaxParticles = 5000000;

	this->scene = iScene;
	this->resolution = scene.resolution;
	this->containerBounds = scene.containerBounds;
	this->cellSize = max(max(containerBounds.x / resolution.x, containerBounds.y / resolution.y),
		containerBounds.z / resolution.z);

	grid = new MACGrid(this->resolution, this->containerBounds, this->cellSize);
}

FluidSolver::~FluidSolver() {
	delete grid;
}

void FluidSolver::constructMACGrid() {
	resolution = scene.resolution;
	containerBounds = scene.containerBounds;

	num_cells = (resolution.x * resolution.y * resolution.z);
}

void FluidSolver::genParticles() {
	Particle p;
	int iter;
	static const int translate = 0;
	//counter for setting max particles in grid seed
	for (int k = 1; k < (int)scene.particleBounds.z; k++) {
		for (int j = 1 + translate; j < (int)scene.particleBounds.y + translate; j++) {
			for (int i = 1; i < (int)scene.particleBounds.x; i++) {
				iter = 0;
				while (iter < scene.seed) {
					p.gridIdx = ivec3(i, j, k);
					p.pos = vec3(rand(i, i + 1), rand(j, j + 1), rand(k, k + 1)) * this->cellSize;
					p.size = scene.displaySize;
					ParticlesContainer.push_back(p);
					iter++;
				}
			}
		}
	}
}

void FluidSolver::FlipSolve() {
	int x = resolution.x;
	int y = resolution.y;
	int z = resolution.z;
	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x + 1; ++i) {
				grid->flip_vel_U->setCell(i, j, k, (*grid->vel_U)(i, j, k) - (*grid->flip_vel_U)(i, j, k));
			}
		}
	}
	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y + 1; ++j) {
			for (int i = 0; i < x; ++i) {
				grid->flip_vel_V->setCell(i, j, k, (*grid->vel_V)(i, j, k) - (*grid->flip_vel_V)(i, j, k));
			}
		}
	}
	for (int k = 0; k < z + 1; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x; ++i) {
				grid->flip_vel_W->setCell(i, j, k, (*grid->vel_W)(i, j, k) - (*grid->flip_vel_W)(i, j, k));
			}
		}
	}

	particle_save = ParticlesContainer;
	//for every particle, set the change in velocity + current particle velocity
	//interpolate
	for (int i = 0; i < particle_save.size(); i++) {
		particle_save.at(i).speed.x += grid->flip_vel_U->interpolate(particle_save.at(i).pos);
		particle_save.at(i).speed.y += grid->flip_vel_V->interpolate(particle_save.at(i).pos);
		particle_save.at(i).speed.z += grid->flip_vel_W->interpolate(particle_save.at(i).pos);
	}

}

void FluidSolver::PicSolve() {
	particle_save_pic = ParticlesContainer;
	//for every particle, set the new grid velocity
	//interpolate
	for (int i = 0; i < particle_save_pic.size(); i++) {
		particle_save_pic.at(i).speed.x = grid->vel_U->interpolate(particle_save_pic.at(i).pos);
		particle_save_pic.at(i).speed.y = grid->vel_V->interpolate(particle_save_pic.at(i).pos);
		particle_save_pic.at(i).speed.z = grid->vel_W->interpolate(particle_save_pic.at(i).pos);
	}

}

void FluidSolver::CalculateGravityToCell(float delta) {
	vec3 speed;
	for (int k = 0; k < resolution.z; k++) {
		for (int j = 0; j < resolution.y; j++) {
			for (int i = 0; i < resolution.x; i++) {
				speed = glm::vec3(0.f, -scene.gravity, 0.f) * delta;
				grid->vel_V->setCellAdd(i, j, k, speed.y);
			}
		}
	}
}

void FluidSolver::storeParticleVelocityToGrid() {
	//for all the grid indices, calculate vel_u, vel_v, vel_w
	ivec3 index;
	vec3 pos, r;
	int x, y, z;
	float h = this->cellSize, weight = 0.f;
	float fract_partx, fract_party, fract_partz;

	for (int i = 0; i < ParticlesContainer.size(); ++i) {

		const Particle& par = ParticlesContainer.at(i);
		// this particle may have been advected in the previous step and it's gridIdx should be updated now

		vec3 w_pos = par.pos;
		index = ivec3(floor(par.pos.x / cellSize), floor(par.pos.y / cellSize), floor(par.pos.z / cellSize));
		ParticlesContainer.at(i).gridIdx = index;

		grid->P->setCellMark(index.x, index.y, index.z, FLUID);

		//################## X direction #############################################
		pos = grid->vel_U->worldToLocal(w_pos);
		x = floor(pos.x / cellSize), y = floor(pos.y / cellSize), z = floor(pos.z / cellSize);
		fract_partx = (pos[0] - x * cellSize);
		fract_party = (pos[1] - y * cellSize);
		fract_partz = (pos[2] - z * cellSize);

		//splatting to all neighbors that are/will be involved in trilerp
		weight = (1 - fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_U->setCellAdd(x, y, z, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x, y, z, weight);

		weight = (fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_U->setCellAdd(x + 1, y, z, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x + 1, y, z, weight);

		weight = (fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_U->setCellAdd(x + 1, y + 1, z, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x + 1, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_U->setCellAdd(x, y + 1, z, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (fract_partz);
		grid->vel_U->setCellAdd(x, y + 1, z + 1, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x, y + 1, z + 1, weight);

		weight = (1 - fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_U->setCellAdd(x, y, z + 1, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x, y, z + 1, weight);

		weight = (fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_U->setCellAdd(x + 1, y, z + 1, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x + 1, y, z + 1, weight);

		weight = (fract_partx) * (fract_party) * (fract_partz);
		grid->vel_U->setCellAdd(x + 1, y + 1, z + 1, par.speed.x * weight);
		grid->save_kernel_wt_U->setCellAdd(x + 1, y + 1, z + 1, weight);

		//######################## y direction splatting #####################################
		pos = grid->vel_V->worldToLocal(w_pos);
		x = floor(pos.x / cellSize), y = floor(pos.y / cellSize), z = floor(pos.z / cellSize);
		fract_partx = (pos[0] - x * cellSize);
		fract_party = (pos[1] - y * cellSize);
		fract_partz = (pos[2] - z * cellSize);

		//splatting to all neighbors that are/will be involved in trilerp
		float speedY = par.speed.y;
		weight = (1 - fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_V->setCellAdd(x, y, z, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x, y, z, weight);

		weight = (fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_V->setCellAdd(x + 1, y, z, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x + 1, y, z, weight);

		weight = (fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_V->setCellAdd(x + 1, y + 1, z, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x + 1, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_V->setCellAdd(x, y + 1, z, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (fract_partz);
		grid->vel_V->setCellAdd(x, y + 1, z + 1, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x, y + 1, z + 1, weight);

		weight = (1 - fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_V->setCellAdd(x, y, z + 1, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x, y, z + 1, weight);

		weight = (fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_V->setCellAdd(x + 1, y, z + 1, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x + 1, y, z + 1, weight);

		weight = (fract_partx) * (fract_party) * (fract_partz);
		grid->vel_V->setCellAdd(x + 1, y + 1, z + 1, speedY * weight);
		grid->save_kernel_wt_V->setCellAdd(x + 1, y + 1, z + 1, weight);

		//############### z direction #######################################################
		pos = grid->vel_W->worldToLocal(w_pos);
		x = floor(pos.x / cellSize), y = floor(pos.y / cellSize), z = floor(pos.z / cellSize);
		fract_partx = (pos[0] - x * cellSize);
		fract_party = (pos[1] - y * cellSize);
		fract_partz = (pos[2] - z * cellSize);

		//splatting to all neighbors that are/will be involved in trilerp
		weight = (1 - fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_W->setCellAdd(x, y, z, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x, y, z, weight);

		weight = (fract_partx) * (1 - fract_party) * (1 - fract_partz);
		grid->vel_W->setCellAdd(x + 1, y, z, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x + 1, y, z, weight);

		weight = (fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_W->setCellAdd(x + 1, y + 1, z, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x + 1, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (1 - fract_partz);
		grid->vel_W->setCellAdd(x, y + 1, z, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x, y + 1, z, weight);

		weight = (1 - fract_partx) * (fract_party) * (fract_partz);
		grid->vel_W->setCellAdd(x, y + 1, z + 1, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x, y + 1, z + 1, weight);

		weight = (1 - fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_W->setCellAdd(x, y, z + 1, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x, y, z + 1, weight);

		weight = (fract_partx) * (1 - fract_party) * (fract_partz);
		grid->vel_W->setCellAdd(x + 1, y, z + 1, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x + 1, y, z + 1, weight);

		weight = (fract_partx) * (fract_party) * (fract_partz);
		grid->vel_W->setCellAdd(x + 1, y + 1, z + 1, par.speed.z * weight);
		grid->save_kernel_wt_W->setCellAdd(x + 1, y + 1, z + 1, weight);

	}

	x = grid->P->resolution.x;
	y = grid->P->resolution.y;
	z = grid->P->resolution.z;
	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x + 1; ++i) {
				grid->vel_U->setCell(i, j, k, (*grid->vel_U)(i, j, k) / max((*grid->save_kernel_wt_U)(i, j, k), EPSILON));
			}
		}
	}
	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y + 1; ++j) {
			for (int i = 0; i < x; ++i) {
				grid->vel_V->setCell(i, j, k, (*grid->vel_V)(i, j, k) / max((*grid->save_kernel_wt_V)(i, j, k), EPSILON));
			}
		}
	}
	for (int k = 0; k < z + 1; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x; ++i) {
				grid->vel_W->setCell(i, j, k, (*grid->vel_W)(i, j, k) / max((*grid->save_kernel_wt_W)(i, j, k), EPSILON));
			}
		}
	}
}

void FluidSolver::initializeMarkerGrid() {
	int x = resolution.x;
	int y = resolution.y;
	int z = resolution.z;

	//neighborhood of 0, based on distance

	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x; ++i) {
				if (i == 0 || i == x - 1 || j == 0 || j == y - 1 || k == 0 || k == z - 1)
				{
					grid->P->setCellMark(i, j, k, SOLID);
				}
				else {
					grid->P->setCellMark(i, j, k, AIR);
				}
			}
		}
	}
}


void FluidSolver::insertCoefficient(int id, int i, int j, int k, double w, std::vector<Eigen::Triplet<double>>& coeffs) {
	int id1 = k * resolution.x * resolution.y + j * resolution.x + i;
	coeffs.push_back(Eigen::Triplet<double>(id, id1, w));
}

void FluidSolver::buildMatrixA(std::vector<Eigen::Triplet<double>>& coefficients, long n) {
	for (int k = 0; k < resolution.z; ++k) {
		for (int j = 0; j < resolution.y; ++j) {
			for (int i = 0; i < resolution.x; ++i) {
				int id = k * resolution.x * resolution.y + j * resolution.x + i; //id for matrix
				float scale = 1.f / (cellSize * cellSize);
				float Adiag = 0.f;
				if (grid->P->getCellMark(i, j, k) == FLUID) {
					//x
					if (i > 0 && grid->P->getCellMark(i - 1, j, k) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i - 1, j, k, -scale, coefficients);
					}
					if (i < n && grid->P->getCellMark(i + 1, j, k) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i + 1, j, k, -scale, coefficients);
					}
					if (i > 0 && grid->P->getCellMark(i - 1, j, k) == AIR) {
						Adiag += scale;
					}
					if (i < n && grid->P->getCellMark(i + 1, j, k) == AIR) { //checking Aplusi, so both sides of the cell
						Adiag += scale;
					}

					//y
					if (j > 0 && grid->P->getCellMark(i, j - 1, k) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i, j - 1, k, -scale, coefficients);
					}
					if (j < n && grid->P->getCellMark(i, j + 1, k) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i, j + 1, k, -scale, coefficients);
					}
					if (j > 0 && grid->P->getCellMark(i, j - 1, k) == AIR) {
						Adiag += scale;
					}
					if (j < n && grid->P->getCellMark(i, j + 1, k) == AIR) {
						Adiag += scale;
					}

					//z
					if (k > 0 && grid->P->getCellMark(i, j, k - 1) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i, j, k - 1, -scale, coefficients);
					}
					if (k < n && grid->P->getCellMark(i, j, k + 1) == FLUID) {
						Adiag += scale;
						insertCoefficient(id, i, j, k + 1, -scale, coefficients);
					}
					if (k > 0 && grid->P->getCellMark(i, j, k - 1) == AIR) {
						Adiag += scale;
					}
					if (k < n && grid->P->getCellMark(i, j, k + 1) == AIR) {
						Adiag += scale;
					}
				}
				insertCoefficient(id, i, j, k, Adiag, coefficients);
			}
		}
	}
}

void FluidSolver::buildDivergences(Eigen::VectorXd& rhs) {
	float scale = 1 / this->cellSize;
	double divergence = 0.f;
	for (int k = 0; k < resolution.z; ++k) {
		for (int j = 0; j < resolution.y; ++j) {
			for (int i = 0; i < resolution.x; ++i) {
				if (grid->P->getCellMark(i, j, k) == FLUID) {
					int id = k * resolution.x * resolution.y + j * resolution.x + i;
					divergence = scale * (
						(*grid->vel_U)(i + 1, j, k) -
						(*grid->vel_U)(i, j, k)
						+
						(*grid->vel_V)(i, j + 1, k) -
						(*grid->vel_V)(i, j, k)
						+
						(*grid->vel_W)(i, j, k + 1) -
						(*grid->vel_W)(i, j, k)
						);
					rhs[id] = -divergence;
				}
			}
		}
	}
}

void FluidSolver::fillPressureGrid(Eigen::VectorXd x) {
	int iter = 0; int id = 0;
	for (int k = 0; k < resolution.z; ++k) {
		for (int j = 0; j < resolution.y; ++j) {
			for (int i = 0; i < resolution.x; ++i) {
				if (grid->P->getCellMark(i, j, k) == FLUID) {
					id = k * resolution.x * resolution.y + j * resolution.x + i;
					if (std::isnan(x[id]))
						std::cout << "ERROR with PCG" << std::endl;
					grid->P->setCell(i, j, k, ((float)x[id]));
				}
			}
		}
	}
}

void FluidSolver::ProjectPressure() {
	int x = resolution.x, y = resolution.y, z = resolution.z;
	long m = x * y * z;


	Eigen::VectorXd p(m);
	p.setZero();

	Eigen::VectorXd rhs(m);
	rhs.setZero();
	buildDivergences(rhs);

	std::vector<Eigen::Triplet<double>> coefficients;
	buildMatrixA(coefficients, m);
	Eigen::SparseMatrix<double> A(m, m);
	A.setZero();
	A.setFromTriplets(coefficients.begin(), coefficients.end());


	// solve
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> pcg(A);  // performs a Cholesky factorization of A
	p = pcg.solve(rhs);

	// fill eigen vector into grid.P
	fillPressureGrid(p);
	SubtractPressureGradient();
}

void FluidSolver::SubtractPressureGradient() {
	float dx = cellSize;
	float scale = 1.f / (dx);

	int x = resolution.x;
	int y = resolution.y;
	int z = resolution.z;

	//loop over i,j,k
	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y; ++j) {
			for (int i = 0; i < x; ++i) {
				if (grid->P->getCellMark(i - 1, j, k) == FLUID || grid->P->getCellMark(i, j, k) == FLUID) {
					if (grid->P->getCellMark(i - 1, j, k) == SOLID || grid->P->getCellMark(i, j, k) == SOLID) {
						//set cell value as usolid at i,j,k
						grid->vel_U->setCell(i, j, k, 0.f);
					}
					else {
						grid->vel_U->setCellAdd(i, j, k, -(scale * ((*grid->P)(i, j, k) - (*grid->P)(i - 1, j, k))));
					}
				}
				if (grid->P->getCellMark(i, j - 1, k) == FLUID || grid->P->getCellMark(i, j, k) == FLUID) {
					if (grid->P->getCellMark(i, j - 1, k) == SOLID || grid->P->getCellMark(i, j, k) == SOLID) {
						//set cell value as usolid at i,j,k
						grid->vel_V->setCell(i, j, k, 0.f);
					}
					else {
						grid->vel_V->setCellAdd(i, j, k, -(scale * ((*grid->P)(i, j, k) - (*grid->P)(i, j - 1, k))));
					}
				}
				if (grid->P->getCellMark(i, j, k - 1) == FLUID || grid->P->getCellMark(i, j, k) == FLUID) {
					if (grid->P->getCellMark(i, j, k - 1) == SOLID || grid->P->getCellMark(i, j, k) == SOLID) {
						//set cell value as usolid at i,j,k
						grid->vel_W->setCell(i, j, k, 0.f);
					}
					else {
						grid->vel_W->setCellAdd(i, j, k, -(scale * ((*grid->P)(i, j, k) - (*grid->P)(i, j, k - 1))));
					}
				}
			}
		}
	}
}

vec3 EulerStep(const vec3 pos, const vec3 speed, float time_step) {
	return pos + speed * time_step;
}

vec3 FluidSolver::integratePos(const vec3& pos, const vec3& speed, const float& time_step, bool RK2) {
	vec3 new_pos(0.f);
	if (RK2) {
		//RK2 integration
		vec3 k1 = speed * (pos)*time_step / 2.f;
		vec3 k2 = speed * (pos + k1) * time_step;
		//        new_pos = pos + k2;
		new_pos = pos + speed * time_step;
	}
	else {
		//RK4
		vec3 k1 = speed * (pos)*time_step / 2.f;
		vec3 k2 = speed * (pos + k1) * time_step;
		vec3 k3 = speed * (pos + k2) * time_step;
		vec3 k4 = speed * (pos + k3);
		new_pos = pos + (0.1666666f) * (k1 + 2.f * k2 + 2.f * k3 + k4);
	}

	new_pos.x = glm::clamp(new_pos.x + EPSILON, cellSize * 1.f, cellSize * (resolution.x - 1) - EPSILON);
	new_pos.y = glm::clamp(new_pos.y + EPSILON, cellSize * 1.f, cellSize * (resolution.y - 1) - EPSILON);
	new_pos.z = glm::clamp(new_pos.z + EPSILON, cellSize * 1.f, cellSize * (resolution.z - 1) - EPSILON);

	return new_pos;
}

void FluidSolver::ExtrapolateVelocity() {
	int x = resolution.x;
	int y = resolution.y;
	int z = resolution.z;

	MACGridData extrapolate_grid_U(resolution, containerBounds, cellSize);
	MACGridData extrapolate_grid_V(resolution, containerBounds, cellSize);
	MACGridData extrapolate_grid_W(resolution, containerBounds, cellSize);

	//making three grid structures for all axes, so that we store the (fluid or near fluid) and (solid/air or near solid/(solid or air)) mark
	//this is to extrapolate in to the right cells so that tri lerp does not lerp with an empty cell

	//setting a temporary marker grid to get extrapolated velocities

	//x direction
	for (int k = 1; k < z; ++k) {
		for (int j = 1; j < y; ++j) {
			for (int i = 1; i < x; ++i) {
				//marking fluid or near fluid cells as FLUID
				if ((i < x && grid->P->getCellMark(i, j, k) == FLUID) || (i > 0 && grid->P->getCellMark(i - 1, j, k) == FLUID)) {
					extrapolate_grid_U.setCellMark(i, j, k, FLUID);
				}
				//marking solid for extrapolation
				else if (((i < x && (grid->P->getCellMark(i, j, k) == SOLID)) &&
					(i > 0 && (grid->P->getCellMark(i - 1, j, k) == SOLID))) //check if me solid and near solid
					|| ((i < x && (grid->P->getCellMark(i, j, k) == AIR)) &&
					(i > 0 && (grid->P->getCellMark(i - 1, j, k) == AIR))) // bridson
					) //me air, and near (solid or air)
				{
					extrapolate_grid_U.setCellMark(i, j, k, SOLID);
				}
				if ((j < y && grid->P->getCellMark(i, j, k) == FLUID) || (j > 0 && grid->P->getCellMark(i, j - 1, k) == FLUID)) {
					extrapolate_grid_V.setCellMark(i, j, k, FLUID);
				}
				//marking solid for extrapolation
				else if (((j < y && (grid->P->getCellMark(i, j, k) == SOLID)) &&
					(j > 0 && (grid->P->getCellMark(i, j - 1, k) == SOLID))) //check if me solid and near solid
					|| ((j < y && (grid->P->getCellMark(i, j, k) == AIR)) &&
					//me air, and near (solid or air)
					(j > 0 && (grid->P->getCellMark(i, j - 1, k) == AIR)))) //me air, and near (solid or air)

				{
					extrapolate_grid_V.setCellMark(i, j, k, SOLID);
				}
				if ((k < z && grid->P->getCellMark(i, j, k) == FLUID) || (k > 0 && grid->P->getCellMark(i, j, k - 1) == FLUID)) {
					extrapolate_grid_W.setCellMark(i, j, k, FLUID);
				}
				//marking solid for extrapolation
				else if (((k < z && (grid->P->getCellMark(i, j, k) == SOLID)) &&
					(k > 0 && (grid->P->getCellMark(i, j, k - 1) == SOLID))) //check if me solid and near solid
					|| ((k < z && (grid->P->getCellMark(i, j, k) == AIR)) &&
					//me air, and near (solid or air)
					(k > 0 && (grid->P->getCellMark(i, j, k - 1) == AIR)))) //me air, and near (solid or air)
				{
					extrapolate_grid_W.setCellMark(i, j, k, SOLID);
				}


			}
		}
	}

	//neighborhood of 0, based on distance
	for (int k = 0; k < z + 1; ++k) {
		for (int j = 0; j < y + 1; ++j) {
			for (int i = 0; i < x + 1; ++i) {
				ivec3 q[6] = { ivec3(i - 1,j,k), ivec3(i + 1,j,k),
					ivec3(i,j - 1,k), ivec3(i,j + 1,k),
					ivec3(i,j,k - 1), ivec3(i,j,k + 1) };
				//x faces
				if ((j < y && k < z) && extrapolate_grid_U.getCellMark(i, j, k) == SOLID) {
					unsigned int wsum = 0;
					float sum = 0.0f;
					for (unsigned int neighbor = 0; neighbor < 6; ++neighbor) {
						if (q[neighbor][0] >= 0 && q[neighbor][0] < x + 1 && q[neighbor][1] >= 0 &&
							q[neighbor][1] < y && q[neighbor][2] >= 0 && q[neighbor][2] < z) {

							if (extrapolate_grid_U.getCellMark(q[neighbor]) == FLUID) {
								wsum++;
								sum += (*grid->vel_U)(q[neighbor][0], q[neighbor][1], q[neighbor][2]);
							}
						}
					}
					if (wsum) {
						grid->vel_U->setCell(i, j, k, sum / wsum);
					}
				}
				//y faces
				if ((i < x && k < z) && extrapolate_grid_V.getCellMark(i, j, k) == SOLID) {
					unsigned int wsum = 0;
					float sum = 0.0f;

					for (unsigned int neighbor = 0; neighbor < 6; ++neighbor) {
						if (q[neighbor][0] >= 0 && q[neighbor][0] < x && q[neighbor][1] >= 0 &&
							q[neighbor][1] < y + 1 && q[neighbor][2] >= 0 && q[neighbor][2] < z) {

							if (extrapolate_grid_V.getCellMark(q[neighbor]) == FLUID) {
								wsum++;
								sum += (*grid->vel_V)(q[neighbor][0], q[neighbor][1], q[neighbor][2]);
							}
						}
					}
					if (wsum) {
						grid->vel_V->setCell(i, j, k, sum / wsum);
					}
				}
				//z faces
				if ((j < y && i < x) && extrapolate_grid_W.getCellMark(i, j, k) == SOLID) {
					unsigned int wsum = 0;
					float sum = 0.0f;

					for (unsigned int neighbor = 0; neighbor < 6; ++neighbor) {
						if (q[neighbor][0] >= 0 && q[neighbor][0] < x && q[neighbor][1] >= 0 &&
							q[neighbor][1] < y && q[neighbor][2] >= 0 && q[neighbor][2] < z + 1) {

							if (extrapolate_grid_W.getCellMark(q[neighbor]) == FLUID) {
								wsum++;
								sum += (*grid->vel_W)(q[neighbor][0], q[neighbor][1], q[neighbor][2]);
							}
						}
					}
					if (wsum) {
						grid->vel_W->setCell(i, j, k, sum / wsum);
					}
				}
			}
		}
	}

	// on the boundary there's only 1 FLUID neighbor (+1 in the dimension) - more efficient to set just that neighbor
}


void FluidSolver::setBoundaryVelocitiesToZero() {

	int x = resolution.x;
	int y = resolution.y;
	int z = resolution.z;
	//neighborhood of 0, based on distance

	for (int k = 0; k < z; ++k) {
		for (int j = 0; j < y; ++j) {
			grid->vel_U->setCell(0, j, k, 0.f);
			grid->vel_U->setCell(x, j, k, 0.f);
			grid->vel_U->setCell(1, j, k, 0.f);
			grid->vel_U->setCell(x - 1, j, k, 0.f);
		}
	}
	for (int k = 0; k < z; ++k) {
		for (int i = 0; i < x; ++i) {
			grid->vel_V->setCell(i, 0, k, 0.f);
			grid->vel_V->setCell(i, y, k, 0.f);
			//set both 0 and 1 to zero
			grid->vel_V->setCell(i, 1, k, 0.f);
			grid->vel_V->setCell(i, y - 1, k, 0.f);
		}
	}
	for (int j = 0; j < y; ++j) {
		for (int i = 0; i < x; ++i) {
			grid->vel_W->setCell(i, j, 0, 0.f);
			grid->vel_W->setCell(i, j, z, 0.f);
			grid->vel_W->setCell(i, j, 1, 0.f);
			grid->vel_W->setCell(i, j, z - 1, 0.f);
		}
	}
}

void FluidSolver::storeCurrentGridVelocities() {
	*this->grid->flip_vel_U = *this->grid->vel_U;
	*this->grid->flip_vel_V = *this->grid->vel_V;
	*this->grid->flip_vel_W = *this->grid->vel_W;

}

void FluidSolver::clearGrid() {

	std::fill(grid->save_kernel_wt_U->data.begin(), grid->save_kernel_wt_U->data.end(), 0.f);
	std::fill(grid->save_kernel_wt_V->data.begin(), grid->save_kernel_wt_V->data.end(), 0.f);
	std::fill(grid->save_kernel_wt_W->data.begin(), grid->save_kernel_wt_W->data.end(), 0.f);

	std::fill(grid->vel_U->data.begin(), grid->vel_U->data.end(), 0.f);
	std::fill(grid->vel_V->data.begin(), grid->vel_V->data.end(), 0.f);
	std::fill(grid->vel_W->data.begin(), grid->vel_W->data.end(), 0.f);

	std::fill(grid->flip_vel_U->data.begin(), grid->flip_vel_U->data.end(), 0.f);
	std::fill(grid->flip_vel_V->data.begin(), grid->flip_vel_V->data.end(), 0.f);
	std::fill(grid->flip_vel_W->data.begin(), grid->flip_vel_W->data.end(), 0.f);

	std::fill(grid->P->data.begin(), grid->P->data.end(), 0.f);
	std::fill(grid->P->mData.begin(), grid->P->mData.end(), 0);
}

void FluidSolver::step(const float& dt) {
	delta = dt;
	this->initializeMarkerGrid();

	// Step 3 - Store Particle Velocity at current time step to MACGrid
	this->storeParticleVelocityToGrid();

	// Step 4 - Store a temporary copy of Grid Velocities for FLIP
	this->storeCurrentGridVelocities();

	// Step 5 - Add Body Forces like Gravity to MACGrid
	this->CalculateGravityToCell(delta);

	this->setBoundaryVelocitiesToZero();

	//pressure solve
	this->ProjectPressure();

	this->ExtrapolateVelocity();
	this->setBoundaryVelocitiesToZero();

	// Step  - Calculate new flip & pic velocities for each particle
	this->FlipSolve();
	this->PicSolve();

	// Step - Lerp(FLIPVelocity, PICVelocity, 0.95)
	for (int i = 0; i < this->ParticlesContainer.size(); ++i) {
		this->ParticlesContainer.at(i).speed = (1.f - scene.viscosity) * this->particle_save_pic.at(i).speed +
			scene.viscosity * this->particle_save.at(i).speed;
		this->ParticlesContainer.at(i).pos = this->integratePos(this->ParticlesContainer.at(i).pos,
			this->ParticlesContainer.at(i).speed, delta, scene.rk2);
	}
	this->clearGrid();
}
