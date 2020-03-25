//make
// Created by Niklaus on 12.11.19.
//https://github.com/kbladin/Fluid_Simulation This project helped us a lot to get started, and some function are based from there.
//

#ifndef MAIN_FLUIDSOLVER_H
#define MAIN_FLUIDSOLVER_H

#include "FluidDomain.h"
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "GuiData.h"
#include <algorithm>
#include "ViscositySolver.h"
#include "PressureSolver.h"
#include <random>

class FluidSolverMemoryPool
{
public:
    FluidSolverMemoryPool(int size_x, int size_y);
    FluidSolverMemoryPool(const FluidDomain& fluid_domain);
    ~FluidSolverMemoryPool();

    void AssignValidTmpBuffer();

    // Used for velocity extension
    using MatrixXuc = Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>;
    MatrixXuc valid_u_tmp_front;
    MatrixXuc valid_u_tmp_back;
    MatrixXuc valid_v_tmp_front;
    MatrixXuc valid_v_tmp_back;

    Eigen::MatrixXi validU;
	Eigen::MatrixXi validV;

    // Used for PIC solve
    MatrixXt u_sum;
	MatrixXt v_sum;
	MatrixXt u_weight;
	MatrixXt v_weight;
};


/*
 * This class can calculate fluid steps
 */
class FluidSolver
{
public:
    FluidSolver(const FluidSolverMemoryPool& mem_pool);

    ~FluidSolver();

    void stepPICFLIP(FluidDomain& fluid_domain, doubleT dt);

	int* m_iter;
	float* m_acc;
private:

    void advectVelocityField();
    void applyViscosity(doubleT dt);
    void project(doubleT dt);
    void advectFluidParticles(doubleT dt);

    void updateViscosity();

    void addExternalAcceleration(doubleT a_x, doubleT a_y, doubleT dt);
    void enforceDirichlet();
    void solvePoissonCorrectVelocity(doubleT dt);
    void solvePressure(doubleT dt);
    void applyPressure(MatrixXt &pressureGrid, doubleT dt);
    void extrapolateVelocityField(int n_iterations);

    void transferVelocityToGridSpread(Eigen::MatrixXi &isSetU, Eigen::MatrixXi &isSetV);
    void transferTemperatureParticlesToGrid();
    void transferVelocityToParticlesPIC();
    void transferVelocityToParticlesPICFLIP(doubleT pic_ratio, bool dampBorder);
    void transferTemperatureGridToParticles();
    void reseedParticles();

    doubleT cfl_timestep();

    void temperatureSolve(doubleT dt);

    FluidDomain* m_fluid_domain;
    FluidSolverMemoryPool m_mem_pool;

    doubleT m_CFLConditionNumber = 2.0;

    std::mt19937 rng;
};

#endif //MAIN_FLUIDSOLVER_H
