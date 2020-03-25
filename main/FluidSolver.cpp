//
// Created by Niklaus on 12.11.19.
//https://github.com/kbladin/Fluid_Simulation This project helped us a lot to get started, and some function are based from there.
//
#include "FluidSolver.h"
#include <iostream>

using namespace Eigen;

FluidSolverMemoryPool::FluidSolverMemoryPool(int size_x, int size_y) :
        u_sum(size_x + 1, size_y), v_sum(size_x, size_y + 1),
        u_weight(size_x + 1, size_y), v_weight(size_x, size_y + 1),
        valid_u_tmp_front(size_x + 1, size_y), valid_u_tmp_back(size_x + 1, size_y),
        valid_v_tmp_front(size_x, size_y + 1), valid_v_tmp_back(size_x, size_y + 1),
        validU(size_x + 1, size_y), validV(size_x, size_y + 1)
{}

FluidSolverMemoryPool::FluidSolverMemoryPool(const FluidDomain& fluid_domain) :
        u_sum(fluid_domain.getSizeX() + 1, fluid_domain.getSizeY()),
        v_sum(fluid_domain.getSizeX(), fluid_domain.getSizeY() + 1),
        u_weight(fluid_domain.getSizeX() + 1, fluid_domain.getSizeY()),
        v_weight(fluid_domain.getSizeX(), fluid_domain.getSizeY() + 1),
        valid_u_tmp_front(fluid_domain.getSizeX() + 1, fluid_domain.getSizeY()), valid_u_tmp_back(fluid_domain.getSizeX() + 1, fluid_domain.getSizeY()),
        valid_v_tmp_front(fluid_domain.getSizeX(), fluid_domain.getSizeY() + 1), valid_v_tmp_back(fluid_domain.getSizeX(), fluid_domain.getSizeY() + 1),
        validU(fluid_domain.getSizeX() + 1, fluid_domain.getSizeY()), validV(fluid_domain.getSizeX(), fluid_domain.getSizeY() + 1)
{}

FluidSolverMemoryPool::~FluidSolverMemoryPool() = default;

void FluidSolverMemoryPool::AssignValidTmpBuffer()
{
    valid_u_tmp_front = valid_u_tmp_back;
    valid_v_tmp_front = valid_v_tmp_back;
}


FluidSolver::FluidSolver(const FluidSolverMemoryPool& mem_pool)
        : m_mem_pool(mem_pool), rng(124) {}

FluidSolver::~FluidSolver() = default;

void FluidSolver::stepPICFLIP(FluidDomain& fluid_domain, doubleT dt)
{

    m_fluid_domain = &fluid_domain;

    doubleT t = 0.0;
    while (t < dt)
    {
        //ensures that a particle moves at most m_CFLConditionNumber(=5) of cells per iteration
        doubleT substep = cfl_timestep();

        //enforce stability condition with viscosity present -> here you see why we need implicity
        //substep = std::min(substep, pow(fluid_domain.getDeltaX(),2) / (4.0 * viscosity) - 1e-6);
        if (t + substep > dt)
            substep = dt - t;

        std::cout << "substep=" << substep << std::endl;

        // Classify the cells of the domain (AIR, LIQUID or SOLID)
        fluid_domain.classifyCells();

        advectVelocityField();

        // Add gravity
        addExternalAcceleration(0, -9.82, substep);

        //Diffusion
        applyViscosity(substep);

        //Pressure solve
        project(substep);

	    transferTemperatureGridToParticles();

	    advectFluidParticles(substep);

        if(GuiData::instance->m_reseedingOn)
            reseedParticles();

	    transferTemperatureParticlesToGrid();
	    temperatureSolve(dt);

        t += substep;
    }
}

void FluidSolver::advectVelocityField()
{
	Eigen::MatrixXi isValueSetU = Eigen::MatrixXi::Zero(m_fluid_domain->getSizeX() + 1, m_fluid_domain->getSizeY());
	Eigen::MatrixXi isValueSetV = Eigen::MatrixXi::Zero(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY() + 1);

    transferVelocityToGridSpread(isValueSetU, isValueSetV);

    m_mem_pool.validU.setZero();
    m_mem_pool.validV.setZero();

    //Advect U
    for (int j = 0; j < m_fluid_domain->getSizeY(); ++j) {
        for (int i = 0; i < m_fluid_domain->getSizeX() + 1; ++i) {
            bool isEdge = i == 0 || i == m_fluid_domain->getSizeX();
            if(	!isEdge && (m_fluid_domain->cellType(i, j) == LIQUID || m_fluid_domain->cellType(i - 1, j) == LIQUID)
                && isValueSetU(i,j))
            {
                m_mem_pool.validU(i,j) = 1;
                m_fluid_domain->set_uijHalfInd_tmp(i, j, m_fluid_domain->uijHalfInd(i,j));
            } else {
                m_fluid_domain->set_uijHalfInd_tmp(i, j, 0);
            }
        }
    }
    //Advect V
    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; ++j) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); ++i) {
            bool isEdge = j == 0 || j == m_fluid_domain->getSizeY();
            if(!isEdge && (m_fluid_domain->cellType(i, j) == LIQUID || m_fluid_domain->cellType(i, j - 1) == LIQUID)
                && isValueSetV(i,j))
            {
                m_mem_pool.validV(i,j) = 1;
                m_fluid_domain->set_vijHalfInd_tmp(i, j, m_fluid_domain->vijHalfInd(i, j));
            } else {
                m_fluid_domain->set_vijHalfInd_tmp(i, j, 0);
            }
        }
    }
    m_fluid_domain->assignVelocityTmp();
    extrapolateVelocityField((int)std::ceil(m_CFLConditionNumber)+2);
    m_fluid_domain->updatePreviousVelocityBuffer();
}

void FluidSolver::applyViscosity(doubleT dt)
{
    updateViscosity(); //set cellwise viscosity values
    ViscositySolver vis_solver = ViscositySolver();
    vis_solver.applyViscosityToVelocityField(m_fluid_domain, dt);
}


void FluidSolver::project(doubleT dt)
{
    enforceDirichlet();
    solvePoissonCorrectVelocity(dt);
    extrapolateVelocityField((int)std::ceil(m_CFLConditionNumber)+2);
    enforceDirichlet();
}

/**
 * Advect the fluid particles.
 * There is an issue that the velocity field close to the vertical boundary doens't come to rest.
 * To avoid this issue a proposed "hack" was to set pic_ratio to 1 at the boundary. But this overdampens the simulation
 * As a fix, the method calculates the max velocity of all particles, and only if it is low enough, the boundary has adjusted
 * pic-ratio. As an effect the simulation comes succesfully to rest.
 */
void FluidSolver::advectFluidParticles(doubleT dt)
{
    //get diff field
    m_fluid_domain->updateVelocityDiffBuffer();
    // Transfer to particles
    doubleT max_vel = 0.0;
    for (auto & p : m_fluid_domain->getParticles())
    {
        doubleT curr_vel = sqrt(p.velX() * p.velX() + p.velY() * p.velY());
        max_vel = std::max(curr_vel, max_vel);
    }
    bool dampBorder = max_vel < m_fluid_domain->getDeltaX();
    transferVelocityToParticlesPICFLIP(0.01, dampBorder);
    // Advect
    m_fluid_domain->advectAndEnsureOutsideObstacles(dt);
}

/**
 * Sets the viscosity of the cell based on temperature.
 */
void FluidSolver::updateViscosity()
{
    doubleT max_viscosity = GuiData::instance->m_viscosity;
    doubleT maxT = 80;
    doubleT minT = 50;
    doubleT min_viscosity = 5;
    for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {

            doubleT t = m_fluid_domain->Tij(i,j);
            doubleT v;
            if (t < minT)
                v = max_viscosity;
            else if (t > maxT)
                v = min_viscosity;
            else {
                t -= minT;
                t *= (GuiData::instance->m_brushTemperature / (maxT - minT));
                doubleT factor = 1.0 / std::max(pow((t), 2), 1.0);
                v = std::max(max_viscosity * factor, min_viscosity);
            }
            m_fluid_domain->setViscosity(i, j, v);
        }
    }
}

/**
 * Transfer particle velocity onto the grid - build velocity field.
 */
void FluidSolver::transferVelocityToGridSpread(Eigen::MatrixXi &isSetU, Eigen::MatrixXi &isSetV)
{
    if (!m_fluid_domain)
        return;

    m_mem_pool.u_sum.setZero();
    m_mem_pool.v_sum.setZero();
    m_mem_pool.u_weight.setZero();
    m_mem_pool.v_weight.setZero();

    for (auto & p : m_fluid_domain->getParticles())
    {
        m_fluid_domain->addToValueInterpolated(p.posX(), p.posY() - 0.5 * m_fluid_domain->getDeltaY(), p.velX(), m_mem_pool.u_sum);
        m_fluid_domain->addToValueInterpolated(p.posX() - 0.5 * m_fluid_domain->getDeltaX(), p.posY(), p.velY(), m_mem_pool.v_sum);
        m_fluid_domain->addToValueInterpolated(p.posX(), p.posY() - 0.5 * m_fluid_domain->getDeltaY(), 1.0, m_mem_pool.u_weight);
        m_fluid_domain->addToValueInterpolated(p.posX() - 0.5 * m_fluid_domain->getDeltaX(), p.posY(), 1.0, m_mem_pool.v_weight);
    }

    doubleT eps = 1e-9;
    //Construct U velocity
    for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX() + 1; i++) {
            if (m_mem_pool.u_weight(i, j) < eps)
                continue;
            doubleT new_u = m_mem_pool.u_sum(i, j) / m_mem_pool.u_weight(i, j);
            m_fluid_domain->set_uijHalfInd_tmp(i, j, new_u);
            isSetU(i,j) = 1;
        }
    }
    //Construct V velocity
    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {
            if (m_mem_pool.v_weight(i, j) < eps)
                continue;
            doubleT new_v = m_mem_pool.v_sum(i, j) / m_mem_pool.v_weight(i, j);
            m_fluid_domain->set_vijHalfInd_tmp(i, j, new_v);
            isSetV(i,j) = 1;
        }
    }
    m_fluid_domain->assignVelocityTmp();
}

/**
 * Iteratively adds velocity to non-fluid cells that have fluid neighbors. Therefore extending the velocity field at the boundary
 *
 */
void FluidSolver::extrapolateVelocityField(int n_iterations)
{
    if (!m_fluid_domain)
        return;

   char UNKNOWN = 0x00;
   char WAITING = 0x01;
   char KNOWN = 0x02;
   char DONE = 0x03;

   using MatrixXc = Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic>;
   MatrixXc statusU(m_fluid_domain->getSizeX() + 1, m_fluid_domain->getSizeY());
   for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
       for (int i = 0; i < m_fluid_domain->getSizeX() + 1; i++) {
           statusU(i,j) = m_mem_pool.validU(i,j) ? KNOWN : UNKNOWN;
           bool isEdge = i == 0 || i == m_fluid_domain->getSizeX() || j == 0 || j == m_fluid_domain->getSizeY() - 1;
           if (isEdge && statusU(i,j) == UNKNOWN)
               statusU(i,j) = DONE;
       }
   }

    MatrixXc statusV(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY() + 1);
    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {
            statusV(i,j) = m_mem_pool.validV(i,j) ? KNOWN : UNKNOWN;
            bool isEdge = i == 0 || i == m_fluid_domain->getSizeX() - 1 || j == 0 || j == m_fluid_domain->getSizeY();
            if (isEdge && statusV(i,j) == UNKNOWN)
                statusV(i,j) = DONE;
        }
    }

   struct GridIndex {
       int i; int j;
       GridIndex() = default;
       GridIndex(int _i, int _j) : i(_i), j(_j) {}
   };

   std::vector<GridIndex> extrapolationCellsU;
   std::vector<GridIndex> extrapolationCellsV;
   for (int iter = 0; iter < n_iterations; iter++)
   {
       extrapolationCellsU.clear();
       for (int j = 1; j < m_fluid_domain->getSizeY() - 1; j++) {
           for (int i = 1; i < m_fluid_domain->getSizeX(); i++) {
               if (statusU(i,j) != KNOWN)
                   continue;

               int count = 0;
               if (statusU(i - 1, j) == UNKNOWN) {
                   extrapolationCellsU.emplace_back(i - 1, j);
                   statusU(i - 1, j) = WAITING;
                   count++;
               } else if (statusU(i - 1, j) == WAITING) {
                   count++;
               }

               if (statusU(i + 1, j) == UNKNOWN) {
                   extrapolationCellsU.emplace_back(i + 1, j);
                   statusU(i + 1, j) = WAITING;
                   count++;
               } else if (statusU(i + 1, j) == WAITING) {
                   count++;
               }

               if (statusU(i, j - 1) == UNKNOWN) {
                   extrapolationCellsU.emplace_back(i, j - 1);
                   statusU(i, j - 1) = WAITING;
                   count++;
               } else if (statusU(i, j - 1) == WAITING) {
                   count++;
               }

               if (statusU(i, j + 1) == UNKNOWN) {
                   extrapolationCellsU.emplace_back(i, j + 1);
                   statusU(i, j + 1) = WAITING;
                   count++;
               } else if (statusU(i, j + 1) == WAITING) {
                   count++;
               }

               if (count == 0)
                   statusU(i, j) = DONE;
           }
       }

       extrapolationCellsV.clear();
       for (int j = 1; j < m_fluid_domain->getSizeY(); j++) {
           for (int i = 1; i < m_fluid_domain->getSizeX() - 1; i++) {
               if (statusV(i,j) != KNOWN)
                   continue;

               int count = 0;
               if (statusV(i - 1, j) == UNKNOWN) {
                   extrapolationCellsV.emplace_back(i - 1, j);
                   statusV(i - 1, j) = WAITING;
                   count++;
               } else if (statusV(i - 1, j) == WAITING) {
                   count++;
               }

               if (statusV(i + 1, j) == UNKNOWN) {
                   extrapolationCellsV.emplace_back(i + 1, j);
                   statusV(i + 1, j) = WAITING;
                   count++;
               } else if (statusV(i + 1, j) == WAITING) {
                   count++;
               }

               if (statusV(i, j - 1) == UNKNOWN) {
                   extrapolationCellsV.emplace_back(i, j - 1);
                   statusV(i, j - 1) = WAITING;
                   count++;
               } else if (statusV(i, j - 1) == WAITING) {
                   count++;
               }

               if (statusV(i, j + 1) == UNKNOWN) {
                   extrapolationCellsV.emplace_back(i, j + 1);
                   statusV(i, j + 1) = WAITING;
                   count++;
               } else if (statusV(i, j + 1) == WAITING) {
                   count++;
               }

               if (count == 0)
                   statusV(i, j) = DONE;
           }
       }

       for (auto & g : extrapolationCellsU)
       {
           doubleT sum = 0.0; int count = 0;
           if (statusU(g.i - 1, g.j) == KNOWN) { sum += m_fluid_domain->uijHalfInd(g.i - 1, g.j); count++; }
           if (statusU(g.i + 1, g.j) == KNOWN) { sum += m_fluid_domain->uijHalfInd(g.i + 1, g.j); count++; }
           if (statusU(g.i, g.j - 1) == KNOWN) { sum += m_fluid_domain->uijHalfInd(g.i, g.j - 1); count++; }
           if (statusU(g.i, g.j + 1) == KNOWN) { sum += m_fluid_domain->uijHalfInd(g.i, g.j + 1); count++; }
           m_fluid_domain->set_uijHalfInd(g.i, g.j, sum / count);
       }
       for (auto & g : extrapolationCellsU)
           statusU(g.i, g.j) = KNOWN;

       for (auto & g : extrapolationCellsV)
       {
           doubleT sum = 0.0; int count = 0;
           if (statusV(g.i - 1, g.j) == KNOWN) { sum += m_fluid_domain->vijHalfInd(g.i - 1, g.j); count++; }
           if (statusV(g.i + 1, g.j) == KNOWN) { sum += m_fluid_domain->vijHalfInd(g.i + 1, g.j); count++; }
           if (statusV(g.i, g.j - 1) == KNOWN) { sum += m_fluid_domain->vijHalfInd(g.i, g.j - 1); count++; }
           if (statusV(g.i, g.j + 1) == KNOWN) { sum += m_fluid_domain->vijHalfInd(g.i, g.j + 1); count++; }
           m_fluid_domain->set_vijHalfInd(g.i, g.j, sum / count);
       }
       for (auto & g : extrapolationCellsV)
           statusV(g.i, g.j) = KNOWN;
   }
}

/**
 * Adds forces. e.g gravity
 */
void FluidSolver::addExternalAcceleration(doubleT a_x, doubleT a_y, doubleT dt)
{
    if (!m_fluid_domain)
        return;

    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; ++j)
    {
        for (int i = 0; i < m_fluid_domain->getSizeX(); ++i)
        {
            bool isEdge = j == 0 || j == m_fluid_domain->getSizeY();
            if (!isEdge && !(m_fluid_domain->isBorder(i, j) || m_fluid_domain->isBorder(i, j - 1)) && (m_fluid_domain->phi(i, j) < 0.0 || m_fluid_domain->phi(i, j - 1) < 0.0))
            { // Only add force to the liquid cells
                m_fluid_domain->set_vijHalfInd(i,j, m_fluid_domain->vijHalfInd(i,j) + a_y * dt);
            }
        }
    }
}

/**
 * Enforce free-slip boundary condition. E.g. velocities orthogonal to boundary are set to 0.
 */
void FluidSolver::enforceDirichlet()
{
    if (!m_fluid_domain)
        return;

    //U Velocity
    for (int j = 0; j < m_fluid_domain->getSizeY(); j++)
    {
        for (int i = 0; i < m_fluid_domain->getSizeX() + 1; i++)
        {
            bool isEdge = i == 0 || i == m_fluid_domain->getSizeX();
            // X velocity
            if(isEdge || (m_fluid_domain->isBorder(i - 1, j) && m_fluid_domain->uijHalfInd(i, j) < 0) ||
                (m_fluid_domain->isBorder(i, j) && m_fluid_domain->uijHalfInd(i, j) > 0))
            {
                m_fluid_domain->set_uijHalfInd(i, j, 0);
                m_fluid_domain->set_uijHalfInd_tmp(i, j, 0);
            }

        }
    }

    //V Velocity
    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; j++)
    {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++)
        {
            bool isEdge = j == 0 || j == m_fluid_domain->getSizeY();
            if(isEdge || (m_fluid_domain->isBorder(i, j - 1) && m_fluid_domain->vijHalfInd(i, j) < 0) ||
                   (m_fluid_domain->isBorder(i, j) && m_fluid_domain->vijHalfInd(i, j) > 0))
            {
                m_fluid_domain->set_vijHalfInd(i, j, 0);
                m_fluid_domain->set_vijHalfInd_tmp(i, j, 0);
            }

        }
    }
}

/**
 * Solves pressure with Gauss-Seidel and applies it. Assumes dirichlet boundary set
 * Ghost pressures are calculated according to liquid SDF at air/fluid boundary
 */
void FluidSolver::solvePoissonCorrectVelocity(doubleT dt)
{
    if (!m_fluid_domain)
        return;

    int it = 0;
    doubleT residual = GuiData::instance->m_acc + 1;

    ParticleLevelSet& SDF = m_fluid_domain->getLiquidSDF();

    m_fluid_domain->clearPressure();
    for (it = 0; residual > GuiData::instance->m_acc && it < GuiData::instance->m_iter; ++it) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {
            for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {

                if (m_fluid_domain->isBorder(i,j) || SDF(i,j) >= 0)
                {
                    m_fluid_domain->set_pij(i, j, 0);
                    continue;
                }

                doubleT b = -(m_fluid_domain->uij_div(i, j) + m_fluid_domain->vij_div(i, j)) / dt; // right-hand
                doubleT Aij = 0.0; //left_hand
                int n_non_solid_neighbors = 0;

                if (i - 1 >= 0 && !m_fluid_domain->isBorder(i - 1,j)) {
                    if (SDF(i - 1, j) < 0.0) {
                        Aij += m_fluid_domain->pij(i - 1, j);
                    }
                    else {
                        Aij -= SDF.getEdgeWeightU(i, j) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }

                if (i + 1 < m_fluid_domain->getSizeX() && !m_fluid_domain->isBorder(i + 1,j)) {
                    if (SDF(i + 1, j) < 0.0) {
                        Aij += m_fluid_domain->pij(i + 1, j);
                    } else {
                        Aij -= SDF.getEdgeWeightU(i + 1, j) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }
                if (j - 1 >= 0 && !m_fluid_domain->isBorder(i, j - 1)) {
                    if (SDF(i, j - 1) < 0.0) {
                        Aij += m_fluid_domain->pij(i, j - 1);
                    } else {
                        Aij -= SDF.getEdgeWeightV(i, j) * m_fluid_domain->pij(i, j);
                    }
                    n_non_solid_neighbors++;
                }
                if (j + 1 < m_fluid_domain->getSizeY() && !m_fluid_domain->isBorder(i, j + 1)) {
                    if (SDF(i, j + 1) < 0.0) {
                        Aij += m_fluid_domain->pij(i,j + 1);
                    } else {
                        Aij -= SDF.getEdgeWeightV(i, j + 1) * m_fluid_domain->pij(i, j);
                    }
                    n_non_solid_neighbors++;
                }
                if (n_non_solid_neighbors == 0)
                    m_fluid_domain->set_pij(i, j,0);
                else
                {
                    doubleT new_p = (pow(m_fluid_domain->getDeltaX(),2) * b + Aij) / n_non_solid_neighbors;
                    m_fluid_domain->set_pij(i, j, new_p);
                }
            }
        }
        // Compute the new residual, i.e. the sum of the squares of the individual residuals (squared L2-norm)
        ///Residual computation commented, bc we use a fixed number of iterations. By experience the desired accuracy is never met anyway
        /*
         residual = 0;
         for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
            for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {

                if (m_fluid_domain->isBorder(i,j) || SDF(i,j) >= 0)
                    continue;

                doubleT b = -(m_fluid_domain->uij_div(i, j) + m_fluid_domain->vij_div(i, j)) / dt; // right-hand
                doubleT Aij = 0.0; //left_hand
                int n_non_solid_neighbors = 0;
                if (i - 1 >= 0 && !m_fluid_domain->isBorder(i - 1,j)) {
                    if (SDF(i - 1, j) < 0.0) {
                        Aij -= m_fluid_domain->pij(i - 1, j);
                    } else {
                        Aij += SDF.getEdgeWeightU(i, j) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }
                if (i + 1 < m_fluid_domain->getSizeX() && !m_fluid_domain->isBorder(i + 1,j)) {
                    if (SDF(i + 1, j) < 0.0) {
                        Aij -= m_fluid_domain->pij(i + 1, j);
                    } else {
                        Aij += SDF.getEdgeWeightU(i + 1, j) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }
                if (j - 1 >= 0 && !m_fluid_domain->isBorder(i,j - 1)) {
                    if (SDF(i, j - 1) < 0.0) {
                        Aij -= m_fluid_domain->pij(i, j - 1);
                    } else {
                        Aij += SDF.getEdgeWeightV(i, j) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }
                if (j + 1 < m_fluid_domain->getSizeY() && !m_fluid_domain->isBorder(i, j + 1)) {
                    if (SDF(i, j + 1) < 0.0) {
                        Aij -= m_fluid_domain->pij(i, j + 1);
                    } else {
                        Aij += SDF.getEdgeWeightV(i, j + 1) * m_fluid_domain->pij(i,j);
                    }
                    n_non_solid_neighbors++;
                }
                doubleT cellResidual = b - (n_non_solid_neighbors * m_fluid_domain->pij(i, j) - Aij) / pow(m_fluid_domain->getDeltaX(), 2);

                residual += cellResidual * cellResidual;
            }

        }

        // Get the L2-norm of the residual
        residual = sqrt(residual);

        // We assume the accuracy is meant for the average L2-norm per grid cell
        residual /= (m_fluid_domain->getSizeX() - 2) * (m_fluid_domain->getSizeY() - 2);
         */
    }
    std::cout << "\t Pressure iterations: " << it<<
              "\n\t Pressure error: " << residual << "\n\n";

    //Correct velocity U
    for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX() + 1; i++) {
            bool isEdge = i == 0 || i == m_fluid_domain->getSizeX();
            if (!isEdge) {
                doubleT p0 = m_fluid_domain->pij(i - 1, j);
                doubleT p1 = m_fluid_domain->pij(i,j);
                if (!m_fluid_domain->isBorder(i-1,j) && SDF(i-1,j) >= 0.0 && !m_fluid_domain->isBorder(i,j) && SDF(i,j) < 0.0)
                    p0 = -SDF.getEdgeWeightU(i,j) * p1;
                if (!m_fluid_domain->isBorder(i,j) && SDF(i,j) >= 0.0 && !m_fluid_domain->isBorder(i - 1,j) && SDF(i - 1,j) < 0.0)
                    p1 = -SDF.getEdgeWeightU(i,j) * p0;
                doubleT u_new = m_fluid_domain->uijHalfInd(i, j) - dt * (p1 - p0) / (m_fluid_domain->getDeltaX());
                m_fluid_domain->set_uijHalfInd_tmp(i, j, u_new);
            } else
                m_fluid_domain->set_uijHalfInd_tmp(i, j, 0.0);
        }
    }
    //Correct velocity V
    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {
            bool isEdge = j == 0 || j == m_fluid_domain->getSizeY();
            if (!isEdge) {
                doubleT p0 = m_fluid_domain->pij(i, j - 1);
                doubleT p1 = m_fluid_domain->pij(i,j);
                if (!m_fluid_domain->isBorder(i,j-1) && SDF(i,j-1) >= 0.0 && !m_fluid_domain->isBorder(i,j) && SDF(i,j) < 0.0)
                    p0 = -SDF.getEdgeWeightV(i,j) * p1;
                if (!m_fluid_domain->isBorder(i,j) && SDF(i,j) >= 0.0 && !m_fluid_domain->isBorder(i ,j -1) && SDF(i ,j -1) < 0.0)
                    p1 = -SDF.getEdgeWeightV(i,j) * p0;
                doubleT v_new = m_fluid_domain->vijHalfInd(i, j) - dt * (p1 - p0) / (m_fluid_domain->getDeltaY());
                m_fluid_domain->set_vijHalfInd_tmp(i, j, v_new);
            } else
                m_fluid_domain->set_vijHalfInd_tmp(i, j, 0.0);
        }
    }
    m_fluid_domain->assignVelocityTmp();
}

//Sets particle velocity based on grid velocity for PIC step.
void FluidSolver::transferVelocityToParticlesPIC()
{
    for (auto p = m_fluid_domain->getParticles().begin(); p != m_fluid_domain->getParticles().end(); p++)
    {
        p->setVelocity(m_fluid_domain->u_interpolate(p->posX(), p->posY()), m_fluid_domain->v_interpolate(p->posX(), p->posY()));
    }
}

/**
 * @param pic_ratio
 * @param dampBorder "hack" to avoid ongoing velocity at boundary.
 */
void FluidSolver::transferVelocityToParticlesPICFLIP(doubleT pic_ratio, bool dampBorder)
{
    if (!m_fluid_domain)
        return;

    for (auto & p : m_fluid_domain->getParticles())
    {
        if (dampBorder) {
            int i = m_fluid_domain->doubleT_to_floor(p.posX() / m_fluid_domain->getDeltaX());
            int j = m_fluid_domain->doubleT_to_floor(p.posY() / m_fluid_domain->getDeltaY());
            i = m_fluid_domain->CLAMP(i, 1, m_fluid_domain->getSizeX() - 2);
            j = m_fluid_domain->CLAMP(j, 1, m_fluid_domain->getSizeY() - 2);

            if (m_fluid_domain->cellType(i - 1, j) == SOLID || m_fluid_domain->cellType(i + 1, j) == SOLID) {
                pic_ratio = 1.0;
            }
        }

        doubleT pic_u = m_fluid_domain->u_interpolate(p.posX(), p.posY());
        doubleT pic_v = m_fluid_domain->v_interpolate(p.posX(), p.posY());

        doubleT flip_u = p.velX() + m_fluid_domain->u_diff_interpolate(p.posX(), p.posY());
        doubleT flip_v = p.velY() + m_fluid_domain->v_diff_interpolate(p.posX(), p.posY());

        p.setVelocity(pic_u * pic_ratio + flip_u * (1.0 - pic_ratio), pic_v * pic_ratio + flip_v * (1.0 - pic_ratio));
    }
}

/**
 * @return maximal timestep to stay within velocity field
 */
doubleT FluidSolver::cfl_timestep()
{
    if (!m_fluid_domain)
        return 0.1;

    doubleT max_velocity = 0.0;

    for (int j = 0; j < m_fluid_domain->getSizeY() + 1; j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX(); i++) {
            max_velocity = std::max(max_velocity, m_fluid_domain->vijHalfInd(i, j));
        }
    }

    for (int j = 0; j < m_fluid_domain->getSizeY(); j++) {
        for (int i = 0; i < m_fluid_domain->getSizeX() + 1; i++) {
            max_velocity = std::max(max_velocity, m_fluid_domain->uijHalfInd(i, j));
        }
    }
    return (m_CFLConditionNumber * m_fluid_domain->getDeltaX()) / max_velocity;
}

//TODO: clamp particle temp to min/max of nearby particles
//problem is that there is a small error and it accumulates
void FluidSolver::transferTemperatureGridToParticles()
{
	for (auto & p : m_fluid_domain->getParticles())
	{
		p.T = m_fluid_domain->valueInterpolated(p.posX(), p.posY(), m_fluid_domain->getTemperature());
		m_fluid_domain->addToValueInterpolated(p.posX(), p.posY(), -p.T * GuiData::instance->m_particleTemperatureTransfer, m_fluid_domain->getTemperature());
	}

}

void FluidSolver::transferTemperatureParticlesToGrid()
{
	for (auto & p : m_fluid_domain->getParticles())
	{
		m_fluid_domain->addToValueInterpolated(p.posX(), p.posY(), p.T * GuiData::instance->m_particleTemperatureTransfer, m_fluid_domain->getTemperature());
	}

	// Last minute note 18-12-19: tried to clamp the values to the mix/max of the *grid* neighborhood,
	// but gives strange effect that temperature moves to the left, not sure why
	// This clamping would only reduce errors that create a local extrema and wouldn't necessarily conserve energy like the current method does.
//	const auto& Clamp = [](doubleT& v, doubleT lo, doubleT hi){
//		return std::max(lo, std::min(hi, v));
//	};
//
//	MatrixXt tmp = m_fluid_domain->getTemperature();
//	for (int j = 1; j < m_fluid_domain->getSizeY() -1; j++) {
//		for (int i = 1; i < m_fluid_domain->getSizeX() -1; i++) {
//			doubleT T = m_fluid_domain->Tij(i, j);
//			std::vector<doubleT> neighbors = {tmp(i +1, j), tmp(i -1, j), tmp(i, j +1), tmp(i, j -1)};
//
//			doubleT min = *std::min(neighbors.begin(), neighbors.end());
//			doubleT max = *std::max(neighbors.begin(), neighbors.end());
//
//			m_fluid_domain->set_Tij(i, j, Clamp(T, min, max));
//		}
//	}
}

/**
 * Attempts to solve Poisson for heat diffusion and timestep T field
 * Using approach from paper Melting and Flowing (melf.pdf), 
 * set up an LSE in each dimension (operator splitting) and solve it with conjugate gradient as its SPD
 */
void FluidSolver::temperatureSolve(doubleT dt)
{
    if (!m_fluid_domain)
        return;

//	auto start = std::chrono::high_resolution_clock::now();
	MatrixXt& TMat = m_fluid_domain->getTemperature();
	auto T = VectorXt::Map(TMat.data(), TMat.cols() * TMat.rows()); //raw reference to std::vector<d> as an Eigen VectorXd
//	std::cout << T.transpose() << std::endl;

	const auto& solver = m_fluid_domain->getTemperatureSolver();
	T = solver.solveWithGuess(T, T);

//	auto end = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<doubleT> diff = end-start;
//	std::cout << "Ttime: " << diff.count() << " s\n";

//	std::cout << T.transpose() << std::endl;
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
}

/**
 * Took the smart reseeding algorithm described in https://pdfs.semanticscholar.org/a1bb/ba8ad75b4ffdaebfe56ce1aec35414247d14.pdf
 * Check that each cell has at most 12 particles (remove other) and at least 3 particles (add 1 particle if too few)
 * Adapted 2 things to prevent fluid expansion: cells at boundary and cells with too high velocity
 * dont receive new particles.
 * The velocity of the added particle is copied from a random particle in the same cell.
 */
void FluidSolver::reseedParticles() {
    auto & particles = m_fluid_domain->getParticles();
	MatrixXi num_particles = MatrixXi::Constant(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY(), 0);
    MatrixXt minimaX = MatrixXt::Constant(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY(), m_fluid_domain->getDeltaX());
	MatrixXt maximaX = MatrixXt::Zero(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY());
	MatrixXt minimaY = MatrixXt::Constant(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY(), m_fluid_domain->getDeltaY());
	MatrixXt maximaY = MatrixXt::Zero(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY());

	MatrixXi indices = MatrixXi::Constant(m_fluid_domain->getSizeX(), m_fluid_domain->getSizeY(), -1);
    for (int k = 0; k < particles.size(); k++)
    {
        MarkerParticle& p = particles[k];
        // Find the particles indices in the grid
        int i = m_fluid_domain->doubleT_to_floor(p.posX() / m_fluid_domain->getDeltaX());
        int j = m_fluid_domain->doubleT_to_floor(p.posY() / m_fluid_domain->getDeltaY());
        i = m_fluid_domain->CLAMP(i, 0, m_fluid_domain->getSizeX() - 1);
        j = m_fluid_domain->CLAMP(j, 0, m_fluid_domain->getSizeY() - 1);
        num_particles(i, j)++;
        if (indices(i, j) == -1)
            indices(i,j) = k;

        if (num_particles(i, j) > 12)
        {
            p.markToDelete();
            num_particles(i,j)--;
        }
        doubleT i_frac = p.posX() / m_fluid_domain->getDeltaX() - i;
        doubleT j_frac = p.posY() / m_fluid_domain->getDeltaY() - j;

        minimaX(i,j) = std::min(i_frac, minimaX(i,j));
        maximaX(i,j) = std::max(i_frac, maximaX(i,j));
        minimaY(i,j) = std::min(j_frac, minimaY(i,j));
        maximaY(i,j) = std::max(j_frac, minimaY(i,j));
    }
    //remove all marked particles, such that all cells have at most 12 particles.
    particles.erase( std::remove_if(particles.begin(), particles.end(),
            [](MarkerParticle & p) { return p.isMarkedToDelete(); }), particles.end());

    doubleT dx = m_fluid_domain->getDeltaX();
    for (int j = 1; j < m_fluid_domain->getSizeY() - 1; j++) {
        for (int i = 1; i < m_fluid_domain->getSizeX() - 1; i++) {
            doubleT vel = sqrt(pow(m_fluid_domain->uij(i,j),2) + pow(m_fluid_domain->vij(i,j),2));
            if (m_fluid_domain->cellType(i,j) == LIQUID && num_particles(i,j) < 3 && vel < dx / 2.0)  // add some.
            {
                doubleT boxMinX = m_fluid_domain->cellType(i - 1,j) == LIQUID ? 0.0 : minimaX(i, j);
                doubleT boxMaxX = m_fluid_domain->cellType(i + 1, j) == LIQUID ? dx : maximaX(i, j);
                doubleT boxMinY = m_fluid_domain->cellType(i, j - 1) == LIQUID ? 0.0 : minimaY(i, j);
                doubleT boxMaxY = m_fluid_domain->cellType(i, j + 1) == LIQUID ? dx : maximaY(i, j);
                if (m_fluid_domain->cellType(i-1,j) != LIQUID) continue;
                if (m_fluid_domain->cellType(i+1,j) != LIQUID) continue;
                if (m_fluid_domain->cellType(i,j-1) != LIQUID) continue;
                if (m_fluid_domain->cellType(i,j+1) != LIQUID) continue;

                if (num_particles(i,j) == 1)
                {
                    boxMinX -= 5e-5;
                    boxMaxX += 5e-5;
                    boxMinY -= 5e-5;
                    boxMaxY += 5e-5;
                }
                std::uniform_real_distribution<> distrX(boxMinX, boxMaxX);
                std::uniform_real_distribution<> distrY(boxMinY, boxMaxY);

                doubleT posX = i * dx + distrX(rng);
                doubleT posY = j * dx + distrY(rng);

                if (indices(i,j) != -1) {
                    doubleT vel_u = particles[indices(i,j)].velX(); //m_fluid_domain->u_interpolate(posX, posY);
                    doubleT vel_v = particles[indices(i,j)].velY(); //m_fluid_domain->v_interpolate(posX, posY);
                    MarkerParticle p(posX, posY, vel_u, vel_v);
                    particles.push_back(p);
                } else {
                    std::cout << "error";
                }
            }
        }
    }
}
