/*
    Created by Niklaus on 28.11.19.
    Method is based on https://github.com/rlguy/FLIPViscosity3D/
*/

#include "ViscositySolver.h"
#include <limits>
#include <iostream>

using namespace Eigen;

void ViscositySolver::applyViscosityToVelocityField(FluidDomain* fluid_domain, doubleT dt) {

    m_size_x = fluid_domain->getSizeX();
    m_size_y = fluid_domain->getSizeY();
    m_delta_x = fluid_domain->getDeltaX();
    m_delta_y = fluid_domain->getDeltaY();
    m_dt = dt;
    m_fluid_domain = fluid_domain;

    doubleT minimum = std::numeric_limits<doubleT>::max();

    for (int j = 0; j < m_size_y; j++) {
        for (int i = 0; i < m_size_x; i++) {
            minimum = std::min(minimum, m_fluid_domain->viscosity(i, j));
        }
    }
    if (minimum <= 0)
        return;

    computeEdgeStateGrid(); //Assign each edge a status
    computePlaneGrid(); //precompute stress samples
    computeMatrixIndexTable(); //grid index to matrix index
    int n_rows = m_matrixIndex.matrixSize;

    std::vector<Eigen::Triplet<doubleT>> triplets;
    triplets.reserve(5 * n_rows);
    VectorXt rhs(n_rows); //holds current velocities plus boundary values.

    Eigen::SparseMatrix<doubleT> A;
    A.resize(n_rows, n_rows);
    A.reserve(5 * n_rows);

    initializeLinearSystem(rhs, triplets); //construct A* U_new = U (U is rhs)
    planes.destroy(); // free some memory.

    //Set up system AU_new = U
    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<doubleT>, Eigen::Lower|Eigen::Upper> m_cg_solver;
    m_cg_solver.compute(A);
    VectorXt U_new(n_rows);

    //Solve it
    U_new = m_cg_solver.solve(rhs);

    std::cout << "\n\t Viscosity iterations: " << m_cg_solver.iterations() <<
                 "\n\t Viscosity error: " << m_cg_solver.error() << "\n\n";

    //Assign updated velocities
    for (int j = 0; j < m_size_y; j++) {
        for (int i = 0; i < m_size_x + 1; i++) {
            int index = m_matrixIndex.U(i, j);
            if (index != -1)
                m_fluid_domain->set_uijHalfInd_tmp(i, j, U_new[index]);
            else
                m_fluid_domain->set_uijHalfInd_tmp(i,j, 0.0);
        }
    }

    for (int j = 0; j < m_size_y + 1; j++) {
        for (int i = 0; i < m_size_x; i++) {
            int index = m_matrixIndex.V(i, j);
            if (index != -1)
                m_fluid_domain->set_vijHalfInd_tmp(i, j, U_new[index]);
            else
                m_fluid_domain->set_vijHalfInd_tmp(i,j, 0.0);
        }
    }
    m_fluid_domain->assignVelocityTmp();
}

/**
 * Classifies all edges to be fluid or solids.
 */
void ViscositySolver::computeEdgeStateGrid()
{
    m_state = EdgeStateGrid(m_size_x, m_size_y);

    for (int j = 0; j < m_state.U.cols(); j++) {
        for (int i = 0; i < m_state.U.rows(); i++)
        {
            bool isEdge = i == 0 || i == m_state.U.rows() - 1;
            if (isEdge || (m_fluid_domain->isBorder(i - 1, j) || m_fluid_domain->isBorder(i,j)))
                m_state.U(i,j) = EdgeState::solid;
            else
                m_state.U(i,j) = EdgeState ::fluid;
        }
    }

    for (int j = 0; j < m_state.V.cols(); j++) {
        for (int i = 0; i < m_state.V.rows(); i++)
        {
            bool isEdge = j == 0 || j == m_state.V.cols() - 1;
            if (isEdge || (m_fluid_domain->isBorder(i, j - 1) || m_fluid_domain->isBorder(i,j)))
                m_state.V(i,j) = EdgeState::solid;
            else
                m_state.V(i,j) = EdgeState ::fluid;
        }
    }
}

/**
 * Determines the cells that will be affected by the viscosity solve
 *  (all the fluid cells + 2 layers of extrapolation around the boundary)
 * Precomputes the stress samples on the mac grid described in https://cs.uwaterloo.ca/~c2batty/papers/BattyBridson08.pdf
 */
void ViscositySolver::computePlaneGrid() {
    planes = ViscosityPlaneGrid(m_size_x, m_size_y);

    MatrixXi validCells(m_size_x + 1, m_size_y + 1);
    validCells.setZero();
    for (int j = 0; j < m_size_y; j++) {
        for (int i = 0; i < m_size_x; i++) {
            if (m_fluid_domain->phi(i,j) < 0.0)
            {
                validCells(i,j) = 1;
            }
        }
    }
    int layers = 2;
    for (int layer = 0; layer < layers; layer++) {
        MatrixXi tempValid = validCells;
        for (int j = 0; j < m_size_y + 1; j++) {
            for (int i = 0; i < m_size_x + 1; i++) {
                if (validCells(i,j)) {
                    if (i + 1 <= m_size_x)  tempValid(i + 1, j) = 1;
                    if (i - 1 >= 0)         tempValid(i - 1, j) = 1;
                    if (j + 1 <= m_size_y)  tempValid(i, j + 1) = 1;
                    if (j - 1 >= 0)         tempValid(i, j - 1) = 1;
                }
            }
        }
        validCells = tempValid;
    }

    doubleT half_dx = 0.5 * m_delta_x;
    estimatePlaneFractions(planes.center, half_dx, half_dx, validCells); //center of a cell
    estimatePlaneFractions(planes.U, 0, half_dx, validCells); // vertical edges of a cell
    estimatePlaneFractions(planes.V, half_dx, 0, validCells); // horizontal edges of a cell
    estimatePlaneFractions(planes.pointU, 0, 0,validCells); // cell corner
    estimatePlaneFractions(planes.pointV, 0, 0, validCells); // cell corner
}

/**
 * Computes fractional liquid SDF values for the (centerStart_x and centerStart_y) for all validCells and stores them in the plane
 */
void ViscositySolver::estimatePlaneFractions(MatrixXt &plane, doubleT centerStart_x, doubleT centerStart_y, MatrixXi &validCells) {
    MatrixXt nodalPhi(plane.rows() + 1, plane.cols() + 1);
	MatrixXi isNodalSet(plane.rows() + 1, plane.cols() + 1);
	isNodalSet.setZero();

    plane.fill(0);
    doubleT half_dx = 0.5 * m_delta_x;

    for (int j = 0; j < plane.cols(); j++) {
        for (int i = 0; i < plane.rows(); i++) {
            if (!validCells(i,j))
                continue;

            doubleT center_x = centerStart_x + i * m_delta_x;
            doubleT center_y = centerStart_y + j * m_delta_y;

            if (!isNodalSet(i,j))
            {
                doubleT n = m_fluid_domain->getLiquidSDF().bilinearInterpolate(center_x - half_dx, center_y - half_dx);
                nodalPhi(i,j) = n;
                isNodalSet(i,j) = 1;
            }
            doubleT phi00 = nodalPhi(i,j);

            if (!isNodalSet(i,j + 1))
            {
                doubleT n = m_fluid_domain->getLiquidSDF().bilinearInterpolate(center_x - half_dx, center_y + half_dx);
                nodalPhi(i,j + 1) = n;
                isNodalSet(i,j + 1) = 1;
            }
            doubleT phi01 = nodalPhi(i, j + 1);

            if (!isNodalSet(i + 1,j))
            {
                doubleT n = m_fluid_domain->getLiquidSDF().bilinearInterpolate(center_x + half_dx, center_y - half_dx);
                nodalPhi(i + 1,j) = n;
                isNodalSet(i + 1,j) = 1;
            }
            doubleT phi10 = nodalPhi(i + 1, j);

            if (!isNodalSet(i + 1,j + 1))
            {
                doubleT n = m_fluid_domain->getLiquidSDF().bilinearInterpolate(center_x + half_dx, center_y + half_dx);
                nodalPhi(i + 1,j + 1) = n;
                isNodalSet(i + 1,j + 1) = 1;
            }
            doubleT phi11 = nodalPhi(i + 1, j + 1);

            if (phi00 < 0.0 && phi01 < 0.0 && phi10 < 0.0 && phi11 < 0.0)
                plane(i,j) = 1.0;
            else if (phi00 >= 0.0 && phi01 >= 0.0 && phi10 >= 0.0 && phi11 >= 0.0)
                plane(i,j) = 0.0;
            else
                plane(i,j) = areaFraction(phi00, phi10, phi01, phi11);
        }
    }
}

/**
 * Computes the matrix index for fluid edges that will be in the solver.
 */
void ViscositySolver::computeMatrixIndexTable()
{
    int dim = (m_size_x + 1) * m_size_y + m_size_x * (m_size_y + 1);
    EdgeIndexer fidx(m_size_x, m_size_y);

    std::vector<bool> isIndexInMatrix(dim, false);
    for (int j = 1; j < m_size_y; j++) {
        for (int i = 1; i < m_size_x; i++) {
            if (m_state.U(i,j) != EdgeState::fluid)
                continue;

            doubleT v = planes.U(i,j);
            doubleT vRight = planes.center(i,j);
            doubleT vLeft = planes.center(i - 1, j);
            doubleT vTop = planes.pointV(i, j + 1);
            doubleT vBottom = planes.pointV(i, j);

            if (v > 0.0 || vRight > 0.0 || vLeft > 0.0 || vTop > 0.0 || vBottom > 0.0)
            {
                int index = fidx.U(i,j);
                isIndexInMatrix[index] = true;
            }
        }
    }
    for (int j = 1; j < m_size_y; j++) {
        for (int i = 1; i < m_size_x; i++) {
            if (m_state.V(i,j) != EdgeState::fluid)
                continue;

            doubleT v = planes.V(i,j);
            doubleT vRight = planes.pointU(i + 1, j);
            doubleT vLeft = planes.pointU(i , j);
            doubleT vTop = planes.center(i, j);
            doubleT vBottom = planes.center(i, j - 1);

            if (v > 0.0 || vRight > 0.0 || vLeft > 0.0 || vTop > 0.0 || vBottom > 0.0)
            {
                int index = fidx.V(i,j);
                isIndexInMatrix[index] = true;
            }
        }
    }

    std::vector<int> gridToMatrix(dim, -1);
    int matIndex = 0;
    for (size_t i = 0; i < isIndexInMatrix.size(); i++) {
        if (isIndexInMatrix[i]) {
            gridToMatrix[i] = matIndex++;
        }
    }
    m_matrixIndex = MatrixIndexer(m_size_x, m_size_y, gridToMatrix);
}

/**
 * Builds the linear system using the fractional planes
 */
void ViscositySolver::initializeLinearSystem(VectorXt &rhs, std::vector<Eigen::Triplet<doubleT>>& triplets) {

    doubleT invDx = m_delta_x;
    doubleT factor = m_dt * invDx * invDx;

    for (int j = 1; j < m_size_y; j++) {
        for (int i = 1; i < m_size_x; i++) {
            if (m_state.U(i,j) != EdgeState::fluid)
                continue;

            int row = m_matrixIndex.U(i, j);
            if (row == -1)
                continue;

           doubleT visc_Right = m_fluid_domain->viscosity(i, j);
           doubleT visc_Left = m_fluid_domain->viscosity(i - 1, j);

           doubleT visc_Top = 0.25 * (m_fluid_domain->viscosity(i - 1, j + 1) +
                                     m_fluid_domain->viscosity(i - 1, j) +
                                     m_fluid_domain->viscosity(i, j + 1) +
                                     m_fluid_domain->viscosity(i, j));

            doubleT visc_Bottom = 0.25 * (m_fluid_domain->viscosity(i - 1, j) +
                                      m_fluid_domain->viscosity(i - 1, j - 1) +
                                      m_fluid_domain->viscosity(i, j) +
                                      m_fluid_domain->viscosity(i, j - 1));

            doubleT v = planes.U(i,j);
            doubleT vRight = planes.center(i,j);
            doubleT vLeft = planes.center(i - 1, j);
            doubleT vTop = planes.pointV(i, j + 1);
            doubleT vBottom = planes.pointV(i, j);

            doubleT factorRight = 2 * factor * visc_Right * vRight;
            doubleT factorLeft = 2 * factor * visc_Left * vLeft;
            doubleT factorTop = factor * visc_Top * vTop;
            doubleT factorBottom = factor * visc_Bottom * vBottom;

            doubleT diag = v + factorRight + factorLeft + factorTop + factorBottom;
            triplets.emplace_back(row,row, diag);
            if (m_state.U(i + 1,j) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.U(i + 1, j), -factorRight); }
            if (m_state.U(i - 1,j) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.U(i - 1, j), -factorLeft); }
            if (m_state.U(i,j + 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.U(i , j + 1), -factorTop); }
            if (m_state.U(i,j - 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.U(i , j - 1), -factorBottom); }

            if (m_state.V(i,j + 1) == EdgeState::fluid)        { triplets.emplace_back(row, m_matrixIndex.V(i, j + 1), -factorTop); }
            if (m_state.V(i - 1,j + 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.V(i - 1, j + 1), factorTop); }
            if (m_state.V(i, j) == EdgeState::fluid)              { triplets.emplace_back(row, m_matrixIndex.V(i , j), factorBottom); }
            if (m_state.V(i - 1, j) == EdgeState::fluid)       { triplets.emplace_back(row, m_matrixIndex.V(i - 1 , j), -factorBottom); }

            doubleT rval = v * m_fluid_domain->uijHalfInd(i,j);
            if (m_state.U(i + 1,j) == EdgeState::solid) { rval -= -factorRight * m_fluid_domain->uijHalfInd(i + 1,j); }
            if (m_state.U(i - 1,j) == EdgeState::solid) { rval -= -factorLeft * m_fluid_domain->uijHalfInd(i - 1,j); }
            if (m_state.U(i,j + 1) == EdgeState::solid) { rval -= -factorTop * m_fluid_domain->uijHalfInd(i,j + 1); }
            if (m_state.U(i,j - 1) == EdgeState::solid) { rval -= -factorBottom * m_fluid_domain->uijHalfInd(i,j - 1); }

            if (m_state.V(i,j + 1) == EdgeState::solid)        { rval -= -factorTop * m_fluid_domain->vijHalfInd(i, j + 1); }
            if (m_state.V(i - 1,j + 1) == EdgeState::solid) { rval -= factorTop * m_fluid_domain->vijHalfInd(i - 1,j + 1); }
            if (m_state.V(i, j) == EdgeState::solid)              { rval -= factorBottom * m_fluid_domain->vijHalfInd(i, j); }
            if (m_state.V(i - 1,j) == EdgeState::solid)        { rval -= -factorBottom * m_fluid_domain->vijHalfInd(i - 1,j); }
            rhs[row] = rval;

        }
    }

    for (int j = 1; j < m_size_y; j++) {
        for (int i = 1; i < m_size_x; i++) {
            if (m_state.V(i,j) != EdgeState::fluid)
                continue;

            int row = m_matrixIndex.V(i, j);
            if (row == -1)
                continue;

            doubleT visc_Right = 0.25 * (m_fluid_domain->viscosity(i, j - 1) +
                                      m_fluid_domain->viscosity(i + 1, j - 1) +
                                      m_fluid_domain->viscosity(i, j) +
                                      m_fluid_domain->viscosity(i + 1, j));

            doubleT visc_Left = 0.25 * (m_fluid_domain->viscosity(i, j - 1) +
                                         m_fluid_domain->viscosity(i - 1, j - 1) +
                                         m_fluid_domain->viscosity(i, j) +
                                         m_fluid_domain->viscosity(i - 1, j));

            doubleT visc_Top = m_fluid_domain->viscosity(i, j);
            doubleT visc_Bottom = m_fluid_domain->viscosity(i, j - 1);

            doubleT v = planes.V(i,j);
            doubleT vRight = planes.pointU(i + 1,j);
            doubleT vLeft = planes.pointU(i , j);
            doubleT vTop = planes.center(i, j);
            doubleT vBottom = planes.center(i, j - 1);

            doubleT factorRight = factor * visc_Right * vRight;
            doubleT factorLeft = factor * visc_Left * vLeft;
            doubleT factorTop =  2 * factor * visc_Top * vTop;
            doubleT factorBottom = 2 * factor * visc_Bottom * vBottom;

            doubleT diag = v + factorRight + factorLeft + factorTop + factorBottom;

            triplets.emplace_back(row,row, diag);
            if (m_state.V(i + 1,j) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.V(i + 1, j), -factorRight); }
            if (m_state.V(i - 1,j) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.V(i - 1, j), -factorLeft); }
            if (m_state.V(i,j + 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.V(i , j + 1), -factorTop); }
            if (m_state.V(i,j - 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.V(i , j - 1), -factorBottom); }

            if (m_state.U(i + 1,j) == EdgeState::fluid)        { triplets.emplace_back(row, m_matrixIndex.U(i + 1, j), -factorRight); }
            if (m_state.U(i + 1,j - 1) == EdgeState::fluid) { triplets.emplace_back(row, m_matrixIndex.U(i + 1, j - 1), factorRight); }
            if (m_state.U(i, j) == EdgeState::fluid)              { triplets.emplace_back(row, m_matrixIndex.U(i , j), factorLeft); }
            if (m_state.U(i ,j - 1) == EdgeState::fluid)       { triplets.emplace_back(row, m_matrixIndex.U(i , j - 1), -factorLeft); }

            doubleT rval = v * m_fluid_domain->vijHalfInd(i,j);
            if (m_state.V(i + 1,j) == EdgeState::solid) { rval -= -factorRight * m_fluid_domain->vijHalfInd(i + 1,j); }
            if (m_state.V(i - 1,j) == EdgeState::solid) { rval -= -factorLeft * m_fluid_domain->vijHalfInd(i - 1,j); }
            if (m_state.V(i,j + 1) == EdgeState::solid) { rval -= -factorTop * m_fluid_domain->vijHalfInd(i,j + 1); }
            if (m_state.V(i,j - 1) == EdgeState::solid) { rval -= -factorBottom * m_fluid_domain->vijHalfInd(i,j - 1); }

            if (m_state.U(i + 1,j) == EdgeState::solid)         { rval -= -factorRight * m_fluid_domain->uijHalfInd(i + 1, j); }
            if (m_state.U(i + 1,j - 1) == EdgeState::solid)  { rval -= factorRight * m_fluid_domain->uijHalfInd(i + 1,j - 1); }
            if (m_state.U(i, j) == EdgeState::solid)               { rval -= factorLeft * m_fluid_domain->uijHalfInd(i, j); }
            if (m_state.U(i,j - 1) == EdgeState::solid)         { rval -= -factorLeft * m_fluid_domain->uijHalfInd(i,j - 1); }
            rhs[row] = rval;
        }
    }

    std::vector<Triplet<doubleT>> goodTriplets;
    for (auto & trip : triplets)
    {
        if ((trip.col() == -1 || trip.row() == -1) && abs(trip.value()) > 1e-9)
            std::cout << "naughtyTriplet " << trip.row() << " " << trip.col() << " " << trip.value() << "\n";
        else if (trip.col() != -1 and trip.row() != -1)
            goodTriplets.push_back(trip);
    }
    triplets = goodTriplets;
}










