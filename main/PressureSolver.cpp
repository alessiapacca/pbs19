//
// Created by Niklaus on 29.11.19.
//
#include "PressureSolver.h"
#include <iostream>

using namespace Eigen;

/*
 * MatrixCoeffs
 */

MatrixCoeffs::MatrixCoeffs() = default;
MatrixCoeffs::MatrixCoeffs(int size) : cells(size, MatrixCell()) {}
MatrixCoeffs::~MatrixCoeffs() = default;

const MatrixCell MatrixCoeffs::operator [](int i) const {
    assert(i >= 0 && i < cells.size()); return cells[i];
}
MatrixCell& MatrixCoeffs::operator [](int i) {
    assert(i >= 0 && i < cells.size()); return cells[i];
}

/*
 *  GridIndexMap
 */

GridIndexKeyMap::GridIndexKeyMap() = default;
GridIndexKeyMap::GridIndexKeyMap(int i, int j) : i_size(i), j_size(j) {
    indices = std::vector<int>(i*j, notFoundValue);
}
GridIndexKeyMap::~GridIndexKeyMap() = default;

void GridIndexKeyMap::clear() {
    for (size_t i = 0; i < indices.size(); i++)
        indices[i] = notFoundValue;
}

void GridIndexKeyMap::insert (GridIndex g, int key) {
    int index = getFlatIndex(g);
    indices[index] = key;
}
int GridIndexKeyMap::find(GridIndex g){
    return find(g.i, g.j);
}

int GridIndexKeyMap::find(int i, int j) {
    if (indices.size() == 0)
        return notFoundValue;
    int index = getFlatIndex(i,j);
    return indices[index];
}

/*
 * PressureSolver
 */
PressureSolver::PressureSolver() = default;
PressureSolver::~PressureSolver() = default;

Eigen::MatrixXd PressureSolver::solve(FluidDomain *fluid_domain, double dt) {
    m_fluid_domain = fluid_domain;
    m_size_x = m_fluid_domain->getSizeX();
    m_size_y = m_fluid_domain->getSizeY();
    m_dx = m_fluid_domain->getDeltaX();
    m_dt = dt;

    initialize();

    Eigen::VectorXd b = Eigen::VectorXd::Zero(m_matrix_size);
    calculateNegativeDivergenceVector(b);
    if (b.cwiseAbs().maxCoeff() < m_pressureSolveTolerance)
        return Eigen::MatrixXd::Zero(m_size_x, m_size_y);

    MatrixCoeffs A(m_matrix_size);
    calculateMatrixCoeffs(A);

    Eigen::VectorXd precon = Eigen::VectorXd::Zero(m_matrix_size);
    calculatePreconditionerVector(A, precon);

    Eigen::VectorXd pressure = Eigen::VectorXd::Zero(m_matrix_size);
    solvePressureSystem(A, b, precon, pressure);

    MatrixXd pressureGrid(m_size_x, m_size_y);
    pressureGrid.setZero();
    for (size_t idx = 0; idx < m_pressureCells.size(); idx++)
    {
        GridIndex g = m_pressureCells[idx];
        pressureGrid(g.i, g.j) = pressure[idx];
    }

    return pressureGrid;
   // return Array2d();
}

void PressureSolver::initialize()
{
    m_pressureCells = std::vector<GridIndex>();
    for (int j = 0; j < m_size_y; j++) {
        for (int i = 0; i < m_size_x; i++) {
            if (m_fluid_domain->cellType(i,j) == LIQUID)
                m_pressureCells.emplace_back(i,j);
        }
    }
    m_matrix_size = m_pressureCells.size();

    keyMap = GridIndexKeyMap(m_size_x, m_size_y);
    for (size_t i = 0; i < m_pressureCells.size(); i++)
    {
        keyMap.insert(m_pressureCells[i], i);
    }
}

void PressureSolver::calculateNegativeDivergenceVector(Eigen::VectorXd &b) {

    for (auto & g : m_pressureCells)
    {
        int i = g.i; int j = g.j;

        double divergence = 0.0;
        divergence -= m_fluid_domain->uijHalfInd(i + 1, j);
        divergence += m_fluid_domain->uijHalfInd(i,j);
        divergence -= m_fluid_domain->vijHalfInd(i,j + 1);
        divergence += m_fluid_domain->vijHalfInd(i, j);
        divergence /= m_dx;
        b[gridToVectorIndex(g)] = divergence;
    }
}

void PressureSolver::calculateMatrixCoeffs(MatrixCoeffs &A) {

    double scale = m_dt / (m_dx * m_dx);
    double half_dx = 0.5 * m_dx;
    GridIndex g;
    for (size_t idx = 0; idx < m_pressureCells.size(); idx++)
    {
        g = m_pressureCells[idx];
        int i = g.i; int j = g.j;
        int index = gridToVectorIndex(g);

        //right neighbor
        double term = m_fluid_domain->cellType(i + 1, j) == SOLID ? 0.0 : scale;
        if (m_fluid_domain->cellType(i + 1, j) == LIQUID) {
            A[index].diag += term;
            A[index].i_plus -= term;
        } else {
            A[index].diag += term;
        }

        //left neighbor
        term = m_fluid_domain->cellType(i, j) == SOLID ? 0.0 : scale;
        if (m_fluid_domain->cellType(i - 1, j) == LIQUID) {
            A[index].diag += term;
        } else {
            A[index].diag += term;
        }

        //top neighbor
        term = m_fluid_domain->cellType(i, j + 1) == SOLID ? 0.0 : scale;
        if (m_fluid_domain->cellType(i, j + 1) == LIQUID) {
            A[index].diag += term;
            A[index].j_plus -= term;
        } else {
            A[index].diag += term;
        }

        //bottom neighbor
        term = m_fluid_domain->cellType(i, j) == SOLID ? 0.0 : scale;
        if (m_fluid_domain->cellType(i, j - 1) == LIQUID) {
            A[index].diag += term;
        } else {
            A[index].diag += term;
        }
    }
}

void PressureSolver::calculatePreconditionerVector(MatrixCoeffs &A, Eigen::VectorXd &precon)
{
    assert(A.size() == precon.size());

    double tau = 0.97; //Tuning constant
    double sigma = 0.25; //Safety constant

    GridIndex g;
    for (size_t idx = 0; idx < m_pressureCells.size(); idx++)
    {
        g = m_pressureCells[idx];
        int i = g.i; int j = g.j;

        int v_idx = gridToVectorIndex(i,j);
        int v_idxI = gridToVectorIndex(i - 1, j);
        int v_idxJ = gridToVectorIndex(i, j - 1);

        double diag = A[v_idx].diag;

        double i_plusI = v_idxI != -1 ? A[v_idxI].i_plus : 0.0;
        double i_plusJ = v_idxJ != -1 ? A[v_idxJ].i_plus : 0.0;

        double j_plusI = v_idxI != -1 ? A[v_idxI].j_plus : 0.0;
        double j_plusJ = v_idxJ != -1 ? A[v_idxJ].j_plus : 0.0;

        double preconI = v_idxI != -1 ? precon[v_idxI] : 0.0;
        double preconJ = v_idxJ != -1 ? precon[v_idxJ] : 0.0;

        double v1 = i_plusI * preconI;
        double v2 = j_plusJ * preconJ;
        double v3 = preconI * preconI;
        double v4 = preconJ * preconJ;

        double e = diag - v1*v1 - v2*v2 -
                tau * (i_plusI * j_plusI * v3 +
                       j_plusJ * i_plusJ * v4);
        if (e < sigma * diag)
            e = diag;

        if (abs(e) > 10e-9) {
            assert(e > 0.0);
            precon[v_idx] = 1 / sqrt(e);
        }
    }
}

void PressureSolver::applyPreconditioner(MatrixCoeffs &A, Eigen::VectorXd &precon, Eigen::VectorXd &residual, Eigen::VectorXd &vect) {

    //Solve A*q = residual
    Eigen::VectorXd q = Eigen::VectorXd::Zero(m_matrix_size);

    GridIndex g;
    for (size_t idx = 0; idx < m_pressureCells.size(); idx++) {
        g = m_pressureCells[idx];
        int i = g.i; int j = g.j;
        int v_idx = gridToVectorIndex(i,j);

        int v_idxI = gridToVectorIndex(i - 1, j);
        int v_idxJ = gridToVectorIndex(i, j - 1);

        double i_plusI = 0.0;
        double preconI = 0.0;
        double qI = 0.0;

        if (v_idxI != -1)
        {
            i_plusI = A[v_idxI].i_plus;
            preconI = precon[v_idxI];
            qI = q[v_idxI];
        }

        double j_plusJ = 0.0;
        double preconJ = 0.0;
        double qJ = 0.0;
        if (v_idxJ != -1)
        {
            j_plusJ = A[v_idxJ].j_plus;
            preconJ = precon[v_idxJ];
            qJ = q[v_idxJ];
        }

        double t = residual[v_idx] - i_plusI * preconI * qI - j_plusJ * preconJ * qJ;
        t = t * precon[v_idx];
        q[v_idx] = t;
    }

    //Solve transpose(A) * z = q;
  for (int idx = (int)m_pressureCells.size() - 1; idx >= 0; idx--)
    {
        g = m_pressureCells[idx];
        int i = g.i; int j = g.j;

        int v_idx = gridToVectorIndex(i,j);
        int v_idxI = gridToVectorIndex(i + 1, j);
        int v_idxJ = gridToVectorIndex(i, j + 1);

        double vectI = v_idxI != - 1 ? vect[v_idxI] : 0.0;
        double vectJ = v_idxJ != - 1 ? vect[v_idxJ] : 0.0;

        double i_plus = A[v_idx].i_plus;
        double j_plus = A[v_idx].j_plus;

        double preconval = precon[v_idx];
        double t = q[v_idx] - i_plus * preconval * vectI - j_plus * preconval * vectJ;
        t = t * preconval;
        vect[v_idx] = t;
    }
}

void PressureSolver::applyMatrix(MatrixCoeffs &A, Eigen::VectorXd &x, Eigen::VectorXd &result)
{
    assert(A.size() == x.size() && x.size() == result.size());

    GridIndex g;
    for (size_t idx = 0; idx < m_pressureCells.size(); idx++)
    {
        g = m_pressureCells[idx];
        int i = g.i; int j = g.j;
        int r_idx = gridToVectorIndex(i,j);

        // val = dot product of column vector x and idxth row of matrix A
        double val = 0.0;
        int v_idx = gridToVectorIndex(i - 1, j);
        if (v_idx != -1) { val += x[v_idx] * A[v_idx].i_plus; }

        v_idx = gridToVectorIndex(i + 1, j);
        if (v_idx != -1) { val += x[v_idx] * A[r_idx].i_plus; }

        v_idx = gridToVectorIndex(i, j - 1);
        if (v_idx != -1) { val += x[v_idx] * A[v_idx].j_plus; }

        v_idx = gridToVectorIndex(i, j + 1);
        if (v_idx != -1) { val += x[v_idx] * A[r_idx].j_plus; }

        val += x[r_idx] * A[r_idx].diag;
        result[r_idx] = val;
    }
}

void PressureSolver::solvePressureSystem(MatrixCoeffs &A, Eigen::VectorXd &b, Eigen::VectorXd &precon, Eigen::VectorXd &pressure)
{
    double tol = m_pressureSolveTolerance;
    if (b.cwiseAbs().maxCoeff() < tol)
        return;

    Eigen::VectorXd residual = b;
    Eigen::VectorXd auxillary = Eigen::VectorXd::Zero(m_matrix_size);
    applyPreconditioner(A, precon, residual, auxillary);

    Eigen::VectorXd search = auxillary;

    double alpha = 0.0;
    double beta = 0.0;
    double sigma = auxillary.dot(residual);
    double sigmaNew = 0.0;
    int iterationNumber = 0;

    while (iterationNumber < m_maxCGIterations)
    {
        applyMatrix(A, search, auxillary);
        alpha = sigma / auxillary.dot(search);
        pressure += search * alpha;
        residual += auxillary * (-alpha);

        if (residual.cwiseAbs().maxCoeff() < tol)
        {
            std::cout << "\n\tPressure Solver Iterations: " << iterationNumber <<
                      "\n\tEstimated Error: " << residual.cwiseAbs().maxCoeff() << "\n\n";
            return;
        }

        applyPreconditioner(A, precon, residual, auxillary);
        sigmaNew = auxillary.dot(residual);
        beta = sigmaNew / sigma;
        search = auxillary + beta * search;
        sigma = sigmaNew;

        iterationNumber++;
    }

    std::cout << "\n\tPressure Solver FAILED: " <<
                "\n\tPressure Solver Iterations: " << iterationNumber <<
                "\n\tEstimated Error: " << residual.cwiseAbs().maxCoeff() << "\n\n";

}
