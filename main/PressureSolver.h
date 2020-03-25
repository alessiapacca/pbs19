//
// Created by Niklaus on 29.11.19.
//

#ifndef MAIN_PRESSURESOLVER_H
#define MAIN_PRESSURESOLVER_H

#include "FluidDomain.h"

/*
 *  MatrixCoeffs
 */

struct MatrixCell {
    double diag;
    double i_plus;
    double j_plus;

    MatrixCell() : diag(0.0), i_plus(0.0), j_plus(0.0) {}
};

class MatrixCoeffs
{
public:
    MatrixCoeffs();
    MatrixCoeffs(int size);
    ~MatrixCoeffs();

    const MatrixCell operator [](int i) const;
    MatrixCell& operator [](int i);

    inline size_t size() { return cells.size(); }

    std::vector<MatrixCell> cells;
};

/*
 * GridIndexKeyMap
 */

struct GridIndex
{
    int i; int j;
    GridIndex() = default;
    GridIndex(int _i, int _j) : i(_i), j(_j) {}
};

class GridIndexKeyMap
{
public:


    GridIndexKeyMap();
    GridIndexKeyMap(int i, int j);
    ~GridIndexKeyMap();

    void clear();
    void insert (GridIndex g, int key);
    int find(GridIndex g);
    int find(int i, int j);

private:
    inline int getFlatIndex(GridIndex g) { return g.i + i_size * g.j; }
    inline int getFlatIndex(int i, int j) { return i + i_size * j; }

    int i_size = 0;
    int j_size = 0;
    std::vector<int> indices;
    int notFoundValue = -1;

};

/*
 * PressureSolver
 */

class PressureSolver
{
public:
    PressureSolver();
    ~PressureSolver();

    Eigen::MatrixXd solve(FluidDomain* fluid_domain, double dt);

private:


    inline int gridToVectorIndex(GridIndex g) { return keyMap.find(g); }
    inline int gridToVectorIndex(int i, int j) { return keyMap.find(i, j); }

    void initialize();
    void calculateNegativeDivergenceVector(Eigen::VectorXd &b);
    void calculateMatrixCoeffs(MatrixCoeffs &A);
    void calculatePreconditionerVector(MatrixCoeffs &A, Eigen::VectorXd &precon);
    void applyMatrix(MatrixCoeffs &A, Eigen::VectorXd &x, Eigen::VectorXd &result);
    void applyPreconditioner(MatrixCoeffs &A, Eigen::VectorXd &precon, Eigen::VectorXd &residual, Eigen::VectorXd &vect);
    void solvePressureSystem(MatrixCoeffs &A, Eigen::VectorXd &b, Eigen::VectorXd &precon, Eigen::VectorXd &pressure);

    GridIndexKeyMap keyMap;

    std::vector<GridIndex> m_pressureCells;
    FluidDomain* m_fluid_domain{};
    int m_size_x = 0;
    int m_size_y = 0;
    double m_dx = 0.0;
    double m_dt = 0.0;
    int m_matrix_size = 0;

    double m_pressureSolveTolerance = 1e-9;
    int m_maxCGIterations = 200;
    double m_minFraction = 0.1;

};





#endif //MAIN_PRESSURESOLVER_H
