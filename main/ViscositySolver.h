/*
    Created by Niklaus on 28.11.19.
    Method is based on https://github.com/rlguy/FLIPViscosity3D/
*/

#ifndef MAIN_VISCOSITYSOLVER_H
#define MAIN_VISCOSITYSOLVER_H

#include "FluidDomain.h"
#include "GuiData.h"

#include <utility>

class ViscositySolver {

public:
    ViscositySolver() {};

    void applyViscosityToVelocityField(FluidDomain* fluid_domain, doubleT dt);

private:

    enum class EdgeState : char {
        air = 0x00,
        fluid = 0x01,
        solid = 0x02
    };

    struct EdgeStateGrid
    {
	    Eigen::Matrix<EdgeState, Eigen::Dynamic, Eigen::Dynamic> U;
	    Eigen::Matrix<EdgeState, Eigen::Dynamic, Eigen::Dynamic> V;
        EdgeStateGrid() {};
        EdgeStateGrid(int i_size, int j_size) :
            U(i_size + 1, j_size),
            V(i_size, j_size + 1) {
        	U.setConstant(EdgeState::air);
        	V.setConstant(EdgeState::air);
        }
    };

    struct EdgeIndexer
    {
        int i_size, j_size;
        EdgeIndexer() {};
        EdgeIndexer(int i, int j) : i_size(i), j_size(j) { v_offset = (i_size + 1) * j_size; }
        int U(int i, int j) { return i + (i_size + 1) * j; }
        int V(int i, int j) { return v_offset + i + i_size * j; }

        private:
            int v_offset;
    };

    struct MatrixIndexer {
        std::vector<int> indexTable;
        EdgeIndexer edgeIndexer;
        int matrixSize;

        MatrixIndexer() {};
        MatrixIndexer(int i, int j, std::vector<int> matrixIndexTable) :
            indexTable(std::move(matrixIndexTable)), edgeIndexer(i,j)
        {
            int matsize = 0;
            for (int k : indexTable)
            {
                if (k != -1)
                    matsize++;
            }
            matrixSize = matsize;
        }

        int U (int i, int j) { return indexTable[edgeIndexer.U(i,j)]; }
        int V (int i, int j) { return indexTable[edgeIndexer.V(i,j)]; }
    };

    struct ViscosityPlaneGrid {
        MatrixXt center;
	    MatrixXt U;
	    MatrixXt V;
	    MatrixXt pointU;
	    MatrixXt pointV;

        ViscosityPlaneGrid() {};
        ViscosityPlaneGrid (int i, int j) {
	        center = MatrixXt::Zero(i, j);
	        U = MatrixXt::Zero(i + 1, j);
	        V = MatrixXt::Zero(i, j + 1);
	        pointU = MatrixXt::Zero(i, j + 1);
	        pointV = MatrixXt::Zero(i + 1, j);
        }

        void destroy() {
            center = MatrixXt::Zero(0,0);
            U = MatrixXt::Zero(0,0);
            V = MatrixXt::Zero(0,0);
            pointU = MatrixXt::Zero(0,0);
            pointV = MatrixXt::Zero(0,0);
        }
    };

    // taken from https://github.com/rlguy/FLIPViscosity3D/blob/master/src/levelsetutils.cpp
    doubleT areaFraction(doubleT phi00, doubleT phi10, doubleT phi01, doubleT phi11)
    {
        doubleT phimid = 0.25 * (phi00 + phi01 + phi10 + phi11);
        return 0.25 * (areaFraction(phi00, phi10, phimid) +
                       areaFraction(phi10, phi11, phimid) +
                       areaFraction(phi11, phi01, phimid) +
                       areaFraction(phi01, phi00, phimid));
    }

    // taken from https://github.com/rlguy/FLIPViscosity3D/blob/master/src/levelsetutils.cpp
    doubleT areaFraction(doubleT phi0, doubleT phi1, doubleT phi2)
    {
        if (phi0 < 0.0) {
            if (phi1 < 0.0) {
                if (phi2 < 0.0)
                    return 0.0;
                else
                    return 1.0 - sortedTriangleFraction(phi2, phi0, phi1);
            } else if (phi2 < 0.0) {
                return 1.0 - sortedTriangleFraction(phi1, phi2, phi0);
            } else {
                return sortedTriangleFraction(phi0, phi1, phi2);
            }
        } else if (phi1 < 0.0) {
            if (phi2 < 0.0)
                return 1.0 - sortedTriangleFraction(phi0, phi1, phi2);
            else
                return sortedTriangleFraction(phi1, phi2, phi0);
        } else if (phi2 < 0.0)
            return  sortedTriangleFraction(phi2, phi0, phi1);
        else
            return 0.0;
    }

    // taken from https://github.com/rlguy/FLIPViscosity3D/blob/master/src/levelsetutils.cpp
    doubleT sortedTriangleFraction(doubleT phi0, doubleT phi1, doubleT phi2)
    {
        return phi0 * phi0 / (2.0 * (phi0 - phi1) * (phi0 - phi2));
    }

    int m_size_x, m_size_y;
    doubleT m_dt;
    doubleT m_delta_x, m_delta_y;
    FluidDomain* m_fluid_domain;
    EdgeStateGrid m_state;
    ViscosityPlaneGrid planes;
    MatrixIndexer m_matrixIndex;


    void computeEdgeStateGrid();
    void computePlaneGrid();
    void estimatePlaneFractions(MatrixXt & plane, doubleT centerStart_x, doubleT centerStart_y, Eigen::MatrixXi &validCells);
    void computeMatrixIndexTable();
    void initializeLinearSystem(VectorXt& rhs, std::vector<Eigen::Triplet<doubleT>>& triplets);

};



#endif //MAIN_VISCOSITYSOLVER_H
