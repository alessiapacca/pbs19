/*
    Created by Niklaus on 04.12.19.
    Method is based on https://github.com/rlguy/FLIPViscosity3D/
*/

#ifndef MAIN_PARTICLELEVELSET_H
#define MAIN_PARTICLELEVELSET_H

#include "MarkerParticle.h"
#include <Eigen/SparseCore>

class ParticleLevelSet {

public:
    ParticleLevelSet();
    ParticleLevelSet(int i, int j, double dx);
    ~ParticleLevelSet();

    double operator()(int i, int j) const;
    double get(int i, int j) const;
    void set(int i, int j, double value) { m_phi(i,j) = value; }
    double getEdgeWeightU(int i, int j);
    double getEdgeWeightV(int i, int j);

    void calculateSignedDistanceField(std::vector<MarkerParticle> &particles, double radius);
    inline double bilinearInterpolate(double x, double y) const
    {
        // Calculate indices
        int i = std::floor(x / m_dx);
        int j = std::floor(y / m_dx);

        i = std::max(i, 0); i = std::min(i, m_size_i - 1);
        j = std::max(j, 0); j = std::min(j, m_size_j - 1);
        double i_frac = x / m_dx - i;
        double j_frac = y / m_dx - j;

        int i_plus1 = std::max(i + 1, 0); i_plus1 = std::min(i_plus1, m_size_i - 1);
        int j_plus1 = std::max(j + 1, 0); j_plus1 = std::min(j_plus1, m_size_j - 1);

        // First interpolate in x, then in y
        double value_00 = m_phi(i, j);
        double value_10 = m_phi(i_plus1, j);
        double value_01 = m_phi(i, j_plus1);
        double value_11 = m_phi(i_plus1, j_plus1);

        // Interpolate x
        double value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
        double value_1 = (1 - i_frac) * value_01 + i_frac * value_11;

        // Interpolate y
        double value = (1 - j_frac) * value_0 + j_frac * value_1;
        return value;
    }

private:

    inline double getMaxDistance() { return 3.0 * m_dx; }

    double fractionInside(double phiLeft, double phiRight) {
        if(phiLeft < 0.0 && phiRight < 0.0) {
            return 1.0;
        }
        if (phiLeft < 0.0 && phiRight >= 0.0) {
            return phiLeft / (phiLeft - phiRight);
        }
        if(phiLeft >= 0.0 && phiRight < 0.0) {
            return phiRight / (phiRight - phiLeft);
        }

        return 0.0;
    }

    void computeSignedDistanceFromParticles(std::vector<MarkerParticle> &particles, double radius);
    void extrapolateSignedDistanceIntoSolids();

    int m_size_i = 0;
    int m_size_j = 0;
    double m_dx = 0.0;
    Eigen::MatrixXd m_phi;

};

#endif //MAIN_PARTICLELEVELSET_H
