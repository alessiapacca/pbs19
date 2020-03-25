/*
    Created by Niklaus on 04.12.19.
    Method is based on https://github.com/rlguy/FLIPViscosity3D/
*/

#include "ParticleLevelSet.h"
#include <iostream>

ParticleLevelSet::ParticleLevelSet() = default;

ParticleLevelSet::ParticleLevelSet(int i, int j, double dx) : m_size_i(i), m_size_j(j), m_dx(dx) {
    m_phi = Eigen::MatrixXd::Constant(i, j, getMaxDistance());
}

ParticleLevelSet::~ParticleLevelSet() = default;

double ParticleLevelSet::operator()(int i, int j) const { return get(i, j); }

double ParticleLevelSet::get(int i, int j) const {
    assert(i >= 0 && i < m_size_i && j >= 0 && j < m_size_j);
    return m_phi(i, j);
}

double ParticleLevelSet::getEdgeWeightU(int i, int j) {
    assert(i >= 0 && i < m_size_i && j >= 0 && j < m_size_j);
    return fractionInside(m_phi(i - 1, j), m_phi(i, j));
}

double ParticleLevelSet::getEdgeWeightV(int i, int j) {
    assert(i >= 0 && i < m_size_i && j >= 0 && j < m_size_j);
    return fractionInside(m_phi(i, j - 1), m_phi(i, j));
}

void ParticleLevelSet::computeSignedDistanceFromParticles(std::vector<MarkerParticle> &particles, double radius) {
    m_phi = Eigen::MatrixXd::Constant(m_size_i, m_size_j, getMaxDistance());

    double half_dx = 0.5 * m_dx;
    for(auto & p : particles) {
        int i = std::floor(p.posX() / m_dx);
        int j = std::floor(p.posY() / m_dx);

        int gmin_i = std::max(0, i - 1);
        int gmin_j = std::max(0, j - 1);
        int gmax_i = std::min(i + 1, m_size_i - 1);
        int gmax_j = std::min(j + 1, m_size_j - 1);

        for(int j = gmin_j; j <= gmax_j; j++) {
            for(int i = gmin_i; i <= gmax_i; i++) {
                double cellX = i * m_dx + half_dx;
                double cellY = j * m_dx + half_dx;
                double dist = sqrt( pow(cellX - p.posX(), 2) + pow(cellY - p.posY(), 2)) - radius;
                if(dist < m_phi(i, j)) {
                    m_phi(i, j) = dist;
                }
            }
        }
    }
}

void ParticleLevelSet::extrapolateSignedDistanceIntoSolids() {
    for(int j = 0; j < m_size_j; j++) {
        for(int i = 0; i < m_size_i; i++) {
            bool isBorder = i == 0 || i == m_size_i - 1 || j == 0 || j == m_size_j - 1;
            if(isBorder && m_phi(i, j) < 0.5 * m_dx) {
                m_phi(i, j) = -0.5 * m_dx;
            }
        }
    }
}

void ParticleLevelSet::calculateSignedDistanceField(std::vector<MarkerParticle> &particles, double radius) {
    computeSignedDistanceFromParticles(particles, radius);
    extrapolateSignedDistanceIntoSolids();
}
