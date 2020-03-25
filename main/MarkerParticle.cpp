//
// Created by Niklaus on 12.11.19.
//
#include "MarkerParticle.h"

MarkerParticle::MarkerParticle()
        : m_posX(0), m_posY(0), m_velX(0), m_velY(0) {}

MarkerParticle::MarkerParticle(doubleT pos_x, doubleT pos_y, doubleT vel_x, doubleT vel_y)
        : m_posX(pos_x), m_posY(pos_y), m_velX(vel_x), m_velY(vel_y) {}

MarkerParticle::~MarkerParticle() = default;


