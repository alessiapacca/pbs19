//
// Created by Niklaus on 12.11.19.
//

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "GuiData.h"


class MarkerParticle {
public:
    MarkerParticle();
    MarkerParticle(doubleT pos_x, doubleT pos_y, doubleT vel_x = 0, doubleT vel_y = 0);
    ~MarkerParticle();

    // Getters
    inline doubleT posX() const { return m_posX; };
    inline doubleT posY() const { return m_posY; };
    inline doubleT velX() const { return m_velX; };
    inline doubleT velY() const { return m_velY; };
    inline bool isMarkedToDelete() { return markedToDelete; }

    // Setters
    inline void setPosition(doubleT pos_x, doubleT pos_y)
    {
        m_posX = pos_x;
        m_posY = pos_y;
    };
    inline void setVelocity(doubleT vel_x, doubleT vel_y)
    {
        m_velX = vel_x;
        m_velY = vel_y;
    };
    inline void advect(doubleT m_dt)
    {
        m_posX += m_velX * m_dt;
        m_posY += m_velY * m_dt;
    };

    inline void markToDelete() { markedToDelete = true; }

    doubleT T = 1;
private:
    doubleT m_posX;
    doubleT m_posY;

    doubleT m_velX;
    doubleT m_velY;

    bool markedToDelete = false;
};


#endif
