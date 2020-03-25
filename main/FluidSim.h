//
// Created by Niklaus on 13.11.19.
//

#ifndef MAIN_FLUIDSIM_H
#define MAIN_FLUIDSIM_H
#include <igl/edges.h>
#include "Simulation.h"
#include "FluidSolver.h"
#include "FluidDomain.h"
#include "GuiData.h"
#include "utils.h"

using namespace std;

/*
 * Simulation of a simple smoke plume rising.
 */
class FluidSim : public Simulation {
public:
    FluidSim() : Simulation() { init(); }

    virtual void init() override {
        m_res_x = GuiData::instance->m_res_x;
        m_res_y = GuiData::instance->m_res_y;
        fluid_domain = new FluidDomain(m_res_x, m_res_y, 1, GuiData::instance->m_length_x);
        FluidSolverMemoryPool mem_pool(*fluid_domain);
        fluid_solver = new FluidSolver(mem_pool);
        //add particles somehow somewhen somewhere <3
        fluid_domain->getMesh(m_renderV, m_renderF);
        m_dt = 0.005 * sqrt((m_res_x + m_res_y) * 0.5);

        reset();
    }

    virtual void resetMembers() override {
        fluid_domain->reset();
//        fluid_domain->addSource(0.5,0.7,0.2,0.90, 0.2);
        fluid_domain->setViscosity(GuiData::instance->m_viscosity);
        fluid_domain->classifyCells();
        std::cout << "resetMembers()" << std::endl;
    }

    virtual void updateRenderGeometry() override {
	    if(GuiData::instance->m_selectedField == GuiData::fluidField)
		    fluid_domain->getColors(m_renderC);
	    else if(GuiData::instance->m_selectedField == GuiData::pressureField)
		    getColors(fluid_domain->getPressure(), m_renderC, true, m_res_x, m_res_y, 1);
	    else if(GuiData::instance->m_selectedField == GuiData::tempField)
		    getColors(fluid_domain->getTemperature(), m_renderC, true, m_res_x, m_res_y, 1);
        if(GuiData::instance->m_velocityOn)
        	fluid_domain->updateEdges(0.4);
    }

    virtual bool advance() override {
        fluid_solver->stepPICFLIP(*fluid_domain, m_dt);
        // advance m_time
        m_time += m_dt;
        m_step++;

        return false;
    }

    virtual void renderRenderGeometry(
            igl::opengl::glfw::Viewer& viewer) override {
        viewer.data().set_mesh(m_renderV, m_renderF);
        viewer.data().set_colors(m_renderC.cast<double>());
        Eigen::MatrixXd particlePoints;
        Eigen::MatrixXd particleColors;
        particlePoints.resize(fluid_domain->getParticles().size(),3);
        particleColors.resize(fluid_domain->getParticles().size(),3);
        int i = 0;
        for (auto it = fluid_domain->getParticles().begin(); it != fluid_domain->getParticles().end() && i < particlePoints.rows(); it++)
        {
            particlePoints.row(i) << it->posX(), it->posY(), 0.0;
            particleColors.row(i) = GuiData::instance->m_colorParticles.cast<double>(); //TODO: add particle temperature/pressure
            i++;
        }
        viewer.data().add_points(particlePoints, particleColors);
	    viewer.data().point_size = GuiData::instance->m_particleSize;

        //viewer.data().add_edges(fluid_domain->s(), fluid_domain->e(), Eigen::RowVector3d(40, 40, 40));
	    if(GuiData::instance->m_velocityOn)
		    viewer.data().add_edges(fluid_domain->vs(), fluid_domain->ve(), fluid_domain->vc());
    }

	/**
	 * Called when you left click with the *nearest vertex* position
	 * @param pos Position in world coordinates
	 */
	virtual void OnMouseClicked(Eigen::Vector3d pos, bool invertTool) override {
		std::cout << "OnMouseClicked: " << pos << std::endl;

		if ((GuiData::instance->m_brushMode == GuiData::fluidBrush && !invertTool) || (GuiData::instance->m_brushMode == GuiData::temperatureBrush && invertTool)) {
			pos(0) /= fluid_domain->getLengthX();
			pos(1) /= fluid_domain->getLengthY();

			fluid_domain->addSource(
					pos(0) - GuiData::instance->m_brushSize / 2,
					pos(0) + GuiData::instance->m_brushSize / 2,
					pos(1) - GuiData::instance->m_brushSize / 2,
					pos(1) + GuiData::instance->m_brushSize / 2,
					0.2);
		}
		else if ((GuiData::instance->m_brushMode == GuiData::temperatureBrush && !invertTool) || (GuiData::instance->m_brushMode == GuiData::fluidBrush && invertTool)) {
			pos(0) *= (double)fluid_domain->getSizeX() / fluid_domain->getLengthX();
			pos(1) *= (double)fluid_domain->getSizeY() / fluid_domain->getLengthY();

			fluid_domain->setTemperature(GuiData::instance->m_brushTemperature,
                         pos(0) - GuiData::instance->m_brushSize * m_res_x / 2,
                         pos(0) + GuiData::instance->m_brushSize * m_res_x / 2,
                         pos(1) - GuiData::instance->m_brushSize * m_res_y / 2,
                         pos(1) + GuiData::instance->m_brushSize * m_res_y / 2);
		}
		else if (GuiData::instance->m_brushMode == GuiData::forceBrush){
//			pos(0) *= (double)fluid_domain->getSizeX() / fluid_domain->getLengthX();
//			pos(1) *= (double)fluid_domain->getSizeY() / fluid_domain->getLengthY();
			fluid_domain->AddPointForce(pos(0), pos(1));
		}
		
		updateRenderGeometry();
	}

#pragma region SettersAndGetters
    void setResX(int r) { m_res_x = r; }
    void setResY(int r) { m_res_y = r; }

    int getResX() const { return m_res_x; }
    FluidDomain* getFluidDomain() {return fluid_domain;}
    int getResY() const { return m_res_y; }
    double getTimestep() const { return m_dt; }
	virtual int getParticleCount() const override {
    	return fluid_domain->getParticles().size();
    }
    virtual int getFluidCellCount() const override {
	    return fluid_domain->getNumOfFluidCells();
	}

#pragma endregion SettersAndGetters

    int m_res_x, m_res_y;
    FluidDomain* fluid_domain;
    FluidSolver* fluid_solver;

    Eigen::MatrixXd m_renderV; // vertex positions,
    Eigen::MatrixXi m_renderF; // face indices
    Eigen::MatrixXf m_renderC; // face colors for rendering
};
#endif //MAIN_FLUIDSIM_H
