//
// Created by Niklaus on 06.11.19.
//

#include <igl/writeOFF.h>
#include <thread>
#include "Gui.h"
#include "Simulator.h"
#include "FluidSim.h"
#include "GuiData.h"

/*
 * GUI for the spring simulation. This time we need additional paramters,
 * e.g. which integrator to use for the simulation and the force applied to the
 * canonball, and we also add some more visualizations (trajectories).
 */
class MeltdownGui : public Gui {
public:
	FluidSim* p_fluidsim = nullptr;
	GuiData* data = nullptr;
	FluidDomain* fluid_domain = nullptr;

    MeltdownGui() {
	    data = new GuiData(); //create an instance of the data and initialise the values
	    p_fluidsim = new FluidSim();
	    fluid_domain = p_fluidsim->getFluidDomain();
	    data->m_dt = p_fluidsim->getTimestep();
        setSimulation(p_fluidsim);
        start();
    }

    virtual void updateSimulationParameters() override {
        // change all parameters of the simulation to the values that are set in the GUI
        p_fluidsim->setTimestep(data->m_dt);
    }

    virtual void clearSimulation() override {}

    virtual void drawSimulationParameterMenu() override {
        ImGui::InputFloat("dt", &data->m_dt, 0.001, 0.005, 5);
        ImGui::InputFloat("Pressure accuracy", &data->m_acc, 1e-6, 1e-5, 8);
		ImGui::InputInt("Pressure iter", &data->m_iter, 100, 1000);

	    if (ImGui::CollapsingHeader("Fluid Properties", ImGuiTreeNodeFlags_DefaultOpen)) {
		    ImGui::InputFloat("Viscosity", &data->m_viscosity, 0.1, 1.0);
		    ImGui::InputFloat("Temp. Diffusivity", &data->m_Tdiffusivity, 0.1, 1.0);
		    ImGui::SliderFloat("Particle Temp. Transfer", &data->m_particleTemperatureTransfer, 0.f, 1.f);
            ImGui::Checkbox("Reseed Particles", &data->m_reseedingOn);
	    }

	    if (ImGui::CollapsingHeader("Visual", ImGuiTreeNodeFlags_DefaultOpen)) {
		    if (ImGui::Combo("Show Field", &data->m_selectedField, data->m_fieldNames.data(),
		                     data->m_fieldNames.size())) {
			    p_fluidsim->updateRenderGeometry();
		    }
		    ImGui::Checkbox("Show Velocity", &data->m_velocityOn);
		    ImGui::InputFloat("Particle Size", &data->m_particleSize, 2.f, 5.f);
		    ImGui::ColorEdit3("Particles", data->m_colorParticles.data());
		    ImGui::ColorEdit4("Fluid", data->m_colorFluid.data());
		    ImGui::ColorEdit4("Air", data->m_colorAir.data());
		    ImGui::ColorEdit4("Boundary", data->m_colorBoundary.data());
	    }

	    if (ImGui::CollapsingHeader("Interaction", ImGuiTreeNodeFlags_DefaultOpen)) {
		    ImGui::Combo("Brush Mode", &data->m_brushMode, data->m_brushNames.data(),
		                 data->m_brushNames.size());
		    ImGui::Combo("Brush Texture", &data->m_brushTexture, data->m_brushTextureNames.data(),
		                 data->m_brushTextureNames.size());
		    ImGui::SliderFloat("Brush Size", &data->m_brushSize, 0.f, 1.f);
		    ImGui::InputFloat("Base Temperature", &data->m_baseTemperature, 1, 2.5);
		    ImGui::InputFloat("Brush Temperature", &data->m_brushTemperature, 1, 2.5);
//		    ImGui::InputFloat("Brush Force", &data->m_brushForce, 1, 2.5);
	    }
    }
};

int main(int argc, char* argv[]) {
    // create a new instance of the GUI for the meltdown simulation
	Eigen::initParallel();
	std::cout << "Eigen Threads: " << Eigen::nbThreads() << std::endl;
    new MeltdownGui();
    return 0;
};