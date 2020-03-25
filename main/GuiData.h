//
// Created by roger on 18.11.19.
//

#ifndef PBS_GUIDATA_H
#define PBS_GUIDATA_H
#include <vector>
#include "Eigen/Dense"

using doubleT = double; //floating point precision used in the simulation
using MatrixXt = Eigen::Matrix<doubleT, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXt = Eigen::Matrix<doubleT, Eigen::Dynamic, 1>;
using RowVector3t = Eigen::Matrix<doubleT, 1, 3>;

/**
 * Stores all the data being manipulated in the GUI
 * Access it from anywhere via GuiData::instance
 */
struct GuiData {
	static GuiData* instance;

	GuiData(){
		instance = this;
//		m_colorFluid	 << 0.9f,   1.f,  0.98f,  0.2f;	//default colors
        m_colorParticles	 << 0.,   0.6,  0.9;
        m_colorFluid	 << 0.f,   0.6f,  0.9f,  0.0f;
		m_colorAir		 << 0.9f,  1.f,   0.98f, 0.0f;
		m_colorBoundary << 0.12f, 0.15f, 0.16f, 1.f;
	}

	//--- Simulation ---
	float m_dt;
	float m_acc = 0.0001;
	int m_iter = 500;

	int m_res_x = 64;
	int m_res_y = 64;
	double m_length_x = 24;

	//--- Fluid Properties ---
	float m_viscosity = 1.0;
	float m_Tdiffusivity = 1.0;
    bool m_reseedingOn = true;

	//Decides how much of the temperature is transferred when the particles advect
	float m_particleTemperatureTransfer = 1./8.; //Assuming on avg around 8 particles per cell

	//--- Interaction ---
	const std::vector<char const*> m_brushNames = {
			"Fluid Brush", "Temperature Brush"//, "Force Brush"
	};
	static const int fluidBrush = 0, temperatureBrush = 1, forceBrush = 2;
	int m_brushMode = 0;
	float m_brushSize = 0.5;

	const std::vector<char const*> m_brushTextureNames = {
			"Rect", "Arrow", "Bunny"
	};
	static const int rectBrushTex = 0, arrowBrushTex = 1, bunnyBrushTex = 2;
	int m_brushTexture = 2;


	float m_baseTemperature = 0;
	float m_brushTemperature = 100;

//	float m_brushForce = 2.5;

	//--- Visuals ---
	const std::vector<char const*> m_fieldNames = {
			"Fluid", "Pressure", "Temperature"
	};
	static const int fluidField = 0, pressureField = 1, tempField = 2;
	int m_selectedField = 0;
	bool m_velocityOn = true;

	float m_particleSize = 3;
	Eigen::Vector3f m_colorParticles;
	Eigen::Vector4f m_colorFluid;
	Eigen::Vector4f m_colorAir;
	Eigen::Vector4f m_colorBoundary;
};


#endif //PBS_GUIDATA_H
