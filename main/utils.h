#pragma once
#include <Eigen/Dense>
#include <igl/colormap.h>
#include "GuiData.h"

/**
	 * from exercise Grid2.h
	 * @param C Color matrix
	 * @param normalize to [0,1]?
	 * @param m_res_x grid resolution, used when C is not yet initialised
	 * @param m_res_y
	 * @param alpha set to 0 to use allow alpha, between 0,1 to use brightness, 1 to not use transparency
	 */
void getColors(MatrixXt field, Eigen::MatrixXf& C, bool normalize, int m_res_x, int m_res_y, float alpha = 0);
