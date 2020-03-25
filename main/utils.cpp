#include "utils.h"

void getColors(MatrixXt field, Eigen::MatrixXf& C, bool normalize, int m_res_x, int m_res_y, float alpha) {
	if (C.rows() == 0) {
		int num_faces = m_res_x * m_res_y * 2; // 2 triangles per cell
		C = Eigen::MatrixXf(num_faces, 4);
	}
	int i = 0;
	double cmin = field(0, 0);
	double cmax = cmin;
	for (int y = 0; y < m_res_y; ++y) {
		for (int x = 0; x < m_res_x; ++x) {
			double c = field(x, y);
			if (normalize) {
				if (c > cmax) cmax = c;
				if (c < cmin) cmin = c;
			} else {
				C.row(i++).setConstant(c);
				C.row(i++).setConstant(c);
			}

		}
	}

	if (alpha == 1)
		C.col(C.cols() - 1).setConstant(alpha);

	if (!normalize) return;
	else if (cmin == cmax) {
		C.setZero();
		C.col(C.cols() - 1).setConstant(alpha);
		return;
	}

//		std::cout << "cmin:" << cmin << " cmax:" << cmax << std::endl;
	for (int y = 0; y < m_res_y; ++y) {
		for (int x = 0; x < m_res_x; ++x) {
			double c = field(x, y);
			c = (c - cmin) / (cmax - cmin); // [0,1]
			double r, g, b;
			igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, c, r, g, b);
			C.row(i++) = Eigen::RowVector4f(r, g, b, 1.);
			C.row(i++) = Eigen::RowVector4f(r, g, b, 1.);
		}
	}
}