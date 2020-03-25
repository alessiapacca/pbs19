//
// Created by Niklaus on 12.11.19.
//https://github.com/kbladin/Fluid_Simulation This project helped us a lot to get started, and some function are based from there.
//

#include "FluidDomain.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <fstream>

using namespace Eigen;

FluidDomain::FluidDomain(int m_size_x, int m_size_y, int m_size_z, doubleT m_length_x)
    : m_size_x(m_size_x), m_size_y(m_size_y), m_size_z(m_size_z), m_length_x(m_length_x), rng(124)
{
	m_length_y = GuiData::instance->m_length_x * (doubleT)m_size_y / (doubleT)m_size_x;
	m_length_z = GuiData::instance->m_length_x * (doubleT)m_size_z / (doubleT)m_size_x;

    m_cell_type = Matrix<CellType, Dynamic, Dynamic>(m_size_x, m_size_y);
    m_u = MatrixXt(m_size_x + 1, m_size_y); //horizontal velocity field
    m_v = MatrixXt(m_size_x, m_size_y + 1); //vertical velocity filed
    m_u_tmp = MatrixXt(m_size_x + 1, m_size_y); //horizontal velocity field
    m_v_tmp = MatrixXt(m_size_x, m_size_y + 1); //vertical velocity filed

    m_viscosity = MatrixXt(m_size_x, m_size_y);

    m_u_prev = MatrixXt(m_size_x + 1, m_size_y); //for FLIP
    m_v_prev = MatrixXt(m_size_x, m_size_y + 1); //for FLIP
    m_u_diff = MatrixXt(m_size_x + 1, m_size_y); //for FLIP
    m_v_diff = MatrixXt(m_size_x, m_size_y + 1); //for FLIP

    m_delta_x = m_length_x / m_size_x; // cell width x
    m_delta_y = m_length_y / m_size_y; // cell width y
    m_delta_z = m_length_z / m_size_z; // cell width z

    m_p = MatrixXt(m_size_x, m_size_y);
    m_T = MatrixXt(m_size_x, m_size_y);
    buildMesh();
    buildGrid();

    particleRadius = m_delta_x * 1.01 * sqrt(3.0) / 2.0;
    liquidSDF = ParticleLevelSet(m_size_x, m_size_y, m_delta_x);

    buildTemperatureMatrix();

    //import brushed from ascii art
	readBrushTexture(arrowBrushImage, "../../data/arrow.txt", 6, 6);
	readBrushTexture(bunnyBrushImage, "../../data/bunny.txt", 345, 337);
}

FluidDomain::~FluidDomain() = default;

void FluidDomain::classifyCells()
{
    liquidSDF.calculateSignedDistanceField(m_particles, particleRadius);
    m_cell_type.fill(AIR);
    int m_num_liquid = 0;

    for (int i = 0; i < m_size_x; ++i) {
        for (int j = 0; j < m_size_y; ++j) {
            if (liquidSDF(i,j) < 0.0)
                m_cell_type(i,j) = LIQUID;
        }
    }
    // Loop through all particles and set cells to liquid if they
    // contain a particle
    /*for (auto & m_particle : m_particles)
    {
        // Find the particles position in the grid
        int x = doubleT_to_floor(m_particle.posX() / m_delta_x);
        int y = doubleT_to_floor(m_particle.posY() / m_delta_y);

        x = CLAMP(x, 0, m_size_x-1);
        y = CLAMP(y, 0, m_size_y-1);

        m_cell_type(x,y) = LIQUID;
    }*/

    // Reset border
    for (int i = 0; i < m_size_x; ++i)
    {
        for (int j = 0; j < m_size_y; ++j)
        {
            if (i == 0 || j == 0 || i == m_size_x - 1 || j == m_size_y - 1) {
                m_cell_type(i,j) = SOLID;
            }
            else if (m_cell_type(i,j) == LIQUID)
                m_num_liquid++;
        }
    }
    m_num_of_fluid_cells = m_num_liquid;
    std::cout << "#LIQUID cells=" << m_num_liquid << std::endl;
}


void FluidDomain::buildMesh() {
    int num_vertices = (m_size_x + 1) * (m_size_y + 1) /** (m_size_z + 1)*/;
    int num_faces = m_size_x * m_size_y * m_size_z * 2; // 6 triangles per cell in 3D

    m_V = Eigen::MatrixXd(num_vertices, 3);
    m_F = Eigen::MatrixXi(num_faces, 3);

    int i = 0;
	for (int z = 0; z <= 0/*m_size_z*/; ++z) {
		for (int y = 0; y <= m_size_y; ++y) {
			for (int x = 0; x <= m_size_x; ++x) {
				m_V.row(i++) = Eigen::RowVector3d(x * m_delta_x, y * m_delta_y, -z * m_delta_z);
			}
		}
	}

	const auto getVertId = [&](int x, int y, int z){
    	return x + y * (m_size_x + 1) + z * (m_size_x + 1) * (m_size_y + 1);
    };

    i = 0;
    for (int z = 0; z < m_size_z; ++z) {
        for (int y = 0; y < m_size_y; ++y) {
            for (int x = 0; x < m_size_x; ++x) {
				const int vert = getVertId(x,y,z);
				//front quad - normal in z axis
	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x+1,y,z), getVertId(x+1,y+1,z));
	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x+1,y+1,z), getVertId(x,y+1,z));
	            //lower quad - y
//	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x+1,y,z), getVertId(x+1,y,z+1));
//	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x+1,y,z+1), getVertId(x,y,z+1));
//	            //left quad - x
//	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x,y+1,z), getVertId(x,y+1,z+1));
//	            m_F.row(i++) = Eigen::RowVector3i(vert, getVertId(x,y+1,z+1), getVertId(x,y,z+1));
            }
        }
    }
}

void FluidDomain::buildGrid() {
    int num_edges = 2 * m_size_x * m_size_y + m_size_x + m_size_y;
    m_start = Eigen::MatrixXd(num_edges, 3);
    m_end = Eigen::MatrixXd(num_edges, 3);

    int i = 0;
    for (int y = 0; y <= m_size_y; ++y)
        for (int x = 0; x <= m_size_x; ++x) {
            if (x < m_size_x) {
                m_start.row(i) = Eigen::RowVector3d(x * m_delta_x, y * m_delta_y, 0);
                m_end.row(i++) = Eigen::RowVector3d((x + 1) * m_delta_x, y * m_delta_y, 0);
            }
            if (y < m_size_y) {
                m_start.row(i) = Eigen::RowVector3d(x * m_delta_x, y * m_delta_y, 0);
                m_end.row(i++) = Eigen::RowVector3d(x * m_delta_x, (y + 1) * m_delta_y, 0);
            }
        }
}

void FluidDomain::updateEdges(doubleT scale = 1) {
    int num_edges = m_size_x * m_size_y;
    if (m_vs.rows() == 0) {
        m_vs = Eigen::MatrixXd(num_edges, 3);
        m_ve = Eigen::MatrixXd(num_edges, 3);
        m_vc = Eigen::MatrixXd(num_edges, 3);
    }

    int k = 0;
    for (int j = 0; j <= m_size_y; ++j)
        for (int i = 0; i <= m_size_x; ++i) {
            doubleT x = (i + 0.5) * m_delta_x;
            doubleT y = (j + 0.5) * m_delta_y;
            /*if (j < m_size_y) {
                m_vc.row(k) = Eigen::RowVector3d(1, 0, 0);
                m_vs.row(k) = Eigen::RowVector3d(i, j + 0.5, 0);
                m_ve.row(k++) = Eigen::RowVector3d(i + m_u(i, j) * scale, j + 0.5, 0);
            }

            if (i < m_size_x) {
                m_vc.row(k) = Eigen::RowVector3d(0, 0, 1);
                m_vs.row(k) = Eigen::RowVector3d(i + 0.5, j, 0);
                m_ve.row(k++) = Eigen::RowVector3d(i + 0.5, j + m_v(i, j) * scale, 0);
            }*/
            if (j < m_size_y && i < m_size_x) {
                m_vc.row(k) = Eigen::RowVector3d(0.12, 0.15, 0.16);
                m_vs.row(k) = Eigen::RowVector3d(x, y, 0);
                doubleT x_e = uij(i,j);
                doubleT y_e = vij(i,j);
                doubleT length = sqrt(pow(x_e,2) + pow(y_e,2));
                if (length > m_delta_x)
                {
                    x_e = (x_e / length) * (m_delta_x * 0.75);
                    y_e = (y_e / length) * (m_delta_y * 0.75);
                }
                m_ve.row(k++) = Eigen::RowVector3d(x + x_e , y + y_e, 0);
            }
        }
}

void FluidDomain::getMesh(MatrixXd& V, Eigen::MatrixXi& F) const {
    V = m_V;
    F = m_F;
}

void FluidDomain::getColors(Eigen::MatrixXf& C, bool normalize) const {
	std::cout << C.rows() << ", " << C.cols() << std::endl;
    if (C.rows() == 0) {
        int num_faces = m_size_x * m_size_y * 2; // 2 triangles per cell
        C = Eigen::MatrixXf(num_faces, 4);
    }
    int i = 0;
    for (int y = 0; y < m_size_y; ++y) {
        for (int x = 0; x < m_size_x; ++x) {
			Eigen::Vector4f col;
            switch (cellType(x,y))
            {
                case LIQUID:
					col = GuiData::instance->m_colorFluid;
                    break;
                case AIR:
					col = GuiData::instance->m_colorAir;
                    break;
                case SOLID:
					col = GuiData::instance->m_colorBoundary;
                    break;
            }
	        C.row(i++) = col;
	        C.row(i++) = col;
//	        C.row(i++) = col;
//	        C.row(i++) = col;
//	        C.row(i++) = col;
//	        C.row(i++) = col;
        }
    }
}

//randomness: 0=grid, 1=random (note: 0 causes assertion)
void FluidDomain::addSource(doubleT xmin, doubleT xmax, doubleT ymin, doubleT ymax, doubleT randomness)
{
	MatrixXb* image = nullptr;

	if (GuiData::instance->m_brushTexture != GuiData::rectBrushTex){ //prevent stretching image brushes
		doubleT ycenter = 0.5 * (ymin + ymax);
		doubleT ywidthNewHalf = (xmax - xmin) * m_size_x / m_size_y / 2;
		ymin = ycenter - ywidthNewHalf;
		ymax = ycenter + ywidthNewHalf;
		if(GuiData::instance->m_brushTexture == GuiData::arrowBrushTex)
			image = &arrowBrushImage;
		else if (GuiData::instance->m_brushTexture == GuiData::bunnyBrushTex)
			image = &bunnyBrushImage;
	}

	const doubleT x_incr = m_delta_x / 2.828;
	const doubleT y_incr = m_delta_y / 2.828;

	const doubleT xsize = xmax - xmin;
	const doubleT ysize = ymax - ymin;

	const auto& clamp01 = [](doubleT value){
		return std::max((doubleT)0, std::min(value, (doubleT)1));
	};

	//returns whether to paint a particular pixel or not using a boolean matrix
	const auto& getPixelValue = [&](doubleT x, doubleT y){
		//convert from world space to pixel space
		//note: first row and last col seem to be skipped somehow
		int xpixel = clamp01(((x / m_length_x) - xmin) / xsize) * (image->rows() -1);
		int ypixel = clamp01(((y / m_length_y) - ymin) / ysize) * (image->cols() -1);
		return (*image)(xpixel, image->cols() - ypixel -1);
	};

	//Add some randomness to the initial positions
	std::normal_distribution<doubleT> gaussian_x(0, x_incr * randomness);
	std::normal_distribution<doubleT> gaussian_y(0, y_incr * randomness);

	for (doubleT y = (ymin * m_length_y); y < (ymax * m_length_y); y += y_incr) {
		for (doubleT x = (xmin * m_length_x); x < (xmax * m_length_x); x += x_incr) {
			const doubleT posx = x + gaussian_x(rng);
			const doubleT posy = y + gaussian_y(rng);
			if (posx > m_delta_x && posy > m_delta_y && posx < m_length_x - m_delta_x &&
			    posy < m_length_y - m_delta_y) //only add particle if its inside
			{
				if(GuiData::instance->m_brushTexture == GuiData::rectBrushTex) {
					addParticle(MarkerParticle(posx, posy, 0, 0));
				}
				else
				{
					if (getPixelValue(posx, posy)) //apply mask
						addParticle(MarkerParticle(posx, posy, 0, 0));
				}
			}
		}
	}

}

void FluidDomain::setTemperature(doubleT T, doubleT xmin, doubleT xmax, doubleT ymin, doubleT ymax){
	for (int y = ymin; y < ymax; y++) {
		for (int x = xmin; x < xmax; x++) {
			if(x >= 0 && y >= 0 && x < m_T.rows() && y < m_T.cols())
				set_Tij(x, y, T);
		}
	}
}

void FluidDomain::advect(doubleT dt)
{
	#pragma omp parallel for
    for (unsigned int i = 0; i < m_particles.size(); ++i)
    {
        m_particles[i].advect(dt);
    }
}

void FluidDomain::advectAndEnsureOutsideObstacles(doubleT dt)
{
    //Valid fluid domain
    doubleT v = -2.0f * m_delta_x - 1e-4;
    doubleT minx = - 0.5 * v;
    doubleT miny = - 0.5 * v;
    doubleT maxx = minx + m_length_x + v;
    doubleT maxy = miny + m_length_y + v;

    for (auto & m_particle : m_particles)
    {
        m_particle.advect(dt);

        int x = doubleT_to_floor(m_particle.posX() / m_delta_x);
        int y = doubleT_to_floor(m_particle.posY() / m_delta_y);
        x = CLAMP(x, 0, m_size_x - 1);
        y = CLAMP(y, 0, m_size_y - 1);

        if (isBorder(x, y))
        {
            //Get nearest point in valid fluid domain
            doubleT posX = std::max(m_particle.posX(), minx);
            doubleT posY = std::max(m_particle.posY(), miny);
            posX = std::min(posX, maxx - (doubleT)1e-6);
            posY = std::min(posY, maxy - (doubleT)1e-6);
            m_particle.setPosition(posX,posY);
        }
    }
}

/**
 * Precomputes the matrix to solve for the temperature implicitly
 * note: to apply different timestep/grid size/diffusivity you need to reset to call this function again
 */
void FluidDomain::buildTemperatureMatrix(){
	using D = doubleT; //floating point data type
	D k = GuiData::instance->m_Tdiffusivity;
	D factorX = k * GuiData::instance->m_dt / getDeltaX() / getDeltaX(); // scalar coefficient
	D factorY = k * GuiData::instance->m_dt / getDeltaY() / getDeltaY();

//	std::vector<D>& Tarr = m_T.data();
	auto T = VectorXt::Map(m_T.data(), m_T.size()); //raw reference to std::vector<d> as an Eigen VectorXd

	//create Eigen sparse matrix from triplets
	m_TImplMatrix = Eigen::SparseMatrix<D>(T.size(), T.size());
	typedef Eigen::Triplet<D> Triplet; //(i, j, value) one element in a SparseMatrix
	std::vector<Triplet> triplets;
	triplets.reserve(3 * m_TImplMatrix.rows());

	//Create matrix A = I - D from paper 1 or (17) from 2dheat.pdf, create one large matrix
	int nx = m_T.cols();
    int ny = m_T.rows();
	for(int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			int t = i*ny + j; //index in T
			if(i == 0 || i == nx -1 || j == 0 || j == ny -1) {
				triplets.push_back({t, t, 1});
			}
			else {
				triplets.push_back({t, t, 1 + 2 * factorX + 2 * factorY});
				triplets.push_back({t, t - 1, -factorX});
				triplets.push_back({t, t + 1, -factorX});
				triplets.push_back({t, t - ny, -factorY});
				triplets.push_back({t, t + ny, -factorY});
			}
		}
	}

	m_TImplMatrix.setFromTriplets(triplets.begin(), triplets.end());

	//precompute the solver
	m_TSolver.compute(m_TImplMatrix);
}

/**
 * Converts an ascii image to a boolean matrix to be used as a mask
 * @param drawChar write true when this character appears in a pixel
 */
void FluidDomain::readBrushTexture(MatrixXb& out, std::string imageFileName, int x_res, int y_res, char drawChar){
	using namespace std;
	out = MatrixXb::Zero(x_res, y_res);

	ifstream fin;
	fin.open(imageFileName, ios::in);
	if (fin.bad() || fin.fail())
		cerr << "Can't open image file: " + imageFileName << endl;

	char c;
	for (int y = 0; y < out.cols(); ++y) {
		for (int x = 0; x < out.rows(); ++x) {
			fin.get(c);
			if(c == '\n')
				fin.get(c);

//			std::cout << c;
			if(c == drawChar)
				out(x, y) = true;
		}
	}
	fin.close();
}

//Still in progress, has some weird artefacts
void FluidDomain::AddPointForce(doubleT x, doubleT y){
//	doubleT strength = GuiData::instance->m_brushForce;
//	for (int j = 0; j < getSizeY() + 1; ++j)
//	{
//		for (int i = 0; i < getSizeX(); ++i)
//		{
//			bool isEdge = j == 0 || j == getSizeY();
//			if (!isEdge && !(isBorder(i, j) || isBorder(i, j - 1)) && (phi(i, j) < 0.0 || phi(i, j - 1) < 0.0))
//			{ // Only add force to the liquid cells
//				doubleT xdistSqr = (x - j)*(x - j)*(j - x);
//				doubleT ydistSqr = (y - i)*(y - i)*(i - y);
//
//				set_uijHalfInd(i,j, uijHalfInd(i,j) + strength / xdistSqr);
//				set_vijHalfInd(i,j, vijHalfInd(i,j) + strength / ydistSqr);
//			}
//		}
//	}
}