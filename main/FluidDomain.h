//
// Created by Niklaus on 12.11.19.
// This class holds all grid data, grid convenience function (e.g. interpolate, get half_indexed..)
//https://github.com/kbladin/Fluid_Simulation This project helped us a lot to get started, and some function are based from there.
//

#ifndef MAIN_FLUIDDOMAIN_H
#define MAIN_FLUIDDOMAIN_H

#include "MarkerParticle.h"
#include <algorithm>
#include <random>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "GuiData.h"
#include "ParticleLevelSet.h"
#include "GuiData.h"

enum CellType
{
    LIQUID, AIR, SOLID
};

class FluidDomain {

public:
    FluidDomain(int m_size_x, int m_size_y, int m_size_z, doubleT m_length_x);
    ~FluidDomain();

    void classifyCells();
    void buildMesh();
    void buildGrid();
    void updateEdges(doubleT scale);
    void getMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const;
    void getColors(Eigen::MatrixXf& C, bool normalize=false) const;
    void addSource(doubleT xmin, doubleT xmax, doubleT ymin, doubleT ymax, doubleT randomness = 0.2);
	void setTemperature(doubleT T, doubleT xmin, doubleT xmax, doubleT ymin, doubleT ymax);
	void buildTemperatureMatrix();
	void clearPressure() { m_p.setZero(); }

    //write buffer to real Array2d
    void assignVelocityTmp() {  m_u = m_u_tmp; m_v = m_v_tmp; };
    void updatePreviousVelocityBuffer() { m_u_prev = m_u; m_v_prev = m_v; }
    void updateVelocityDiffBuffer()
    {
        for (int j = 0; j < m_size_y; ++j)
        {
            for (int i = 0; i < m_size_x; ++i)
            {
                m_u_diff(i, j) = m_u(i, j) - m_u_prev(i, j);
                m_v_diff(i, j) = m_v(i, j) - m_v_prev(i, j);
            }
        }
    }

    const Eigen::MatrixXd& s() { return m_start; }
    const Eigen::MatrixXd& e() { return m_end; }
    const Eigen::MatrixXd& vs() { return m_vs; }
    const Eigen::MatrixXd& ve() { return m_ve; }
    const Eigen::MatrixXd& vc() { return m_vc; }


    /*
     * access function to make coding easier. These function handle the half indices and interpolation.
     * we have two Array2d for u-velocity and v-velocity, where tmp is a buffer to store intermediate changes.
     */
#pragma region ConvenienceFunctions

    template <typename SCALAR>
    inline SCALAR CLAMP(SCALAR d, SCALAR min, SCALAR max) const {
        const SCALAR t = d < min ? min : d;
        return t > max ? max : t;
    }

    inline int doubleT_to_floor(doubleT x) const {
        int a = (int)std::floor(x);
        while (a > x) a--;
        while(a+1 < x) a++;
        return a;
    }

    inline doubleT uij(int i, int j) const { return (m_u(i,j) + m_u(i+1,j)) * 0.5; };
    inline doubleT vij(int i, int j) const { return (m_v(i,j) + m_v(i,j+1)) * 0.5; };
    inline doubleT uijHalfInd(int i, int j) const { return m_u(i,j); };
    inline doubleT vijHalfInd(int i, int j) const { return m_v(i,j); };
    inline doubleT uij_div(int i, int j) const { return (m_u(i + 1, j) - m_u(i, j)) / m_delta_x; }; //divergence
    inline doubleT vij_div(int i, int j) const { return (m_v(i, j + 1) - m_v(i, j)) / m_delta_y; }; //divergence
    inline doubleT pij(int i, int j) const { return m_p(i,j); };
    inline doubleT Tij(int i, int j) const { return m_T(i,j); };
    inline doubleT viscosity(int i, int j) const { return m_viscosity(i,j); };

    inline CellType cellType(int i, int j) const { return m_cell_type(i,j); };
    inline doubleT phi(int i, int j) const { return liquidSDF(i,j); }
    ParticleLevelSet& getLiquidSDF() { return liquidSDF; }

    inline doubleT u_interpolate(doubleT x, doubleT y) const { return valueInterpolated(x, y - m_delta_y * 0.5, m_u); }; //0.5 bc of half index
    inline doubleT v_interpolate(doubleT x, doubleT y) const { return valueInterpolated(x - m_delta_x * 0.5, y, m_v); };


    inline void set_uijHalfInd(int i, int j, doubleT u) { m_u(i, j) = u; };
    inline void set_vijHalfInd(int i, int j, doubleT v) { m_v(i, j) = v; };
    inline void set_uijHalfInd_tmp(int i, int j, doubleT u) { m_u_tmp(i, j) = u; };
    inline void set_vijHalfInd_tmp(int i, int j, doubleT v) { m_v_tmp(i, j) = v; };
    inline void set_pij(int i,int j, doubleT p) { m_p(i, j) = p; };
    inline void set_Tij(int i,int j, doubleT t) { m_T(i, j) = t; };
    inline void increaseCounter() {counter++;};
    inline void setViscosity(int i, int j, doubleT value) { m_viscosity(i,j) = value; };
    inline void setViscosity(float value) { m_viscosity.fill(value); };

    inline void u_tmp_addToInterpolated(doubleT x, doubleT y, doubleT u) { addToValueInterpolated(x, y - 0.5 * m_delta_y, u, m_u_tmp); };
    inline void v_tmp_addToInterpolated(doubleT x, doubleT y, doubleT v) { addToValueInterpolated(x - 0.5 * m_delta_x, y, v, m_v_tmp); };

    inline doubleT u_diff_interpolate(doubleT x, doubleT y) { return valueInterpolated(x, y - m_delta_y * 0.5, m_u_diff); };
    inline doubleT v_diff_interpolate(doubleT x, doubleT y) { return valueInterpolated(x - m_delta_x * 0.5, y, m_v_diff); };

    //interpolates value of world location x,y based on grid values.
    inline doubleT valueInterpolated(doubleT x, doubleT y, const MatrixXt& grid) const
    {
        // Calculate indices
        int i = doubleT_to_floor(x / m_delta_x);
        int j = doubleT_to_floor(y / m_delta_y);
        doubleT i_frac = x / m_delta_x - i;
        doubleT j_frac = y / m_delta_y - j;

        i = CLAMP(i, 0, m_size_x - 1);
        j = CLAMP(j, 0, m_size_y - 1);
        int i_plus1 = CLAMP(i + 1, 0, m_size_x - 1);
        int j_plus1 = CLAMP(j + 1, 0, m_size_y - 1);

        // First interpolate in x, then in y
        doubleT value_00 = grid(i, j);
        doubleT value_10 = grid(i_plus1, j);
        doubleT value_01 = grid(i, j_plus1);
        doubleT value_11 = grid(i_plus1, j_plus1);

        // Interpolate x
        doubleT value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
        doubleT value_1 = (1 - i_frac) * value_01 + i_frac * value_11;

        // Interpolate y
        doubleT value = (1 - j_frac) * value_0 + j_frac * value_1;
        return value;
    }

    //Takes a value at world location x and y and distributes it over the grid
    inline void addToValueInterpolated(doubleT x, doubleT y, doubleT value, MatrixXt& grid)
    {
        // Calculate indices
        int i = doubleT_to_floor(x / m_delta_x);
        int j = doubleT_to_floor(y / m_delta_y);
        int i_plus1 = i + 1;
        int j_plus1 = j + 1;
        doubleT i_frac = x / m_delta_x - i;
        doubleT j_frac = y / m_delta_y - j;

        // Border cases
        i = CLAMP(i, 0, m_size_x - 1);
        j = CLAMP(j, 0, m_size_y - 1);
        i_plus1 = CLAMP(i_plus1, 0, m_size_x - 1);
        j_plus1 = CLAMP(j_plus1, 0, m_size_y - 1);

        // Spread in y
        doubleT value_0 = (1 - j_frac) * value;
        doubleT value_1 = j_frac * value;

        // Spread in x
        doubleT value_00 = (1 - i_frac) * value_0;
        doubleT value_10 = i_frac * value_0;
        doubleT value_01 = (1 - i_frac) * value_1;
        doubleT value_11 = i_frac * value_1;

        // Write data
        grid(i, j) += value_00;
        grid(i_plus1, j) += value_10;
        grid(i, j_plus1) += value_01;
        grid(i_plus1, j_plus1) += value_11;
    }

    doubleT bilinearInterpolateLiquid(doubleT x, doubleT y) {
        // Calculate indices
        int i = doubleT_to_floor(x / m_delta_x);
        int j = doubleT_to_floor(y / m_delta_y);
        doubleT i_frac = x / m_delta_x - i;
        doubleT j_frac = y / m_delta_y - j;

        i = CLAMP(i, 0, m_size_x - 1);
        j = CLAMP(j, 0, m_size_y - 1);
        int i_plus1 = CLAMP(i + 1, 0, m_size_x - 1);
        int j_plus1 = CLAMP(j + 1, 0, m_size_y - 1);

        // First interpolate in x, then in y
        doubleT value_00 = cellType(i,j) == LIQUID ? 1.0 : 0.0;
        doubleT value_10 = cellType(i_plus1,j) == LIQUID ? 1.0 : 0.0;
        doubleT value_01 = cellType(i,j_plus1) == LIQUID ? 1.0 : 0.0;
        doubleT value_11 = cellType(i_plus1,j_plus1) == LIQUID ? 1.0 : 0.0;

        // Interpolate x
        doubleT value_0 = (1 - i_frac) * value_00 + i_frac * value_10;
        doubleT value_1 = (1 - i_frac) * value_01 + i_frac * value_11;

        // Interpolate y
        doubleT value = (1 - j_frac) * value_0 + j_frac * value_1;
        return value;
    }

    bool isBorder(int i, int j) { return i == 0 || i == m_size_x - 1 || j == 0 || j == m_size_y - 1; }
#pragma endregion ConvenienceFunctions

#pragma region MarkerParticleSet
    void addParticle(MarkerParticle p) 	{ m_particles.push_back(p); };
    void reserve(int particle_count)	{ m_particles.reserve(particle_count); };
    void clear() 						{ m_particles.clear(); };
    void advect(doubleT dt);
    void advectAndEnsureOutsideObstacles(doubleT dt);
    std::vector<MarkerParticle>& getParticles() { return m_particles; };
	MatrixXt& getPressure() { return m_p; };
	MatrixXt& getTemperature() { return m_T; };
	Eigen::ConjugateGradient<Eigen::SparseMatrix<doubleT>, Eigen::Lower|Eigen::Upper>& getTemperatureSolver() { return m_TSolver; };

#pragma endregion MarkerParticleSet

    void reset()
    {
        m_cell_type.fill(AIR);
	    m_viscosity.fill(GuiData::instance->m_viscosity);
        m_u.setZero();
        m_v.setZero();
        m_u_tmp.setZero();
        m_v_tmp.setZero();
        m_u_diff.setZero();
        m_v_diff.setZero();
        m_u_prev.setZero();
        m_v_prev.setZero();
        m_particles.clear();
        m_p.setZero();
        m_T.fill(GuiData::instance->m_baseTemperature);
        counter = 1;
        buildTemperatureMatrix();
    }

	inline int getSizeX() const { return m_size_x; }
	inline int getSizeY() const { return m_size_y; }
	inline int getLengthX() const { return m_length_x; }
	inline int getLengthY() const { return m_length_y; }
    inline doubleT getDeltaX() const { return m_delta_x; }
    inline doubleT getDeltaY() const { return m_delta_y; }
    inline int getNumOfFluidCells() const { return m_num_of_fluid_cells; };
    inline int getCounter() const { return counter; }


	using MatrixXb = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
	void readBrushTexture(MatrixXb& out, std::string imageFileName, int x_res, int y_res, char drawChar = 'Q');
	MatrixXb arrowBrushImage;
	MatrixXb bunnyBrushImage;

	void AddPointForce(doubleT x, doubleT y);

private:
    int m_size_x, m_size_y, m_size_z;
    doubleT m_length_x, m_length_y, m_length_z;
    doubleT m_delta_x, m_delta_y, m_delta_z;
    int m_num_of_fluid_cells;

    std::vector<MarkerParticle> m_particles;

    Eigen::Matrix<CellType, Eigen::Dynamic, Eigen::Dynamic> m_cell_type;
    MatrixXt m_u;
	MatrixXt m_v;
	MatrixXt m_u_tmp; //buffer to store velocity changes
	MatrixXt m_v_tmp; //buffer to store velcity changes

	MatrixXt m_u_prev;
	MatrixXt m_v_prev;

	MatrixXt m_u_diff;
	MatrixXt m_v_diff;

    int counter;
	MatrixXt m_p; // pressure
	MatrixXt m_T; // temperature
    Eigen::SparseMatrix<doubleT> m_TImplMatrix;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<doubleT>, Eigen::Lower|Eigen::Upper> m_TSolver;

	//Viscosity
	MatrixXt m_viscosity;

    //for render
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;

	Eigen::MatrixXd m_start;
	Eigen::MatrixXd m_end;
	Eigen::MatrixXd m_vs;
	Eigen::MatrixXd m_ve;
	Eigen::MatrixXd m_vc;

	std::mt19937 rng;

	//Particle SDF
    ParticleLevelSet liquidSDF;
    doubleT particleRadius = 0.0;
};

#endif //MAIN_FLUIDDOMAIN_H
