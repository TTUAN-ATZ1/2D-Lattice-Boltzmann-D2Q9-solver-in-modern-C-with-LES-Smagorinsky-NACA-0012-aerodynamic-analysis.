#pragma once
#include <vector>
#include "D2Q9.hpp"
#include "Grid.hpp"
#include "geometry_naca.hpp"
#include "les_smagorinsky.hpp"
#include "momentum_exchange.hpp"

namespace lbm
{

    class LBMSolver
    {
    public:
        int nx, ny;
        double u_inf;
        double tau_0;

        Grid2D grid;
        LESSmagorinsky les;

        std::vector<NodeType> node_type;
        std::vector<double> wall_dist;

        LBMSolver(int nx_, int ny_, double u_inf_, double Re, double chord_lb, double Cs = 0.16);

        // Cập nhật hình học cánh mới khi đổi góc tấn alpha
        void set_geometry(const NACAGeometry &geom);

        // Chạy 1 bước thời gian LBM
        void step();

        // Chạy nhiều bước lặp
        void run_steps(int num_steps);
    };

}