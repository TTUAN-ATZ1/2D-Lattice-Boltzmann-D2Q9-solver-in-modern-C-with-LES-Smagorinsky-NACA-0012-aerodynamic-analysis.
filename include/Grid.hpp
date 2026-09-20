#pragma once
#include <vector>
#include "D2Q9.hpp"

namespace lbm
{

    class Grid2D
    {
    public:
        int nx, ny; // nut luoi taix,y
        int size;   // o luoi

        std::vector<double> f[Q];  // ham phan bo hat truoc va cham
        std::vector<double> fn[Q]; // ham phan bo hat sau va cham

        // cac bien vi mo
        std::vector<double> rho;
        std::vector<double> ux;
        std::vector<double> uy;

        Grid2D(int nx, int ny);
        // chuyen toa do thanh 1D;
        inline int idx(int x, int y) const
        {
            return y * nx + x;
        }
        // khoi tao truong dong ban dau
        void initialize(double u_inf);
        // ham tinh mat do rho va ux,uy tu ham fi
        void compute_macroscopic();
        // tinh phan bo can bang feq tai mot nut
        static double compute_feq(int i, double rho_val, double u_val, double v_val);
    };
}