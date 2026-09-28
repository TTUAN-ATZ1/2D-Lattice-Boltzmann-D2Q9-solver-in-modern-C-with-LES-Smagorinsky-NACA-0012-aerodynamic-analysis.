#pragma once
#include <vector>
#include "D2Q9.hpp"

namespace lbm
{

    class Grid2D
    {
    public:
        int nx, ny; // nut luoi tai x,y
        int size;   // tong so o luoi

        std::vector<double> f[Q];  // ham phan bo hat truoc va cham
        std::vector<double> fn[Q]; // ham phan bo hat sau va cham

        // cac bien vi mo
        std::vector<double> rho;
        std::vector<double> ux;
        std::vector<double> uy;

        Grid2D(int nx, int ny);
        // chuyen toa do thanh 1D
        inline int idx(int x, int y) const
        {
            return y * nx + x;
        }
        // khoi tao truong dong ban dau
        void initialize(double u_inf);
        // ham tinh mat do rho va ux,uy tu ham fi
        void compute_macroscopic();
        // tinh phan bo can bang feq tai mot nut (inline de vector hoa SIMD)
        static inline double compute_feq(int i, double r, double u, double v)
        {
            const double cu = EX[i] * u + EY[i] * v;
            const double u_2 = u * u + v * v;
            return WEIGHTS[i] * r * (1.0 + INV_CS2 * cu + INV_2CS4 * (cu * cu) - INV_2CS2 * u_2);
        }
        static inline void compute_feq_all(double r, double u, double v, double feq_out[Q])
        {
            const double base = 1.0 - INV_2CS2 * (u * u + v * v);
            for (int i = 0; i < Q; ++i)
            {
                const double cu = EX[i] * u + EY[i] * v;
                feq_out[i] = WEIGHTS[i] * r * (base + INV_CS2 * cu + INV_2CS4 * (cu * cu));
            }
        }
    };
}