#include "Grid.hpp"
#include <cmath>

namespace lbm
{
    // khoi tao luoi
    Grid2D::Grid2D(int nx_, int ny_) : nx(nx_), ny(ny_), size(nx_ * ny_)
    {
        for (int i = 0; i < Q; i++)
        {
            f[i].resize(size, 0.0);
            fn[i].resize(size, 0.0);
        }
        rho.resize(size, 1.0);
        ux.resize(size, 0.0);
        uy.resize(size, 0.0);
    }
    // tinh ham phan bo can bang
    double Grid2D::compute_feq(int i, double r, double u, double v)
    {
        double cu = EX[i] * u + EY[i] * v;
        double u_2 = u * u + v * v;
        return WEIGHTS[i] * r * (1.0 + INV_CS2 * cu + INV_2CS4 * (cu * cu) - INV_2CS2 * u_2);
    }
    // khoi tao gia tri dong vao
    void Grid2D::initialize(double u_inf)
    {
#pragma omp parallel for
        for (int n = 0; n < size; ++n)
        {
            rho[n] = 1.0;
            ux[n] = u_inf;
            uy[n] = 0.0;
            for (int i = 0; i < Q; ++i)
            {
                double feq_val = compute_feq(i, 1.0, u_inf, 0.0);
                f[i][n] = feq_val;
                fn[i][n] = feq_val;
            }
        }
    }
    // tinh cac dai luong vi mo
    void Grid2D::compute_macroscopic()
    {
#pragma omp parallel for
        for (int n = 0; n < size; ++n)
        {
            double r = 0.0;
            double rux = 0.0;
            double ruy = 0.0;
            for (int i = 0; i < Q; ++i)
            {
                double fi = f[i][n];
                r += fi;
                rux += fi * EX[i];
                ruy += fi * EY[i];
            }
            rho[n] = r;
            if (r > 1e-12)
            {
                ux[n] = rux / r;
                uy[n] = ruy / r;
            }
            else
            {
                ux[n] = 0.0;
                uy[n] = 0.0;
            }
        }
    }
}