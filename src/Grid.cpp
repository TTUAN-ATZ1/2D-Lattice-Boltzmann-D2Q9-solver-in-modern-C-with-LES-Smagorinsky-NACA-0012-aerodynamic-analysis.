#include "Grid.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lbm
{
    // khoi tao luoi
    Grid2D::Grid2D(int nx_, int ny_) : nx(nx_), ny(ny_), size(0)
    {
        if (nx<3 || ny<3 || static_cast<long long>(nx)*ny>std::numeric_limits<int>::max())
            throw std::invalid_argument("Grid dimensions must be >=3 and fit an int index");
        size=nx*ny;
        for (int i = 0; i < Q; i++)
        {
            f[i].resize(size, 0.0);
            fn[i].resize(size, 0.0);
        }
        rho.resize(size, 1.0);
        ux.resize(size, 0.0);
        uy.resize(size, 0.0);
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
                ux[n] = std::numeric_limits<double>::quiet_NaN();
                uy[n] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}
