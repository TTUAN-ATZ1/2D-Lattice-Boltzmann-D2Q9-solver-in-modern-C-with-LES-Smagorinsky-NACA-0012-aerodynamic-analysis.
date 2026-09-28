#include "vtk_exporter.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace lbm
{

    VTKExporter::VTKExporter(const PhysicalUnits& physical_units) : units(physical_units) {}

    void VTKExporter::export_field_2d(const std::string &filename,
                                      const Grid2D &grid,
                                      const std::vector<NodeType> &node_type) const
    {
        if (node_type.size()!=static_cast<size_t>(grid.size)) throw std::invalid_argument("VTK mask size mismatch");
        std::ofstream fout;
        fout.exceptions(std::ios::failbit|std::ios::badbit);
        fout.open(filename);
        fout << std::setprecision(17);

        int nx = grid.nx;
        int ny = grid.ny;
        int total = nx * ny;

        fout << "# vtk DataFile Version 3.0\n";
        fout << "LBM 2D Flow Field NACA0012\n";
        fout << "ASCII\n";
        fout << "DATASET STRUCTURED_POINTS\n";
        fout << "DIMENSIONS " << nx << " " << ny << " 1\n";
        fout << "ORIGIN 0 0 0\n";
        fout << "SPACING " << units.dx << ' ' << units.dx << " 1\n";
        fout << "POINT_DATA " << total << "\n";

        // 1. Áp suất thực tế P (Pa)
        fout << "SCALARS Pressure_Pa float 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx = j * nx + i;
                double r = grid.rho[idx];
                double Cp = (2.0 * (r - 1.0)) / (3.0 * units.u_lb * units.u_lb);
                double P_real = units.pressure + Cp * units.q();
                fout << std::scientific << std::setprecision(12) << P_real << "\n";
            }
        }

        // 2. Độ lớn vận tốc thực tế (m/s)
        fout << "SCALARS Velocity_Mag float 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx = j * nx + i;
                double mag = std::sqrt(grid.ux[idx] * grid.ux[idx] + grid.uy[idx] * grid.uy[idx]);
                double vel_real = mag * (units.velocity / units.u_lb);
                fout << std::fixed << std::setprecision(4) << vel_real << "\n";
            }
        }

        // 3. Trường độ xoáy Vorticity (1/s)
        fout << "SCALARS Vorticity_per_s float 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int i_prev = std::max(0, i - 1);
                int i_next = std::min(nx - 1, i + 1);
                int j_prev = std::max(0, j - 1);
                int j_next = std::min(ny - 1, j + 1);

                double duy_dx = (grid.uy[j * nx + i_next] - grid.uy[j * nx + i_prev]) / (i_next - i_prev);
                double dux_dy = (grid.ux[j_next * nx + i] - grid.ux[j_prev * nx + i]) / (j_next - j_prev);
                double vort = (duy_dx - dux_dy)/units.dt;
                fout << std::fixed << std::setprecision(5) << vort << "\n";
            }
        }

        // 4. Mặt nạ ruột cánh Airfoil_Mask (Solid = 1, Fluid/Boundary = 0)
        fout << "SCALARS Airfoil_Mask int 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx = j * nx + i;
                int mask = (node_type[idx] == NODE_SOLID) ? 1 : 0;
                fout << mask << "\n";
            }
        }

        // 5. Vector vận tốc 2D (u, v, 0)
        fout << "VECTORS Velocity float\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx = j * nx + i;
                double u_real = grid.ux[idx] * (units.velocity / units.u_lb);
                double v_real = grid.uy[idx] * (units.velocity / units.u_lb);
                fout << std::fixed << std::setprecision(3) << u_real << " " << v_real << " 0.0\n";
            }
        }

        fout.close();
    }

} // namespace lbm
