#include "vtk_exporter.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>

namespace lbm
{

    VTKExporter::VTKExporter(double u_inf, double rho_r, double V_r, double P_inf_r)
        : u_inf_lb(u_inf), rho_real(rho_r), V_real(V_r), P_inf_real(P_inf_r)
    {
        q_real = 0.5 * rho_real * (V_real * V_real);
    }

    void VTKExporter::export_field_2d(const std::string &filename,
                                      const Grid2D &grid,
                                      const std::vector<NodeType> &node_type) const
    {
        std::ofstream fout(filename);
        if (!fout.is_open())
        {
            std::cerr << "Lỗi: Không thể tạo file " << filename << std::endl;
            return;
        }

        int nx = grid.nx;
        int ny = grid.ny;
        int total = nx * ny;

        fout << "# vtk DataFile Version 3.0\n";
        fout << "LBM 2D Flow Field NACA0012\n";
        fout << "ASCII\n";
        fout << "DATASET STRUCTURED_POINTS\n";
        fout << "DIMENSIONS " << nx << " " << ny << " 1\n";
        fout << "ORIGIN 0 0 0\n";
        fout << "SPACING 1 1 1\n";
        fout << "POINT_DATA " << total << "\n";

        // 1. Áp suất thực tế P (Pa)
        fout << "SCALARS Pressure_Pa float 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int idx = j * nx + i;
                double r = grid.rho[idx];
                double Cp = (2.0 * (r - 1.0)) / (3.0 * u_inf_lb * u_inf_lb);
                double P_real = P_inf_real + Cp * q_real;
                fout << std::fixed << std::setprecision(2) << P_real << "\n";
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
                double vel_real = mag * (V_real / u_inf_lb);
                fout << std::fixed << std::setprecision(4) << vel_real << "\n";
            }
        }

        // 3. Trường độ xoáy Vorticity (1/s)
        fout << "SCALARS Vorticity float 1\nLOOKUP_TABLE default\n";
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
                double vort = duy_dx - dux_dy;
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
                double u_real = grid.ux[idx] * (V_real / u_inf_lb);
                double v_real = grid.uy[idx] * (V_real / u_inf_lb);
                fout << std::fixed << std::setprecision(3) << u_real << " " << v_real << " 0.0\n";
            }
        }

        fout.close();
    }

    void VTKExporter::export_airfoil_surface(const std::string &filename,
                                             const NACAGeometry &geom,
                                             const std::vector<double> &cp_values) const
    {
        std::ofstream fout(filename);
        if (!fout.is_open())
            return;

        int N = static_cast<int>(geom.surface_points.size());

        fout << "# vtk DataFile Version 3.0\n";
        fout << "NACA0012 Surface Profile\n";
        fout << "ASCII\n";
        fout << "DATASET POLYDATA\n";

        // Danh sách điểm
        fout << "POINTS " << N << " float\n";
        for (int i = 0; i < N; ++i)
        {
            fout << std::fixed << std::setprecision(4)
                 << geom.surface_points[i].x << " "
                 << geom.surface_points[i].y << " 0.0\n";
        }

        // Đường viền khép kín
        fout << "LINES 1 " << (N + 1) << "\n";
        fout << N;
        for (int i = 0; i < N; ++i)
        {
            fout << " " << i;
        }
        fout << "\n";

        // Dữ liệu gán lên điểm
        fout << "POINT_DATA " << N << "\n";
        fout << "SCALARS Cp float 1\nLOOKUP_TABLE default\n";
        for (int i = 0; i < N; ++i)
        {
            double cp = (i < static_cast<int>(cp_values.size())) ? cp_values[i] : 0.0;
            fout << std::fixed << std::setprecision(4) << cp << "\n";
        }

        fout.close();
    }

} // namespace lbm