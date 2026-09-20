#include "integrator.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace lbm
{

    Slide32Integrator::Slide32Integrator(double chord, double u_inf, double Re,
                                         double rho_r, double V_r, double P_inf_r)
        : chord_lb(chord), u_inf_lb(u_inf), rho_real(rho_r), V_real(V_r), P_inf_real(P_inf_r)
    {

        q_dyn_lb = 0.5 * 1.0 * (u_inf * u_inf) * chord_lb;
        Cd_friction = 2.0 * (1.328 / std::sqrt(Re));
        q_real = 0.5 * rho_real * (V_real * V_real);
    }

    void Slide32Integrator::build_stencils(int nx, int ny, const std::vector<Point2D> &surface_points)
    {
        stencils.clear();
        stencils.reserve(surface_points.size());

        for (const auto &pt : surface_points)
        {
            int i0 = static_cast<int>(std::floor(pt.x));
            int j0 = static_cast<int>(std::floor(pt.y));
            int i1 = std::min(i0 + 1, nx - 1);
            int j1 = std::min(j0 + 1, ny - 1);

            double dx = pt.x - i0;
            double dy = pt.y - j0;

            StencilPoint sp;
            sp.i0 = i0;
            sp.j0 = j0;
            sp.i1 = i1;
            sp.j1 = j1;
            sp.w00 = (1.0 - dx) * (1.0 - dy);
            sp.w10 = dx * (1.0 - dy);
            sp.w01 = (1.0 - dx) * dy;
            sp.w11 = dx * dy;

            stencils.push_back(sp);
        }
    }

    AeroForce Slide32Integrator::integrate_pressure(const Grid2D &grid, const NACAGeometry &geom) const
    {
        size_t Nb = geom.surface_points.size();
        double total_Fx = 0.0;
        double total_Fy = 0.0;
        double moment_y = 0.0;

        for (size_t k = 0; k < Nb; ++k)
        {
            const auto &sp = stencils[k];

            // 1. Nội suy mật độ rho tại điểm mốc k
            double rho_k = sp.w00 * grid.rho[grid.idx(sp.i0, sp.j0)] + sp.w10 * grid.rho[grid.idx(sp.i1, sp.j0)] + sp.w01 * grid.rho[grid.idx(sp.i0, sp.j1)] + sp.w11 * grid.rho[grid.idx(sp.i1, sp.j1)];

            // 2. Áp suất LBM: P = c_s^2 * (rho - 1.0)
            double P_k = CS2 * (rho_k - 1.0);

            // 3. Vi phân lực Slide 3.2: dF = -P * n * ds
            double dfx = -P_k * geom.normal[k].x * geom.ds[k];
            double dfy = -P_k * geom.normal[k].y * geom.ds[k];

            total_Fx += dfx;
            total_Fy += dfy;

            // Mô-men phục vụ tính tâm áp suất x_cp
            double x_rel = geom.raw_points[k].x / chord_lb;
            moment_y += x_rel * dfy;
        }

        AeroForce result;
        result.Fx = total_Fx;
        result.Fy = total_Fy;
        result.Cl = total_Fy / q_dyn_lb;
        result.Cdp = total_Fx / q_dyn_lb;
        result.Cd = std::max(result.Cdp, 0.0) + Cd_friction;

        if (std::abs(result.Cl) > 0.05)
        {
            double xcp = moment_y / total_Fy;
            result.x_cp = std::max(0.0, std::min(1.0, xcp));
        }
        else
        {
            result.x_cp = 0.25;
        }

        return result;
    }

    void Slide32Integrator::export_csv(const std::string &filename,
                                       const Grid2D &grid,
                                       const NACAGeometry &geom,
                                       double alpha_deg) const
    {
        std::ofstream fout(filename);
        if (!fout.is_open())
        {
            std::cerr << "Lỗi: Không thể tạo file " << filename << std::endl;
            return;
        }

        fout << "k,surface,X_norm,Y_norm,nx,ny,ds,Cp,P_real_Pa,dFx_N,dFy_N\n";

        size_t Nb = geom.surface_points.size();
        size_t mid_point = Nb / 2;

        for (size_t k = 0; k < Nb; ++k)
        {
            const auto &sp = stencils[k];
            double rho_k = sp.w00 * grid.rho[grid.idx(sp.i0, sp.j0)] + sp.w10 * grid.rho[grid.idx(sp.i1, sp.j0)] + sp.w01 * grid.rho[grid.idx(sp.i0, sp.j1)] + sp.w11 * grid.rho[grid.idx(sp.i1, sp.j1)];

            double P_k = CS2 * (rho_k - 1.0);
            double Cp = (2.0 * (rho_k - 1.0)) / (3.0 * u_inf_lb * u_inf_lb);
            double P_real = P_inf_real + Cp * q_real;

            double dfx = -P_k * geom.normal[k].x * geom.ds[k];
            double dfy = -P_k * geom.normal[k].y * geom.ds[k];

            // Đổi lực vi phân sang Newton thực tế
            double force_scale = q_real * 1.0 / q_dyn_lb;
            double dfx_real = dfx * force_scale;
            double dfy_real = dfy * force_scale;

            std::string surf = (k <= mid_point) ? "Upper" : "Lower";

            fout << k << ","
                 << surf << ","
                 << std::fixed << std::setprecision(6)
                 << (geom.raw_points[k].x / chord_lb) << ","
                 << (geom.raw_points[k].y / chord_lb) << ","
                 << geom.normal[k].x << ","
                 << geom.normal[k].y << ","
                 << geom.ds[k] << ","
                 << Cp << ","
                 << std::setprecision(2) << P_real << ","
                 << std::setprecision(5) << dfx_real << ","
                 << dfy_real << "\n";
        }

        fout.close();
    }

}