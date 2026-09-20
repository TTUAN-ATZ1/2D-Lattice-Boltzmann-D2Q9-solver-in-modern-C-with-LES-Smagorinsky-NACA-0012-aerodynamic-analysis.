#pragma once
#include <vector>
#include <string>
#include "geometry_naca.hpp"
#include "Grid.hpp"
#include "momentum_exchange.hpp" // Sử dụng struct AeroForce

namespace lbm
{

    struct StencilPoint
    {
        int i0, j0, i1, j1;
        double w00, w10, w01, w11;
    };

    class Slide32Integrator
    {
    public:
        double chord_lb;
        double u_inf_lb;
        double q_dyn_lb;
        double Cd_friction;

        // Thông số thứ nguyên thực tế phục vụ xuất áp suất Pa
        double rho_real;
        double V_real;
        double P_inf_real;
        double q_real;

        std::vector<StencilPoint> stencils;

        Slide32Integrator(double chord, double u_inf, double Re,
                          double rho_r = 1.225, double V_r = 30.0, double P_inf_r = 101325.0);

        // Xây dựng ma trận nội suy song tuyến tính cho 120 điểm mốc
        void build_stencils(int nx, int ny, const std::vector<Point2D> &surface_points);

        // Tính tích phân lực áp suất Slide 3.2
        AeroForce integrate_pressure(const Grid2D &grid, const NACAGeometry &geom) const;

        // Xuất bảng số liệu 120 điểm mốc (X, Y, Cp, P_Pa, dFx, dFy) ra file CSV
        void export_csv(const std::string &filename,
                        const Grid2D &grid,
                        const NACAGeometry &geom,
                        double alpha_deg) const;
    };

}