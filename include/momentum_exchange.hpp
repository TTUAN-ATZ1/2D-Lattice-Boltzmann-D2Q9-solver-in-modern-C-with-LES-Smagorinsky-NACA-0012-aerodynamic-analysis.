#pragma once
#include <vector>
#include "D2Q9.hpp"
#include "Grid.hpp"
#include "geometry_naca.hpp"

namespace lbm
{

    // Cấu trúc 1 liên kết cắt qua biên cánh
    struct BoundaryLink
    {
        int x, y;      // Tọa độ nút chất lưu (NODE_BOUNDARY)
        int dir;       // Hướng hạt c_i đâm vào ruột cánh (NODE_SOLID)
        int opp_dir;   // Hướng phản xạ dội ngược lại (c_{OPP[dir]})
        double rx, ry;  // Cánh tay đòn từ tâm quay tới trung điểm liên kết
        double x_chord; // Tọa độ thực dọc trục dây cung từ LE [0, chord]
    };

    struct AeroForce
    {
        double Fx, Fy;
        double Cl, Cd, Cdp;
        double Mz;
        double x_cp;
    };

    class MomentumExchange
    {
    public:
        double chord_lb;
        double q_dyn_lb;
        double Cd_friction;
        double x_pivot, y_pivot;
        double current_alpha_rad;

        std::vector<BoundaryLink> links;

        MomentumExchange(double chord, double u_inf, double Re, double xp, double yp);

        // Tìm và lưu toàn bộ các liên kết cắt qua biên cánh
        void build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, double alpha_deg = 0.0);

        // Tính lực và hệ số khí động tức thời tại bước lặp hiện tại
        AeroForce compute_force(const Grid2D &grid) const;
    };

}