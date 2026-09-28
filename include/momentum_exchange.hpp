#pragma once
#include "D2Q9.hpp"
#include "Grid.hpp"
#include "geometry_naca.hpp"
#include "aerodynamics.hpp"
#include <vector>

namespace lbm
{

    struct BoundaryLink
    {
        int x, y;         // Tọa độ nút chất lưu Boundary x_f
        int ff_x, ff_y;   // Tọa độ nút chất lưu lân cận x_ff = x_f - e_d
        bool has_ff;      // Cờ kiểm tra nút x_ff là chất lưu hợp lệ để nội suy Bouzidi
        int dir;          // Hướng đâm vào nút Solid (d)
        int opp_dir;      // Hướng ngược lại dội về chất lưu (OPP[d])
        double q;         // Tỷ lệ khoảng cách từ nút chất lưu đến vách cong thực tế q in (0, 1]
        double rx, ry;    // Vector từ pivot (0.25c) đến giao điểm thực trên bề mặt cánh
    };

    // Nội suy Bouzidi bậc 1 (Bouzidi, Firdaouss & Lallemand, 2001) từ phân bố sau va chạm f*
    inline double compute_ibb_reflected(const BoundaryLink &link, const Grid2D &grid)
    {
        const int n_idx = grid.idx(link.x, link.y);
        const double f_in = grid.f[link.dir][n_idx];
        if (!link.has_ff)
            return f_in;
        if (link.q < 0.5)
        {
            const int ff_idx = grid.idx(link.ff_x, link.ff_y);
            return 2.0 * link.q * f_in + (1.0 - 2.0 * link.q) * grid.f[link.dir][ff_idx];
        }
        return (1.0 / (2.0 * link.q)) * f_in + ((2.0 * link.q - 1.0) / (2.0 * link.q)) * grid.f[link.opp_dir][n_idx];
    }

    class MomentumExchange
    {
    public:
        double chord_lb;
        double x_pivot, y_pivot;
        double q_dyn_lb;
        double current_alpha_rad;
        std::vector<BoundaryLink> links;

        MomentumExchange(double chord, double u_inf, double xp, double yp);

        void build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, const NACAGeometry &geom);
        void build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, double alpha_deg);

        // Tính lực MEA cho biên nội suy IBB: F = sum (f_in + f_out) * e_dir
        AeroForce compute_force(const Grid2D &grid, bool is_post_collision = false) const;
    };

}
