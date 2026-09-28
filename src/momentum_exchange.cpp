#include "momentum_exchange.hpp"
#include <cmath>
#include <algorithm>

namespace lbm
{

    MomentumExchange::MomentumExchange(double chord, double u_inf, double xp, double yp)
        : chord_lb(chord), x_pivot(xp), y_pivot(yp), current_alpha_rad(0.0)
    {
        q_dyn_lb = 0.5 * 1.0 * (u_inf * u_inf) * chord_lb;
        if (!std::isfinite(q_dyn_lb) || q_dyn_lb<=0)
            throw std::invalid_argument("Invalid MEA normalization");
    }

    void MomentumExchange::build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, const NACAGeometry &geom)
    {
        links.clear();
        if (node_type.size()!=static_cast<size_t>(nx)*ny)
            throw std::invalid_argument("MEA mask size mismatch");
        current_alpha_rad = geom.alpha_deg * PI / 180.0;
        const size_t N = geom.surface_points.size();

        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                int n_idx = j * nx + i;
                if (node_type[n_idx] == NODE_BOUNDARY)
                {
                    for (int d = 1; d < Q; ++d)
                    {
                        int ni = i + EX[d];
                        int nj = j + EY[d];
                        int neighbor_idx = nj * nx + ni;

                        if (node_type[neighbor_idx] == NODE_SOLID)
                        {
                            BoundaryLink link;
                            link.x = i;
                            link.y = j;
                            link.dir = d;
                            link.opp_dir = OPP[d];

                            int ffx = i - EX[d];
                            int ffy = j - EY[d];
                            link.ff_x = ffx;
                            link.ff_y = ffy;
                            link.has_ff = (ffx >= 0 && ffx < nx && ffy >= 0 && ffy < ny
                                           && node_type[ffy * nx + ffx] != NODE_SOLID);

                            // Tìm tỷ lệ khoảng cách thực q in (0, 1] từ nút chất lưu (i, j) đến đa giác cánh
                            double q_min = 2.0;
                            const double ex = EX[d];
                            const double ey = EY[d];
                            for (size_t k = 0; k < N; ++k)
                            {
                                const auto &a = geom.surface_points[k];
                                const auto &b = geom.surface_points[(k + 1) % N];
                                const double vx = b.x - a.x;
                                const double vy = b.y - a.y;
                                const double wx = a.x - i;
                                const double wy = a.y - j;
                                const double denom = ex * vy - ey * vx;
                                if (std::abs(denom) > 1e-14)
                                {
                                    const double q_val = (wx * vy - wy * vx) / denom;
                                    const double s_val = (wx * ey - wy * ex) / denom;
                                    if (s_val >= -1e-9 && s_val <= 1.0 + 1e-9 && q_val >= -1e-9 && q_val <= 1.0 + 1e-9)
                                    {
                                        q_min = std::min(q_min, std::clamp(q_val, 1e-4, 1.0));
                                    }
                                }
                            }
                            link.q = (q_min <= 1.0) ? q_min : 0.5;
                            link.rx = (i + link.q * ex) - x_pivot;
                            link.ry = (j + link.q * ey) - y_pivot;
                            links.push_back(link);
                        }
                    }
                }
            }
        }
    }

    void MomentumExchange::build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, double alpha_deg)
    {
        NACAGeometry geom(chord_lb, alpha_deg, x_pivot, y_pivot);
        build_boundary_links(nx, ny, node_type, geom);
    }

    AeroForce MomentumExchange::compute_force(const Grid2D &grid, bool is_post_collision) const
    {
        double total_Fx = 0.0;
        double total_Fy = 0.0;
        double total_Mz = 0.0;

        int num_links = static_cast<int>(links.size());

#pragma omp parallel for reduction(+ : total_Fx, total_Fy, total_Mz)
        for (int k = 0; k < num_links; ++k)
        {
            const auto &link = links[k];
            int n_idx = grid.idx(link.x, link.y);

            // Trước swap (is_post_collision = true): f_in nằm ở f[dir], f_out nằm ở fn[opp_dir]
            // Sau swap (is_post_collision = false): f_in nằm ở fn[dir], f_out nằm ở f[opp_dir]
            double f_in  = is_post_collision ? grid.f[link.dir][n_idx]      : grid.fn[link.dir][n_idx];
            double f_out = is_post_collision ? grid.fn[link.opp_dir][n_idx] : grid.f[link.opp_dir][n_idx];

            double momentum_transfer = f_in + f_out;

            double fx = EX[link.dir] * momentum_transfer;
            double fy = EY[link.dir] * momentum_transfer;

            total_Fx += fx;
            total_Fy += fy;
            total_Mz += (link.rx * fy - link.ry * fx);
        }

        return aerodynamic_force(total_Fx,total_Fy,total_Mz,current_alpha_rad,
                                 chord_lb,q_dyn_lb);
    }

}
