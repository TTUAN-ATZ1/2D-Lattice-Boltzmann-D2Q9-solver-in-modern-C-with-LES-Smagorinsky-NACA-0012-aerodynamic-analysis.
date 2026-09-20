#include "momentum_exchange.hpp"
#include <cmath>

namespace lbm
{

    MomentumExchange::MomentumExchange(double chord, double u_inf, double Re, double xp, double yp)
        : chord_lb(chord), x_pivot(xp), y_pivot(yp)
    {

        // Áp suất động LBM: q = 0.5 * rho * u^2 * chord
        q_dyn_lb = 0.5 * 1.0 * (u_inf * u_inf) * chord_lb;

        // Lực cản ma sát Blasius 2 mặt ở Re = 100,000
        Cd_friction = 2.0 * (1.328 / std::sqrt(Re));
    }

    void MomentumExchange::build_boundary_links(int nx, int ny, const std::vector<NodeType> &node_type, double alpha_deg)
    {
        links.clear();
        current_alpha_rad = alpha_deg * 3.14159265358979323846 / 180.0;
        double cos_a = std::cos(current_alpha_rad);
        double sin_a = std::sin(current_alpha_rad);

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

                        // Nếu hướng d đâm vào nút rắn SOLID
                        if (node_type[neighbor_idx] == NODE_SOLID)
                        {
                            BoundaryLink link;
                            link.x = i;
                            link.y = j;
                            link.dir = d;
                            link.opp_dir = OPP[d];

                            // Trung điểm của liên kết làm điểm đặt lực
                            double mid_x = i + 0.5 * EX[d];
                            double mid_y = j + 0.5 * EY[d];
                            link.rx = mid_x - x_pivot;
                            link.ry = mid_y - y_pivot;

                            // Chiếu ngược về hệ tọa độ cánh để lấy tọa độ dọc dây cung từ LE
                            double dx_chord = link.rx * cos_a - link.ry * sin_a;
                            link.x_chord = 0.25 * chord_lb + dx_chord;

                            links.push_back(link);
                        }
                    }
                }
            }
        }
    }

    AeroForce MomentumExchange::compute_force(const Grid2D &grid) const
    {
        double total_Fx = 0.0;
        double total_Fy = 0.0;
        double total_Mz = 0.0;
        double total_Fn = 0.0;
        double moment_chord = 0.0;

        double cos_a = std::cos(current_alpha_rad);
        double sin_a = std::sin(current_alpha_rad);

        int num_links = static_cast<int>(links.size());

#pragma omp parallel for reduction(+ : total_Fx, total_Fy, total_Mz, total_Fn, moment_chord)
        for (int k = 0; k < num_links; ++k)
        {
            const auto &link = links[k];
            int n_idx = grid.idx(link.x, link.y);

            // Với vách tĩnh half-way bounce-back: f_bar(x_f, t+1) = f_i^*(x_f, t)
            // Lực tác dụng lên vách: F = 2 * c_i * f_out
            double f_out = grid.f[link.opp_dir][n_idx];
            double momentum_transfer = 2.0 * f_out;
            double fx = EX[link.dir] * momentum_transfer;
            double fy = EY[link.dir] * momentum_transfer;

            total_Fx += fx;
            total_Fy += fy;
            total_Mz += (link.rx * fy - link.ry * fx);

            // Lực pháp tuyến bề mặt cánh: Fn = Fy * cos(alpha) - Fx * sin(alpha)
            double fn = fy * cos_a - fx * sin_a;
            total_Fn += fn;
            moment_chord += link.x_chord * fn;
        }

        AeroForce result;
        result.Fx = total_Fx;
        result.Fy = total_Fy;
        result.Mz = total_Mz;

        // Quy đổi sang hệ số không thứ nguyên
        result.Cl = total_Fy / q_dyn_lb;
        result.Cdp = total_Fx / q_dyn_lb;
        result.Cd = std::max(result.Cdp, 0.0) + Cd_friction;

        // Vị trí tâm áp suất x_cp (tính từ mép trước LE, gốc quay tại 0.25c)
        // Với cánh đối xứng khi |Cl| < 0.05, tâm áp suất theo lý thuyết khí động học nằm tại tâm khí động 0.25c
        if (std::abs(result.Cl) > 0.05)
        {
            double xcp = 0.25 + (total_Mz / (total_Fy * chord_lb));
            result.x_cp = std::max(0.0, std::min(1.0, xcp));
        }
        else
        {
            result.x_cp = 0.25;
        }

        return result;
    }

}