#include "solver.hpp"
#include <cmath>
#include <algorithm>

namespace lbm
{

    LBMSolver::LBMSolver(int nx_, int ny_, double u_inf_, double Re, double chord_lb, double Cs)
        : nx(nx_), ny(ny_), u_inf(u_inf_),
          grid(nx_, ny_),
          les(Cs, 3.0 * (u_inf_ * chord_lb / Re) + 0.5, (u_inf_ * chord_lb / Re), u_inf_, Re)
    {

        tau_0 = 3.0 * (u_inf * chord_lb / Re) + 0.5;
        grid.initialize(u_inf);
    }

    void LBMSolver::set_geometry(const NACAGeometry &geom)
    {
        std::vector<NodeType> old_types = node_type;
        geom.classify_grid(nx, ny, node_type, wall_dist);

        // Khử shock rác số học: Nếu nút từ trong ruột cánh (SOLID) vừa lộ ra chất lưu
        if (!old_types.empty())
        {
            int total = nx * ny;
            for (int n = 0; n < total; ++n)
            {
                if (old_types[n] == NODE_SOLID && node_type[n] != NODE_SOLID)
                {
                    grid.rho[n] = 1.0;
                    grid.ux[n] = u_inf;
                    grid.uy[n] = 0.0;
                    for (int d = 0; d < Q; ++d)
                    {
                        double feq = Grid2D::compute_feq(d, 1.0, u_inf, 0.0);
                        grid.f[d][n] = feq;
                        grid.fn[d][n] = feq;
                    }
                }
            }
        }
    }

    void LBMSolver::step()
    {
        int total = nx * ny;

        // 1. Thu thập biến vĩ mô rho, ux, uy
        grid.compute_macroscopic();

// 2. Va chạm (Collision) với LES Smagorinsky + Van Driest
#pragma omp parallel for
        for (int n = 0; n < total; ++n)
        {
            if (node_type[n] == NODE_SOLID)
            {
                // Trong ruột cánh: khóa vận tốc = 0
                grid.ux[n] = 0.0;
                grid.uy[n] = 0.0;
                continue;
            }

            double r = grid.rho[n];
            double u = grid.ux[n];
            double v = grid.uy[n];

            double f_node[Q];
            double feq_node[Q];
            for (int i = 0; i < Q; ++i)
            {
                f_node[i] = grid.f[i][n];
                feq_node[i] = Grid2D::compute_feq(i, r, u, v);
            }

            // Tính omega hiệu dụng có tắt vách Van Driest
            double omega_eff = les.compute_omega(r, f_node, feq_node, wall_dist[n]);

            // Thực hiện va chạm BGK
            for (int i = 0; i < Q; ++i)
            {
                grid.f[i][n] = f_node[i] - omega_eff * (f_node[i] - feq_node[i]);
            }
        }

// 3. Phản xạ Bounce-Back tại vách cánh cho các liên kết biên
#pragma omp parallel for
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
                            // Bounce-back: hạt phóng vào vách bị dội ngược lại hướng đối đỉnh
                            grid.fn[OPP[d]][n_idx] = grid.f[d][n_idx];
                        }
                    }
                }
            }
        }

// 4. Truyền hạt không gian (Streaming)
#pragma omp parallel for
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int n_idx = j * nx + i;
                if (node_type[n_idx] == NODE_SOLID)
                    continue;

                for (int d = 0; d < Q; ++d)
                {
                    int src_i = i - EX[d];
                    int src_j = j - EY[d];

                    if (src_i >= 0 && src_i < nx && src_j >= 0 && src_j < ny)
                    {
                        int src_idx = src_j * nx + src_i;
                        if (node_type[src_idx] != NODE_SOLID)
                        {
                            grid.fn[d][n_idx] = grid.f[d][src_idx];
                        }
                    }
                }
            }
        }

// 5. Điều kiện biên hầm gió số
#pragma omp parallel for
        for (int i = 0; i < nx; ++i)
        {
            // Biên đáy: Free-Slip (phản xạ gương trục y)
            int b_idx = 0 * nx + i;
            grid.fn[2][b_idx] = grid.fn[4][b_idx];
            grid.fn[5][b_idx] = grid.fn[8][b_idx];
            grid.fn[6][b_idx] = grid.fn[7][b_idx];

            // Biên đỉnh: Free-Slip
            int t_idx = (ny - 1) * nx + i;
            grid.fn[4][t_idx] = grid.fn[2][t_idx];
            grid.fn[7][t_idx] = grid.fn[6][t_idx];
            grid.fn[8][t_idx] = grid.fn[5][t_idx];
        }

#pragma omp parallel for
        for (int j = 0; j < ny; ++j)
        {
            // Biên vào (x = 0): Inflow cố định vận tốc u_inf
            int in_idx = j * nx + 0;
            for (int d = 0; d < Q; ++d)
            {
                grid.fn[d][in_idx] = Grid2D::compute_feq(d, 1.0, u_inf, 0.0);
            }

            // Biên ra (x = nx - 1): Convective Outflow trôi tự do
            int out_idx = j * nx + (nx - 1);
            int prev_idx = j * nx + (nx - 2);
            for (int d = 0; d < Q; ++d)
            {
                grid.fn[d][out_idx] = grid.fn[d][prev_idx];
            }
        }

        // 6. Cập nhật con trỏ fn -> f
        for (int d = 0; d < Q; ++d)
        {
            grid.f[d].swap(grid.fn[d]);
        }
    }

    void LBMSolver::run_steps(int num_steps)
    {
        for (int s = 0; s < num_steps; ++s)
        {
            step();
        }
    }

} // namespace lbm