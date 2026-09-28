#include "solver.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <string>

namespace lbm
{

    LBMSolver::LBMSolver(int nx_, int ny_, double u_inf_, double Re, double chord_lb_, double Cs, OuterBoundary boundary, CollisionModel collision)
        : nx(nx_), ny(ny_), u_inf(u_inf_), chord_lb(chord_lb_), x_te_lab(0.0), outer_boundary(boundary), collision_model(collision),
          grid(nx_, ny_),
          les(Cs, 3.0 * (u_inf_ * chord_lb_ / Re) + 0.5, (u_inf_ * chord_lb_ / Re), u_inf_, Re, chord_lb_)
    {
        if (collision_model==CollisionModel::MRT && Cs!=0)
            throw std::invalid_argument("MRT currently requires cs=0; the BGK SGS closure is not reused for MRT");
        last_mea_force = {};
        grid.initialize(u_inf);
    }

    void LBMSolver::set_geometry(const NACAGeometry &geom)
    {
        double max_x = 0.0;
        if (!geom.surface_points.empty())
        {
            max_x = geom.surface_points[0].x;
            for (const auto &pt : geom.surface_points)
            {
                max_x = std::max(max_x, pt.x);
            }
        }
        x_te_lab = max_x;

        geom.classify_grid(nx, ny, node_type, wall_dist);
        sgs_coeff.resize(grid.size);
#pragma omp parallel for
        for (int n = 0; n < grid.size; ++n)
        {
            sgs_coeff[n] = (node_type[n] == NODE_SOLID)
                               ? 0.0
                               : les.precompute_sgs_coeff(wall_dist[n], static_cast<double>(n % nx), x_te_lab);
        }

        MomentumExchange mea_builder(chord_lb, u_inf, geom.x_pivot, geom.y_pivot);
        mea_builder.build_boundary_links(nx, ny, node_type, geom);
        ibb_links = std::move(mea_builder.links);
        refresh_macroscopic();
    }

    void LBMSolver::refresh_macroscopic()
    {
        if (node_type.size()!=static_cast<size_t>(grid.size))
            throw std::logic_error("Set geometry before running the solver");
        int bad=grid.size;
#pragma omp parallel for reduction(min:bad)
        for (int n=0;n<grid.size;++n) {
            if (node_type[n]==NODE_SOLID) {
                grid.ux[n]=0.0;
                grid.uy[n]=0.0;
                continue;
            }
            const double f0=grid.f[0][n], f1=grid.f[1][n], f2=grid.f[2][n];
            const double f3=grid.f[3][n], f4=grid.f[4][n], f5=grid.f[5][n];
            const double f6=grid.f[6][n], f7=grid.f[7][n], f8=grid.f[8][n];
            const double r = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            grid.rho[n]=r;
            if (!std::isfinite(r) || r<=1e-12) {
                grid.ux[n]=std::numeric_limits<double>::quiet_NaN();
                grid.uy[n]=std::numeric_limits<double>::quiet_NaN();
                bad=std::min(bad,n);
                continue;
            }
            const double u = ((f1 - f3) + (f5 - f6 - f7 + f8)) / r;
            const double v = ((f2 - f4) + (f5 + f6 - f7 - f8)) / r;
            grid.ux[n]=u;
            grid.uy[n]=v;
            if (!std::isfinite(u) || !std::isfinite(v)) bad=std::min(bad,n);
        }
        if (bad<grid.size) throw std::runtime_error("Invalid fluid state at x="+
            std::to_string(bad%nx)+", y="+std::to_string(bad/nx));
    }

    FlowHealth LBMSolver::health() const
    {
        if (node_type.size()!=static_cast<size_t>(grid.size))
            throw std::logic_error("Set geometry before diagnostics");
        double mass=0,rmin=std::numeric_limits<double>::infinity(),rmax=0,speed=0;
        double fin=0,fout=0,minf=std::numeric_limits<double>::infinity();
#pragma omp parallel for reduction(+:mass,fin,fout) reduction(min:rmin,minf) reduction(max:rmax,speed)
        for (int n=0;n<grid.size;++n) {
            if (node_type[n]==NODE_SOLID) continue;
            mass+=grid.rho[n]; rmin=std::min(rmin,grid.rho[n]); rmax=std::max(rmax,grid.rho[n]);
            speed=std::max(speed,std::hypot(grid.ux[n],grid.uy[n]));
            const double w=(n/nx==0 || n/nx==ny-1)?0.5:1.0;
            if (n%nx==0) fin+=w*grid.rho[n]*grid.ux[n];
            if (n%nx==nx-1) fout+=w*grid.rho[n]*grid.ux[n];
            for (int d=0;d<Q;++d) minf=std::min(minf,grid.f[d][n]);
        }
        return {mass,rmin,rmax,speed/std::sqrt(CS2),fin,fout,minf};
    }

    void LBMSolver::step(const MomentumExchange *mea)
    {
        int total = nx * ny;
        if (wall_dist.size()!=static_cast<size_t>(total))
            throw std::logic_error("Set geometry before stepping");
        if (sgs_coeff.size()!=static_cast<size_t>(total))
        {
            sgs_coeff.resize(total);
#pragma omp parallel for
            for (int n = 0; n < total; ++n)
                sgs_coeff[n] = (node_type[n] == NODE_SOLID)
                                   ? 0.0
                                   : les.precompute_sgs_coeff(wall_dist[n], static_cast<double>(n % nx), x_te_lab);
        }

        int bad = total;
        int num_ibb = static_cast<int>(ibb_links.size());

#pragma omp parallel
        {
            // 1. Gộp tính biến vĩ mô (rho, ux, uy) và Va chạm (Collision LES/MRT) trong 1 lần đọc bộ nhớ
#pragma omp for reduction(min:bad)
            for (int n = 0; n < total; ++n)
            {
                if (node_type[n] == NODE_SOLID)
                    continue;

                double f_node[Q] = {
                    grid.f[0][n], grid.f[1][n], grid.f[2][n],
                    grid.f[3][n], grid.f[4][n], grid.f[5][n],
                    grid.f[6][n], grid.f[7][n], grid.f[8][n]
                };

                const double r = f_node[0] + f_node[1] + f_node[2] + f_node[3] + f_node[4]
                               + f_node[5] + f_node[6] + f_node[7] + f_node[8];
                if (!std::isfinite(r) || r <= 1e-12)
                {
                    bad = std::min(bad, n);
                    continue;
                }
                const double inv_r = 1.0 / r;
                const double u = ((f_node[1] - f_node[3]) + (f_node[5] - f_node[6] - f_node[7] + f_node[8])) * inv_r;
                const double v = ((f_node[2] - f_node[4]) + (f_node[5] + f_node[6] - f_node[7] - f_node[8])) * inv_r;
                grid.rho[n] = r;
                grid.ux[n] = u;
                grid.uy[n] = v;

                double feq_node[Q];
                Grid2D::compute_feq_all(r, u, v, feq_node);

                const double omega_eff = les.compute_omega_fast(r, f_node, feq_node, sgs_coeff[n]);

                if (collision_model == CollisionModel::MRT)
                    collide_mrt(f_node, feq_node, omega_eff);
                else
                    for (int i = 0; i < Q; ++i)
                        f_node[i] -= omega_eff * (f_node[i] - feq_node[i]);

                for (int i = 0; i < Q; ++i)
                {
                    grid.f[i][n] = f_node[i];
                    if (!std::isfinite(f_node[i]))
                        bad = std::min(bad, n);
                }
            }

            // 2. Phản xạ Interpolated Bounce-Back (Bouzidi bậc 1) tại vách cong của cánh
#pragma omp for
            for (int k = 0; k < num_ibb; ++k)
            {
                const auto &link = ibb_links[k];
                int n_idx = grid.idx(link.x, link.y);
                grid.fn[link.opp_dir][n_idx] = compute_ibb_reflected(link, grid);
            }

            // 3. Truyền hạt không gian (Streaming) - tối ưu hóa nhánh cho nút FLUID bên trong miền
#pragma omp for
            for (int j = 0; j < ny; ++j)
            {
                const bool interior_j = (outer_boundary == OuterBoundary::Channel) && (j > 0 && j < ny - 1);
                for (int i = 0; i < nx; ++i)
                {
                    int n_idx = j * nx + i;
                    const NodeType nt = node_type[n_idx];
                    if (nt == NODE_SOLID)
                        continue;

                    if (interior_j && nt == NODE_FLUID && i > 0 && i < nx - 1)
                    {
                        grid.fn[0][n_idx] = grid.f[0][n_idx];
                        grid.fn[1][n_idx] = grid.f[1][n_idx - 1];
                        grid.fn[2][n_idx] = grid.f[2][n_idx - nx];
                        grid.fn[3][n_idx] = grid.f[3][n_idx + 1];
                        grid.fn[4][n_idx] = grid.f[4][n_idx + nx];
                        grid.fn[5][n_idx] = grid.f[5][n_idx - nx - 1];
                        grid.fn[6][n_idx] = grid.f[6][n_idx - nx + 1];
                        grid.fn[7][n_idx] = grid.f[7][n_idx + nx + 1];
                        grid.fn[8][n_idx] = grid.f[8][n_idx + nx - 1];
                        continue;
                    }

                    for (int d = 0; d < Q; ++d)
                    {
                        int src_i = i - EX[d];
                        int src_j = j - EY[d];
                        if (outer_boundary == OuterBoundary::Periodic)
                        {
                            src_i = (src_i + nx) % nx;
                            src_j = (src_j + ny) % ny;
                        }

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

            // 4. Điều kiện biên hầm gió số
            if (outer_boundary == OuterBoundary::Channel)
            {
#pragma omp for
                for (int i = 0; i < nx; ++i)
                {
                    int b_idx = i;
                    grid.fn[2][b_idx] = grid.fn[4][b_idx];
                    grid.fn[5][b_idx] = grid.fn[8][b_idx];
                    grid.fn[6][b_idx] = grid.fn[7][b_idx];

                    int t_idx = (ny - 1) * nx + i;
                    grid.fn[4][t_idx] = grid.fn[2][t_idx];
                    grid.fn[7][t_idx] = grid.fn[6][t_idx];
                    grid.fn[8][t_idx] = grid.fn[5][t_idx];
                }

#pragma omp for
                for (int j = 0; j < ny; ++j)
                {
                    for (int side = 0; side < 2; ++side)
                    {
                        const int boundary = j * nx + (side == 0 ? 0 : nx - 1);
                        const int interior = boundary + (side == 0 ? 1 : -1);
                        double r = 0, mx = 0, my = 0;
                        for (int d = 0; d < Q; ++d)
                        {
                            const double fi = grid.fn[d][interior];
                            r += fi;
                            mx += fi * EX[d];
                            my += fi * EY[d];
                        }
                        const double u = mx / r, v = my / r;
                        const double rb = side == 0 ? r : 1.0;
                        const double ub = side == 0 ? u_inf : u;
                        const double vb = (side == 0 || j == 0 || j == ny - 1) ? 0.0 : v;
                        for (int d = 0; d < Q; ++d)
                            grid.fn[d][boundary] = Grid2D::compute_feq(d, rb, ub, vb)
                                                 + grid.fn[d][interior] - Grid2D::compute_feq(d, r, u, v);
                    }
                }
            }
        } // end parallel region

        if (bad < total)
            throw std::runtime_error("Invalid fluid state or collision at x=" +
                                     std::to_string(bad % nx) + ", y=" + std::to_string(bad / nx));

        // 5. TÍNH LỰC MEA (khi bước hiện tại cần lấy mẫu)
        if (mea)
        {
            last_mea_force = mea->compute_force(grid, /*is_post_collision=*/true);
        }

        // 6. Hoán đổi con trỏ fn <-> f
        for (int d = 0; d < Q; ++d)
        {
            grid.f[d].swap(grid.fn[d]);
        }
    }

    void LBMSolver::run_steps(int num_steps, const MomentumExchange *mea)
    {
        for (int s = 0; s < num_steps; ++s)
        {
            step(mea);
        }
        refresh_macroscopic();
    }

} // namespace lbm
