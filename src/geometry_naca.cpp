#define _USE_MATH_DEFINES
#include <cmath>
#include "geometry_naca.hpp"
#include <algorithm>

namespace lbm
{

    NACAGeometry::NACAGeometry(double chord_len, double alpha, double xp, double yp, int num_points)
        : chord(chord_len), alpha_deg(alpha), x_pivot(xp), y_pivot(yp), Nb(num_points)
    {
        build_profile();
    }

    void NACAGeometry::build_profile()
    {
        surface_points.clear();
        raw_points.clear();
        ds.clear();
        normal.clear();

        int n_half = (Nb / 2) + 1; // 61 diem
        std::vector<double> xc(n_half);
        std::vector<double> yt(n_half);

        // phan bo diem Cosine
        for (int i = 0; i < n_half; ++i)
        {
            double beta = M_PI * i / (n_half - 1);
            xc[i] = 0.5 * chord * (1.0 - std::cos(beta));
            double x_rel = std::max(xc[i] / chord, 1e-6);

            // Khep kin mep sau
            yt[i] = 5.0 * 0.12 * chord * (0.2969 * std::sqrt(x_rel) - 0.1260 * x_rel - 0.3516 * (x_rel * x_rel) + 0.2843 * (x_rel * x_rel * x_rel) - 0.1036 * (x_rel * x_rel * x_rel * x_rel));
        }

        // di theo chieu CGW: mep duoi -> lung canh ->mui -> bung canh
        // Lung canh (xc dao nguoc: c -> 0)
        for (int i = n_half - 1; i >= 0; --i)
        {
            raw_points.push_back({xc[i], yt[i]});
        }
        // Bung canh (xc chay tu sau mui ve truoc duoi: 1 -> n_half-2)
        for (int i = 1; i < n_half - 1; ++i)
        {
            raw_points.push_back({xc[i], -yt[i]});
        }
        // Khep kin ve diem dau
        if (raw_points.size() < (size_t)Nb)
        {
            raw_points.push_back(raw_points[0]);
        }

        //  Ma tran xoay theo goc tan xung quanh 0.25c
        double ar = alpha_deg * M_PI / 180.0;
        double cos_a = std::cos(ar);
        double sin_a = std::sin(ar);
        double c_quarter = 0.25 * chord;

        for (size_t k = 0; k < raw_points.size(); ++k)
        {
            double dx = raw_points[k].x - c_quarter;
            double dy = raw_points[k].y;
            double x_rot = x_pivot + dx * cos_a + dy * sin_a;
            double y_rot = y_pivot - dx * sin_a + dy * cos_a;
            surface_points.push_back({x_rot, y_rot});
        }

        // vector phap tuyen ngoai chuan tac n = (dy/ds, -dx/ds)
        size_t N = surface_points.size();
        ds.resize(N);
        normal.resize(N);

        for (size_t k = 0; k < N; ++k)
        {
            size_t k_next = (k + 1) % N;
            double dx_seg = surface_points[k_next].x - surface_points[k].x;
            double dy_seg = surface_points[k_next].y - surface_points[k].y;
            double len = std::sqrt(dx_seg * dx_seg + dy_seg * dy_seg);
            ds[k] = len;
            normal[k].x = dy_seg / (len + 1e-12);
            normal[k].y = -dx_seg / (len + 1e-12);
        }
    }

    double NACAGeometry::distance_to_segment(double px, double py, double ax, double ay, double bx, double by)
    {
        double abx = bx - ax;
        double aby = by - ay;
        double apx = px - ax;
        double apy = py - ay;
        double ab_len_sq = abx * abx + aby * aby;
        if (ab_len_sq < 1e-12)
            return std::sqrt(apx * apx + apy * apy);

        double t = (apx * abx + apy * aby) / ab_len_sq;
        t = std::max(0.0, std::min(1.0, t));
        double proj_x = ax + t * abx;
        double proj_y = ay + t * aby;
        double dx = px - proj_x;
        double dy = py - proj_y;
        return std::sqrt(dx * dx + dy * dy);
    }

    void NACAGeometry::classify_grid(int nx, int ny, std::vector<NodeType> &node_type, std::vector<double> &wall_dist) const
    {
        int total = nx * ny;
        node_type.assign(total, NODE_FLUID);
        wall_dist.assign(total, 1e6);

        size_t N = surface_points.size();

//  Thuật toán Ray Casting xác định SOLID vs FLUID & tính khoảng cách vách
#pragma omp parallel for
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int n_idx = j * nx + i;
                double px = i;
                double py = j;

                // Kiểm tra trong vùng bounding box của cánh
                double d_min = 1e6;
                int cross_count = 0;

                for (size_t k = 0; k < N; ++k)
                {
                    size_t k_next = (k + 1) % N;
                    double ax = surface_points[k].x;
                    double ay = surface_points[k].y;
                    double bx = surface_points[k_next].x;
                    double by = surface_points[k_next].y;

                    // Tính khoảng cách tới đoạn biên k
                    double d_seg = distance_to_segment(px, py, ax, ay, bx, by);
                    if (d_seg < d_min)
                        d_min = d_seg;

                    // Ray casting theo tia nằm ngang sang phải
                    if ((ay > py) != (by > py))
                    {
                        double x_intersect = (bx - ax) * (py - ay) / (by - ay + 1e-12) + ax;
                        if (px < x_intersect)
                        {
                            cross_count++;
                        }
                    }
                }

                wall_dist[n_idx] = d_min;
                if (cross_count % 2 == 1)
                {
                    node_type[n_idx] = NODE_SOLID;
                }
            }
        }

// Xác định các nút BOUNDARY (FLUID tiếp giáp với SOLID) cho thuật toán MEA
#pragma omp parallel for
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                int n_idx = j * nx + i;
                if (node_type[n_idx] == NODE_FLUID)
                {
                    bool near_solid = false;
                    for (int d = 1; d < Q; ++d)
                    {
                        int ni = i + EX[d];
                        int nj = j + EY[d];
                        int neighbor_idx = nj * nx + ni;
                        if (node_type[neighbor_idx] == NODE_SOLID)
                        {
                            near_solid = true;
                            break;
                        }
                    }
                    if (near_solid)
                    {
                        node_type[n_idx] = NODE_BOUNDARY;
                    }
                }
            }
        }
    }

}