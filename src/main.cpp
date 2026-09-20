#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <sstream>

#include "solver.hpp"
#include "geometry_naca.hpp"
#include "momentum_exchange.hpp"
#include "integrator.hpp"
#include "vtk_exporter.hpp"

namespace fs = std::filesystem;

int main()
{
    // 1. Thiết lập thông số mô phỏng độ phân giải cao
    const int Nx = 640;
    const int Ny = 360;
    const double chord_lb = 150.0; // Dây cung 150 nodes (gấp 3 lần bản cũ)
    const double u_inf_lb = 0.04;  // Vận tốc mạng (Ma = 0.069 << 0.3)
    const double Re = 100000.0;
    const double Cs = 0.16; // Hằng số Smagorinsky 0.16

    const double x_pivot = Nx * 0.25;
    const double y_pivot = Ny * 0.50;

    // Tạo các thư mục lưu kết quả
    fs::create_directories("results_vtk");
    fs::create_directories("results_csv");

    // Khởi tạo bộ giải LBM
    lbm::LBMSolver solver(Nx, Ny, u_inf_lb, Re, chord_lb, Cs);
    lbm::MomentumExchange mea(chord_lb, u_inf_lb, Re, x_pivot, y_pivot);
    lbm::Slide32Integrator slide32(chord_lb, u_inf_lb, Re);
    lbm::VTKExporter vtk(u_inf_lb);

    // Mở file lưu bảng đặc tuyến Polar tổng hợp
    std::ofstream polar_out("results_csv/polar_comparison_Re100k.csv");
    polar_out << "Alpha,Cl_MEA,Cd_MEA,x_cp_MEA,Cl_Slide32,Cd_Slide32,x_cp_Slide32\n";

    std::cout << " Alpha | Cl(MEA) | Cd(MEA) | x_cp(MEA) | Cl(Sl32) | Cd(Sl32) |   TTime  \n";
    std::cout << "--------------------------------------------------------------------------------\n"
              << std::flush;

    auto start_all = std::chrono::high_resolution_clock::now();

    // Quét góc tấn từ 0 đến 20 độ
    for (int a = 0; a <= 30; ++a)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        double alpha_deg = static_cast<double>(a);

        // 1. Tạo hình học cánh và phân loại lưới
        lbm::NACAGeometry geom(chord_lb, alpha_deg, x_pivot, y_pivot, 120);
        solver.set_geometry(geom);

        // 2. Xây dựng liên kết cho MEA và Stencils cho Slide 3.2
        mea.build_boundary_links(Nx, Ny, solver.node_type, alpha_deg);
        slide32.build_stencils(Nx, Ny, geom.surface_points);

        // 3. Số bước lặp thời gian thích nghi (Adaptive Time-Averaging):
        // - Trước thất tốc (a <= 10): Dòng dừng bám dính, hội tụ nhanh (2000 bước, avg 500)
        // - Sau thất tốc (a >= 11): Dòng không dừng nhả xoáy (Strouhal T ~ 5400 steps).
        //   Cần 10000 bước: 4000 bước xả nhiễu quá độ (transient), 6000 bước lấy trung bình trọn chu kỳ xoáy.
        int steps;
        int n_avg;

        if (a <= 10)
        {
            steps = (a == 0) ? 3000 : 2000;
            n_avg = 500;
        }
        else
        {
            steps = 10000; // Đủ để vượt qua giai đoạn quá độ và bao phủ chu kỳ xoáy
            n_avg = 6000;  // Lấy trung bình 6000 bước cuối (> 1 chu kỳ T ~ 5400 bước)
        }

        double sum_cl_mea = 0.0, sum_cd_mea = 0.0, sum_xcp_mea = 0.0;
        double sum_cl_sl = 0.0, sum_cd_sl = 0.0, sum_xcp_sl = 0.0;
        int count_avg = 0;

        for (int s = 0; s < steps; ++s)
        {
            solver.step();

            // Lấy trung bình các bước cuối để triệt tiêu dao động pha của chuỗi xoáy
            if (s >= (steps - n_avg))
            {
                auto f_mea = mea.compute_force(solver.grid);
                auto f_sl = slide32.integrate_pressure(solver.grid, geom);

                sum_cl_mea += f_mea.Cl;
                sum_cd_mea += f_mea.Cd;
                sum_xcp_mea += f_mea.x_cp;

                sum_cl_sl += f_sl.Cl;
                sum_cd_sl += f_sl.Cd;
                sum_xcp_sl += f_sl.x_cp;
                count_avg++;
            }
        }

        double cl_mea = sum_cl_mea / count_avg;
        double cd_mea = sum_cd_mea / count_avg;
        double xcp_mea = sum_xcp_mea / count_avg;

        double cl_sl = sum_cl_sl / count_avg;
        double cd_sl = sum_cd_sl / count_avg;
        double xcp_sl = sum_xcp_sl / count_avg;

        auto t1 = std::chrono::high_resolution_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();

        // In kết quả từng góc ra màn hình
        std::cout << std::fixed << std::setprecision(1) << std::setw(5) << alpha_deg << " | "
                  << std::setprecision(3) << std::setw(7) << cl_mea << " | "
                  << std::setw(7) << cd_mea << " | "
                  << std::setw(8) << xcp_mea << "c | "
                  << std::setw(8) << cl_sl << " | "
                  << std::setw(8) << cd_sl << " | "
                  << std::setprecision(1) << std::setw(8) << sec << "s\n"
                  << std::flush;

        // Ghi vào file CSV đối chứng
        polar_out << alpha_deg << ","
                  << cl_mea << "," << cd_mea << "," << xcp_mea << ","
                  << cl_sl << "," << cd_sl << "," << xcp_sl << "\n";
        polar_out.flush();

        // Xuất file ParaView cho toàn bộ chuỗi 00 -> 30 để tạo Animation
        std::stringstream ss_idx;
        ss_idx << std::setw(2) << std::setfill('0') << a;
        std::string idx_str = ss_idx.str();

        // 1. Trường dòng 2D (Pressure, Velocity, Vorticity, Airfoil Mask)
        std::string vtk_flow = "results_vtk/naca0012_flow_" + idx_str + ".vtk";
        vtk.export_field_2d(vtk_flow, solver.grid, solver.node_type);

        std::string csv_file = "results_csv/slide32_points_a" + std::to_string(a) + ".csv";
        slide32.export_csv(csv_file, solver.grid, geom, alpha_deg);
    }

    polar_out.close();
    auto end_all = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_all - start_all).count();

    return 0;
}