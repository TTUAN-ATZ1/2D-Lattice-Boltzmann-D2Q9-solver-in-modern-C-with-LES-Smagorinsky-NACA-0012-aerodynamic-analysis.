#pragma once
#include <cmath>
#include "D2Q9.hpp"

namespace lbm
{

    class LESSmagorinsky
    {
    public:
        double Cs;        // Hằng số Smagorinsky (chuẩn 0.16)
        double tau_0;     // Thời gian hồi phục phân tử cơ sở
        double inv_tau_0; // 1.0 / tau_0 tính sẵn
        double nu_0;      // Độ nhớt phân tử (lattice units)
        double u_tau;     // Vận tốc ma sát vách u_tau
        double chord_lb;  // Chiều dài dây cung airfoil

        LESSmagorinsky(double cs_val, double tau0, double nu0, double u_inf, double Re, double chord = 150.0);

        // Tính sẵn hệ số tĩnh 18 * (Cs * f_vd)^2 tại mỗi ô lưới 1 lần duy nhất khi đặt hình học cánh
        double precompute_sgs_coeff(double wall_distance, double x_coord, double x_te) const;

        // Tính nhanh tần số va chạm hiệu dụng omega_eff trong vòng lặp thời gian (inline, không gọi exp)
        inline double compute_omega_fast(double rho_val,
                                         const double *f_node,
                                         const double *feq_node,
                                         double sgs_coeff) const
        {
            if (sgs_coeff == 0.0)
                return inv_tau_0;

            double pi_xx = 0.0, pi_yy = 0.0, pi_xy = 0.0;
            for (int i = 1; i < Q; ++i)
            {
                const double f_neq = f_node[i] - feq_node[i];
                pi_xx += EX[i] * EX[i] * f_neq;
                pi_yy += EY[i] * EY[i] * f_neq;
                pi_xy += EX[i] * EY[i] * f_neq;
            }

            const double Q_val = pi_xx * pi_xx + pi_yy * pi_yy + 2.0 * (pi_xy * pi_xy);
            const double C = (sgs_coeff / rho_val) * std::sqrt(2.0 * Q_val);
            const double tau_eff = 0.5 * (tau_0 + std::sqrt(tau_0 * tau_0 + C));
            return 1.0 / tau_eff;
        }

        double compute_omega(double rho_val,
                             const double *f_node,
                             const double *feq_node,
                             double wall_distance,
                             double x_coord,
                             double x_te) const;
    };

} // namespace lbm