#include "les_smagorinsky.hpp"
#include <cmath>
#include <algorithm>

namespace lbm
{

    LESSmagorinsky::LESSmagorinsky(double cs_val, double tau0, double nu0, double u_inf, double Re)
        : Cs(cs_val), tau_0(tau0), nu_0(nu0)
    {

        // Ước tính hệ số ma sát Blasius và vận tốc ma sát vách u_tau
        double cdf = 1.328 / std::sqrt(Re);
        u_tau = u_inf * std::sqrt(0.5 * cdf);
    }

    double LESSmagorinsky::compute_omega(double rho_val,
                                         const double *f_node,
                                         const double *feq_node,
                                         double wall_distance) const
    {
        // 1. Tính các thành phần của tensor ứng suất không cân bằng Pi^(1)
        double pi_xx = 0.0;
        double pi_yy = 0.0;
        double pi_xy = 0.0;

        for (int i = 1; i < Q; ++i)
        {
            double f_neq = f_node[i] - feq_node[i];
            pi_xx += EX[i] * EX[i] * f_neq;
            pi_yy += EY[i] * EY[i] * f_neq;
            pi_xy += EX[i] * EY[i] * f_neq;
        }

        // Bất biến bậc hai Q = Pi_ij * Pi_ij
        double Q_val = pi_xx * pi_xx + pi_yy * pi_yy + 2.0 * (pi_xy * pi_xy);
        double sqrt_Q = std::sqrt(2.0 * Q_val);

        // 2. Hệ số tắt vách Van Driest Damping
        double y_plus = (wall_distance * u_tau) / (nu_0 + 1e-12);
        double f_vd = 1.0 - std::exp(-y_plus / 26.0);
        f_vd = std::max(0.0, std::min(1.0, f_vd)); // Giới hạn [0, 1]

        // Chiều dài lọc Smagorinsky có hiệu chỉnh tắt vách
        double ls = Cs * f_vd;
        double ls_sq = ls * ls;

        // 3. Giải phương trình bậc 2 tìm tau_eff
        // 4K = (18 * ls^2 / rho) * sqrt(2 * Q) = (18 * ls^2 / rho) * sqrt_Q
        constexpr double COEFF = 18.0;
        double C = (COEFF * ls_sq / (rho_val + 1e-12)) * sqrt_Q;

        double discr = tau_0 * tau_0 + C;
        double tau_eff = 0.5 * (tau_0 + std::sqrt(discr));

        // Đảm bảo an toàn số trị: tau_eff >= tau_0 và omega <= 1.99
        tau_eff = std::max(tau_eff, tau_0);
        double omega_eff = 1.0 / tau_eff;
        return std::min(omega_eff, 1.99);
    }

}