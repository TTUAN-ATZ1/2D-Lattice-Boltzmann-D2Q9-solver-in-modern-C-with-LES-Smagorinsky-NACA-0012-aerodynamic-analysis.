#pragma once
#include <vector>
#include "D2Q9.hpp"

namespace lbm
{

    class LESSmagorinsky
    {
    public:
        double Cs;    // Hằng số Smagorinsky (chuẩn 0.16)
        double tau_0; // Thời gian hồi phục phân tử cơ sở
        double nu_0;  // Độ nhớt phân tử (lattice units)
        double u_tau; // Vận tốc ma sát vách u_tau

        LESSmagorinsky(double cs_val, double tau0, double nu0, double u_inf, double Re);

        // Tính tần số va chạm hiệu dụng omega_eff tại 1 nút lưới
        double compute_omega(double rho_val,
                             const double *f_node,
                             const double *feq_node,
                             double wall_distance) const;
    };

} // namespace lbm