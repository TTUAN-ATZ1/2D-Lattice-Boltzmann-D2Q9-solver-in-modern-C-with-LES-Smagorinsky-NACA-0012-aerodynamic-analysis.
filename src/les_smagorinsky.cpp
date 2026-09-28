#include "les_smagorinsky.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace lbm
{

    LESSmagorinsky::LESSmagorinsky(double cs_val, double tau0, double nu0, double u_inf, double Re, double chord)
        : Cs(cs_val), tau_0(tau0), inv_tau_0(1.0 / tau0), nu_0(nu0), chord_lb(chord)
    {
        if (!std::isfinite(Cs)||Cs<0||!std::isfinite(tau_0)||tau_0<=0.5
            ||!std::isfinite(nu_0)||nu_0<=0||!std::isfinite(Re)||Re<=0
            ||!std::isfinite(u_inf)||u_inf<=0||!std::isfinite(chord)||chord<=0)
            throw std::invalid_argument("Invalid LBM viscosity/Smagorinsky parameters");
        double cdf = 1.328 / std::sqrt(Re);
        u_tau = u_inf * std::sqrt(0.5 * cdf);
    }

    double LESSmagorinsky::precompute_sgs_coeff(double wall_distance, double x_coord, double x_te) const
    {
        if (Cs == 0.0)
            return 0.0;
        bool in_airfoil_region = (wall_distance < 0.15 * chord_lb) && (x_coord < x_te + 0.05 * chord_lb);
        double f_vd = 1.0;
        if (in_airfoil_region)
        {
            double y_plus = (wall_distance * u_tau) / (nu_0 + 1e-12);
            f_vd = 1.0 - std::exp(-y_plus / 26.0);
            f_vd = std::max(0.0, std::min(1.0, f_vd));
        }
        double ls = Cs * f_vd;
        return 18.0 * ls * ls;
    }

    double LESSmagorinsky::compute_omega(double rho_val,
                                         const double *f_node,
                                         const double *feq_node,
                                         double wall_distance,
                                         double x_coord,
                                         double x_te) const
    {
        if (!std::isfinite(rho_val)||rho_val<=0)
            return std::numeric_limits<double>::quiet_NaN();
        return compute_omega_fast(rho_val, f_node, feq_node, precompute_sgs_coeff(wall_distance, x_coord, x_te));
    }

} // namespace lbm
