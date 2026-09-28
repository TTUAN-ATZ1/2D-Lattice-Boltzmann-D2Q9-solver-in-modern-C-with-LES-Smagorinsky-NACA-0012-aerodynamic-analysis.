#pragma once
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace lbm {
struct PhysicalUnits {
    double u_lb, rho, velocity, pressure, nu, chord, dx, dt;
    PhysicalUnits(double chord_lb, double lattice_velocity, double re,
                  double real_velocity=10.0, double real_nu=1.5e-5,
                  double real_rho=1.225, double real_pressure=101325.0)
        : u_lb(lattice_velocity), rho(real_rho), velocity(real_velocity),
          pressure(real_pressure), nu(real_nu) {
        for (double v:{chord_lb,u_lb,re,velocity,nu,rho,pressure})
            if (!std::isfinite(v)||v<=0) throw std::invalid_argument("Invalid physical-unit parameter");
        chord=re*nu/velocity;
        dx=chord/chord_lb;
        dt=dx*u_lb/velocity;
    }
    double q() const { return 0.5*rho*velocity*velocity; }
};
}
