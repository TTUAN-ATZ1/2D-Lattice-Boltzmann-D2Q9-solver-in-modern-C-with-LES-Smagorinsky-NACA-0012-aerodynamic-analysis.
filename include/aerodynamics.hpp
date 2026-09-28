#pragma once
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lbm {
constexpr double PI = 3.14159265358979323846;
inline double unavailable() { return std::numeric_limits<double>::quiet_NaN(); }

struct AeroForce {
    double Fx = 0, Fy = 0, Mz = 0; // lattice force per span; moment about quarter chord, +CCW
    double Cl = 0, Cd = 0;
    double Cn = 0, Cm = 0;
    double x_cp = unavailable(); // chord fraction from LE, not clamped
};

inline AeroForce aerodynamic_force(double fx, double fy, double mz, double alpha_rad,
                                    double chord, double q_chord) {
    if (!(chord > 0) || !(q_chord > 0) || !std::isfinite(chord) || !std::isfinite(q_chord))
        throw std::invalid_argument("Invalid force normalization");
    if (!std::isfinite(fx) || !std::isfinite(fy) || !std::isfinite(mz))
        throw std::runtime_error("Non-finite aerodynamic force or moment");
    AeroForce f;
    f.Fx=fx; f.Fy=fy; f.Mz=mz;
    f.Cl=fy/q_chord; f.Cd=fx/q_chord;
    f.Cn=(fx*std::sin(alpha_rad)+fy*std::cos(alpha_rad))/q_chord;
    f.Cm=mz/(q_chord*chord);
    if (std::abs(f.Cn)>1e-8) f.x_cp=0.25+f.Cm/f.Cn;
    return f;
}

struct ACSample { double alpha, cn, cm; };
inline double aerodynamic_center(const std::vector<ACSample>& samples) {
    if (samples.size()<3) return unavailable();
    double a=0,n=0,m=0;
    for (const auto& s:samples) {
        if (!std::isfinite(s.alpha)||!std::isfinite(s.cn)||!std::isfinite(s.cm))
            return unavailable();
        a+=s.alpha; n+=s.cn; m+=s.cm;
    }
    a/=samples.size(); n/=samples.size(); m/=samples.size();
    double va=0,an=0,am=0;
    for (const auto& s:samples) {
        const double da=s.alpha-a;
        va+=da*da; an+=da*(s.cn-n); am+=da*(s.cm-m);
    }
    if (va<=1e-20 || std::abs(an/va)<1e-8) return unavailable();
    return 0.25+am/an;
}
}
