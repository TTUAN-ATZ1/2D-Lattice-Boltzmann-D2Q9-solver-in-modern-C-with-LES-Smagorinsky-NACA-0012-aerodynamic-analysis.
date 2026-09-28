#pragma once
#include "D2Q9.hpp"

namespace lbm {
enum class CollisionModel { BGK, MRT };
constexpr double MRT_ENERGY_RATE=1.1;
constexpr double MRT_EPSILON_RATE=1.0;
// 1.2 admitted growing finite-wavelength modes at nu=6e-5, U=0.04.
// 1.98 passes the documented moving-flow checks; it is NOT a viscosity cap.
constexpr double MRT_ENERGY_FLUX_RATE=1.98;

// Orthogonal D2Q9 raw moments (Lallemand & Luo, 2000).
// Conserved density/momentum are untouched. Only the two shear modes use
// omega=1/(0.5+3*nu); bulk/kinetic relaxation does not clip molecular viscosity.
inline void collide_mrt(double f[Q], const double eq[Q], double omega) {
    double g[Q];
    for (int d=0;d<Q;++d) g[d]=f[d]-eq[d];
    const double axes=g[1]+g[2]+g[3]+g[4];
    const double diag=g[5]+g[6]+g[7]+g[8];
    const double e=MRT_ENERGY_RATE*(-4*g[0]-axes+2*diag)/36;
    const double eps=MRT_EPSILON_RATE*(4*g[0]-2*axes+diag)/36;
    const double qx=MRT_ENERGY_FLUX_RATE*(-2*g[1]+2*g[3]+g[5]-g[6]-g[7]+g[8])/12;
    const double qy=MRT_ENERGY_FLUX_RATE*(-2*g[2]+2*g[4]+g[5]+g[6]-g[7]-g[8])/12;
    const double pxx=omega*(g[1]-g[2]+g[3]-g[4])/4;
    const double pxy=omega*(g[5]-g[6]+g[7]-g[8])/4;
    f[0]-=-4*e+4*eps;
    f[1]-=-e-2*eps-2*qx+pxx;
    f[2]-=-e-2*eps-2*qy-pxx;
    f[3]-=-e-2*eps+2*qx+pxx;
    f[4]-=-e-2*eps+2*qy-pxx;
    f[5]-=2*e+eps+qx+qy+pxy;
    f[6]-=2*e+eps-qx+qy-pxy;
    f[7]-=2*e+eps-qx-qy+pxy;
    f[8]-=2*e+eps+qx-qy-pxy;
}
}
