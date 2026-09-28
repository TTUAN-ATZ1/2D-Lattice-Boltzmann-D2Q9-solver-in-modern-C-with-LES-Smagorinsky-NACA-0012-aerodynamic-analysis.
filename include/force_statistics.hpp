#pragma once
#include "aerodynamics.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace lbm {
struct ForceSummary {
    AeroForce mean;
    double cl_rms=0, se_cl=0, se_cd=0, se_cn=0, se_cm=0;
    bool block_stable=false;
    size_t samples=0;
};

inline ForceSummary summarize_forces(const std::vector<AeroForce>& samples,double alpha,
        double chord,double q_chord,double abs_tol,double rel_tol,
        bool enough_block_time) {
    if (samples.size()<20) throw std::invalid_argument("Need at least 20 averaging samples");
    ForceSummary out;
    out.samples=samples.size();
    double fx=0,fy=0,mz=0;
    for (const auto& f:samples) { fx+=f.Fx; fy+=f.Fy; mz+=f.Mz; }
    const double count=static_cast<double>(samples.size());
    out.mean=aerodynamic_force(fx/count,fy/count,mz/count,alpha,chord,q_chord);
    std::array<std::array<double,4>,5> blocks{}; // Cl, Cd, Cn, Cm
    for (size_t b=0;b<5;++b) {
        size_t begin=b*samples.size()/5,end=(b+1)*samples.size()/5;
        for (size_t i=begin;i<end;++i) {
            blocks[b][0]+=samples[i].Cl; blocks[b][1]+=samples[i].Cd;
            blocks[b][2]+=samples[i].Cn; blocks[b][3]+=samples[i].Cm;
        }
        for (auto& v:blocks[b]) v/=end-begin;
    }
    std::array<double,4> errors{};
    out.block_stable=enough_block_time;
    for (size_t j=0;j<4;++j) {
        double mean=0,lo=blocks[0][j],hi=lo;
        for (const auto& b:blocks) { mean+=b[j]/5; lo=std::min(lo,b[j]); hi=std::max(hi,b[j]); }
        for (const auto& b:blocks) errors[j]+=(b[j]-mean)*(b[j]-mean)/20;
        errors[j]=std::sqrt(errors[j]);
        if (hi-lo>abs_tol+rel_tol*std::abs(mean)) out.block_stable=false;
    }
    out.se_cl=errors[0]; out.se_cd=errors[1]; out.se_cn=errors[2]; out.se_cm=errors[3];
    for (const auto& f:samples) out.cl_rms+=(f.Cl-out.mean.Cl)*(f.Cl-out.mean.Cl)/count;
    out.cl_rms=std::sqrt(out.cl_rms);
    if (std::abs(out.mean.Cn)<=std::max(1e-8,3*out.se_cn)) out.mean.x_cp=unavailable();
    return out;
}
}
