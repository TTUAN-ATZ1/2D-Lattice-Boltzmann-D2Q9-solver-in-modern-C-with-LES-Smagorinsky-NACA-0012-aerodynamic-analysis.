#include "solver.hpp"
#include "vtk_exporter.hpp"
#include "run_config.hpp"
#include "force_statistics.hpp"
#include "physical_units.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace fs=std::filesystem;

static std::ofstream output(const fs::path& path) {
    std::ofstream f;
    f.exceptions(std::ios::failbit|std::ios::badbit);
    f.open(path); f << std::setprecision(17);
    return f;
}
struct Record {
    double alpha;
    lbm::ForceSummary total;
    double ac_total=lbm::unavailable();
};

static double fit_ac(const std::vector<Record>& records,size_t first,size_t end) {
    std::vector<lbm::ACSample> data;
    double lo=1e100,hi=-1e100,error=0;
    for (size_t i=first;i<end;++i) {
        const auto& s=records[i].total;
        if (!s.block_stable) return lbm::unavailable();
        data.push_back({records[i].alpha,s.mean.Cn,s.mean.Cm});
        lo=std::min(lo,s.mean.Cn); hi=std::max(hi,s.mean.Cn); error=std::max(error,s.se_cn);
    }
    if (hi-lo<=std::max(1e-8,6*error)) return lbm::unavailable();
    return lbm::aerodynamic_center(data);
}
static void write_polar(const fs::path& path,const std::vector<Record>& records) {
    auto f=output(path);
    f << "alpha_deg,Fx_mean_lb,Fy_mean_lb,Mz_quarter_mean_lb,Cl,Cd,Cn,Cm_quarter,x_cp_over_c,x_cp_valid,x_ac_local_over_c,block_stable,samples,Cl_rms,Cl_block_SE,Cd_block_SE,Cn_block_SE,Cm_block_SE\n";
    for (const auto& r:records) {
        const auto& s=r.total; const auto& a=s.mean;
        f << r.alpha << ','
          << a.Fx << ',' << a.Fy << ',' << a.Mz << ',' << a.Cl << ',' << a.Cd << ','
          << a.Cn << ',' << a.Cm << ',' << a.x_cp << ',' << std::isfinite(a.x_cp) << ','
          << r.ac_total << ',' << s.block_stable << ',' << s.samples << ','
          << s.cl_rms << ',' << s.se_cl << ',' << s.se_cd << ',' << s.se_cn << ',' << s.se_cm << '\n';
    }
    f.close();
}
int main(int argc,char** argv) {
    const fs::path csv_dir("results_csv");
    const fs::path vtk_dir("results_vtk");
    try {
        std::string config_path;
        for (int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if (arg=="--help") {
                std::cout << "Usage: lbm_naca0012 [--config path/to/config.cfg]\n";
                return 0;
            }
            if (arg=="--config" && i+1<argc) config_path=argv[++i];
            else throw std::invalid_argument("Unknown/missing argument: "+arg);
        }
        RunConfig c;
        if (!config_path.empty()) c=RunConfig::load(config_path);
        else c.validate();

        const double xp=(c.nx-1)*c.x_fraction,yp=(c.ny-1)*c.y_fraction;
        const double nu=c.u_inf*c.chord/c.re,tau=0.5+3*nu,q_chord=0.5*c.u_inf*c.u_inf*c.chord;
        const lbm::PhysicalUnits units(c.chord,c.u_inf,c.re,c.real_velocity,c.real_nu,c.real_rho,c.real_pressure);
        for (double angle:c.angles) {
            lbm::NACAGeometry g(c.chord,angle,xp,yp,c.surface_points);
            for (const auto& p:g.surface_points)
                if(p.x<3||p.x>c.nx-4||p.y<3||p.y>c.ny-4)
                    throw std::invalid_argument("Geometry needs at least three cells of domain margin");
        }
        fs::create_directories(csv_dir);
        if (c.write_vtk) fs::create_directories(vtk_dir);
        { auto state=output(csv_dir/"status.txt"); state << "running\n"; }
        {
            auto meta=output(csv_dir/"metadata.txt");
#ifdef _OPENMP
            meta << "openmp_threads=" << omp_get_max_threads() << '\n';
#else
            meta << "openmp_threads=1\n";
#endif
#define META(field) meta << #field "=" << c.field << '\n'
            META(nx); META(ny); META(chord); META(u_inf); META(re); META(cs); META(surface_points); META(collision);
            META(steps); META(average_steps); META(sample_every); META(health_every);
            META(x_fraction); META(y_fraction); META(real_velocity); META(real_nu); META(real_rho); META(real_pressure);
            META(max_mach); META(max_density_deviation); META(stable_abs); META(stable_rel); META(min_block_times);
            META(ac_min); META(ac_max); META(write_vtk);
#undef META
            meta << "angles="; for (double a:c.angles) meta << a << ' '; meta << '\n';
            meta << "nu_lb=" << nu << "\ntau_molecular=" << tau << "\nomega_molecular=" << 1/tau
                 << "\nx_pivot_lb=" << xp << "\ny_pivot_lb=" << yp
                 << "\nchord_m=" << units.chord << "\ndx_m=" << units.dx << "\ndt_s=" << units.dt
                 << "\ntotal_convective_times=" << c.steps*c.u_inf/c.chord
                 << "\naveraging_convective_times=" << c.average_steps*c.u_inf/c.chord
                 << "\nwall_boundary=interpolated_bounce_back_bouzidi_linear\n"
                 << "force_reference=quarter_chord\nmoment_sign=counterclockwise_positive\n";
            meta.close();
        }
        std::cout << "CSV output: " << fs::absolute(csv_dir).string()
                  << "\nVTK output: " << fs::absolute(vtk_dir).string()
                  << "\nRe=" << c.re << " nu_lb=" << nu << " tau=" << tau << " Cs=" << c.cs
                  << " Collision=" << c.collision << '\n';

        lbm::VTKExporter vtk(units);
        std::vector<Record> records;
        for (size_t ai=0;ai<c.angles.size();++ai) {
            const double angle=c.angles[ai];
            std::ostringstream id; id << "angle_" << std::setfill('0') << std::setw(3) << ai;
            const auto prefix=id.str();
            lbm::NACAGeometry geom(c.chord,angle,xp,yp,c.surface_points);
            lbm::LBMSolver solver(c.nx,c.ny,c.u_inf,c.re,c.chord,c.cs,lbm::OuterBoundary::Channel,
                c.collision=="mrt"?lbm::CollisionModel::MRT:lbm::CollisionModel::BGK);
            solver.set_geometry(geom);
            lbm::MomentumExchange mea(c.chord,c.u_inf,xp,yp);
            mea.build_boundary_links(c.nx,c.ny,solver.node_type,geom);
            auto history=output(csv_dir/(prefix+"_history.csv"));
            auto health=output(csv_dir/(prefix+"_health.csv"));
            history << "step,t_convective,alpha_deg,in_average,Cl,Cd,Cn,Cm\n";
            health << "step,t_convective,mass,rho_min,rho_max,max_mach,inlet_flux,outlet_flux,min_population\n";
            std::vector<lbm::AeroForce> totals;
            const int warmup=c.steps-c.average_steps;
            for (int step=1;step<=c.steps;++step) {
                try { solver.step(step%c.sample_every==0?&mea:nullptr); }
                catch (const std::exception& e) {
                    throw std::runtime_error("alpha="+std::to_string(angle)+", step="+std::to_string(step)+": "+e.what());
                }
                if (step%c.health_every==0 || step==1 || step==c.steps) {
                    solver.refresh_macroscopic();
                    const auto h=solver.health();
                    health << step << ',' << step*c.u_inf/c.chord << ',' << h.mass << ',' << h.rho_min << ','
                           << h.rho_max << ',' << h.max_mach << ',' << h.inlet_flux << ',' << h.outlet_flux << ',' << h.min_population << '\n';
                    if (h.max_mach>c.max_mach || std::max(std::abs(h.rho_min-1),std::abs(h.rho_max-1))>c.max_density_deviation) {
                        health.flush();
                        history.flush();
                        int nmin=-1,nmax=-1,nmach=-1;
                        for (int n=0;n<solver.grid.size;++n) {
                            if (solver.node_type[n]==lbm::NODE_SOLID) continue;
                            if (nmin<0 || solver.grid.rho[n]<solver.grid.rho[nmin]) nmin=n;
                            if (nmax<0 || solver.grid.rho[n]>solver.grid.rho[nmax]) nmax=n;
                            if (nmach<0 || std::hypot(solver.grid.ux[n],solver.grid.uy[n])>
                                std::hypot(solver.grid.ux[nmach],solver.grid.uy[nmach])) nmach=n;
                        }
                        std::ostringstream failure;
                        failure << std::setprecision(17)
                            << "Flow-health limit exceeded at alpha=" << angle << ", step=" << step << '\n';
                        if (1-h.rho_min>c.max_density_deviation)
                            failure << "TRIGGER: minimum density below " << 1-c.max_density_deviation << '\n';
                        if (h.rho_max-1>c.max_density_deviation)
                            failure << "TRIGGER: maximum density above " << 1+c.max_density_deviation << '\n';
                        if (h.max_mach>c.max_mach)
                            failure << "TRIGGER: maximum Mach above " << c.max_mach << '\n';
                        failure << "rho_min=" << h.rho_min << " at x=" << nmin%c.nx << ", y=" << nmin/c.nx << '\n'
                            << "rho_max=" << h.rho_max << " at x=" << nmax%c.nx << ", y=" << nmax/c.nx << '\n'
                            << "max_mach=" << h.max_mach << " at x=" << nmach%c.nx << ", y=" << nmach/c.nx;
                        fs::create_directories(vtk_dir);
                        const auto snapshot=vtk_dir/(prefix+"_failure_step_"+std::to_string(step)+".vtk");
                        try {
                            vtk.export_field_2d(snapshot.string(),solver.grid,solver.node_type);
                            failure << "\nFailure field saved: " << snapshot.string();
                        } catch (const std::exception& e) {
                            failure << "\nCould not save failure field: " << e.what();
                        }
                        throw std::runtime_error(failure.str());
                    }
                }
                if (step%c.sample_every==0) {
                    const auto& f=solver.last_mea_force;
                    history << step << ',' << step*c.u_inf/c.chord << ',' << angle << ',' << (step>warmup) << ','
                            << f.Cl << ',' << f.Cd << ',' << f.Cn << ',' << f.Cm << '\n';
                    if (step>warmup) {
                        totals.push_back(f);
                    }
                }
                if (step%std::max(c.health_every,c.steps/10)==0)
                    std::cout << "alpha=" << angle << " step=" << step << '/' << c.steps << '\n' << std::flush;
            }
            history.close(); health.close();
            const bool enough_time=c.average_steps*c.u_inf/(5*c.chord)>=c.min_block_times;
            auto total=lbm::summarize_forces(totals,angle*lbm::PI/180,c.chord,q_chord,c.stable_abs,c.stable_rel,enough_time);
            records.push_back({angle,total});
            if (c.write_vtk) {
                solver.refresh_macroscopic();
                vtk.export_field_2d((vtk_dir/(prefix+"_flow_final.vtk")).string(),solver.grid,solver.node_type);
            }
            write_polar(csv_dir/"polar.csv",records);
            std::cout << "alpha=" << angle << " mean Cl=" << total.mean.Cl << " Cd=" << total.mean.Cd
                      << " Cm=" << total.mean.Cm << " xcp/c=" << total.mean.x_cp
                      << " block_stable=" << total.block_stable << '\n';
        }
        if (records.size()>=3) {
            for (size_t i=0;i<records.size();++i) {
                const size_t first=i==0?0:i==records.size()-1?i-2:i-1;
                if (records[first].alpha<c.ac_min || records[first+2].alpha>c.ac_max) continue;
                records[i].ac_total=fit_ac(records,first,first+3);
            }
        }
        write_polar(csv_dir/"polar.csv",records);
        auto ac=output(csv_dir/"ac_fit.csv");
        ac << "method,alpha_min_requested,alpha_max_requested,angles_used,x_ac_over_c\n";
        size_t first=0,end=0;
        while (first<records.size()&&records[first].alpha<c.ac_min) ++first;
        end=first;
        while (end<records.size()&&records[end].alpha<=c.ac_max) ++end;
        ac << "MEA_IBB," << c.ac_min << ',' << c.ac_max << ','
           << end-first << ',' << fit_ac(records,first,end) << '\n';
        ac.close();
        auto state=output(csv_dir/"status.txt");
        state << "completed\n";
        state.close();
        std::cout << "Completed. Results written to results_csv/ and results_vtk/.\n";
        return 0;
    } catch (const std::exception& e) {
        try {
            fs::create_directories(csv_dir);
            auto state=output(csv_dir/"status.txt");
            state << "failed\n" << e.what() << '\n';
            state.close();
        } catch (...) {}
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
