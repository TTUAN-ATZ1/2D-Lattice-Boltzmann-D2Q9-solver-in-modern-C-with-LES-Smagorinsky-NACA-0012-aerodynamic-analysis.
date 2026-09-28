#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct RunConfig {
    int nx=640,ny=360,surface_points=120,steps=30000,average_steps=15000;
    int sample_every=10,health_every=100;
    double chord=150,u_inf=0.04,re=100000,cs=0.16,x_fraction=0.25,y_fraction=0.5;
    double real_velocity=3.0,real_nu=1.5e-5,real_rho=1.225,real_pressure=101325;
    double max_mach=0.3,max_density_deviation=0.1;
    double stable_abs=0.005,stable_rel=0.02,min_block_times=1;
    double ac_min=0,ac_max=10;
    bool write_vtk=true;
    std::string collision="bgk";
    std::vector<double> angles{0, 1, 2, 4, 6, 8, 10, 12, 14, 16};

    static std::string trim(const std::string& s) {
        const auto b=s.find_first_not_of(" \t\r\n");
        return b==std::string::npos?"":s.substr(b,s.find_last_not_of(" \t\r\n")-b+1);
    }
    static double number(const std::string& s) {
        size_t n=0; const double v=std::stod(s,&n);
        if (n!=s.size()||!std::isfinite(v)) throw std::invalid_argument("Invalid number: "+s);
        return v;
    }
    static RunConfig load(const std::string& path) {
        RunConfig c;
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open configuration: "+path);
        std::map<std::string,int*> ints{{"nx",&c.nx},{"ny",&c.ny},{"surface_points",&c.surface_points},
            {"steps",&c.steps},{"average_steps",&c.average_steps},{"sample_every",&c.sample_every},
            {"health_every",&c.health_every}};
        std::map<std::string,double*> doubles{{"chord",&c.chord},{"u_inf",&c.u_inf},{"re",&c.re},{"cs",&c.cs},
            {"x_fraction",&c.x_fraction},{"y_fraction",&c.y_fraction},{"real_velocity",&c.real_velocity},
            {"real_nu",&c.real_nu},{"real_rho",&c.real_rho},{"real_pressure",&c.real_pressure},
            {"max_mach",&c.max_mach},{"max_density_deviation",&c.max_density_deviation},
            {"stable_abs",&c.stable_abs},{"stable_rel",&c.stable_rel},{"min_block_times",&c.min_block_times},
            {"ac_min",&c.ac_min},{"ac_max",&c.ac_max}};
        std::set<std::string> seen;
        std::string line;
        while (std::getline(in,line)) {
            line=trim(line.substr(0,line.find('#')));
            if (line.empty()) continue;
            const auto eq=line.find('=');
            if (eq==std::string::npos) throw std::invalid_argument("Expected key=value: "+line);
            const auto k=trim(line.substr(0,eq)),v=trim(line.substr(eq+1));
            if (!seen.insert(k).second) throw std::invalid_argument("Duplicate setting: "+k);
            if (ints.count(k)) {
                const double value=number(v);
                if (value<1||value>1000000000||std::floor(value)!=value)
                    throw std::invalid_argument("Expected positive integer: "+k);
                *ints.at(k)=static_cast<int>(value);
            } else if (doubles.count(k)) *doubles.at(k)=number(v);
            else if (k=="angles") {
                c.angles.clear(); std::istringstream list(v); std::string item;
                while (std::getline(list,item,',')) c.angles.push_back(number(trim(item)));
            } else if (k=="collision") c.collision=v;
            else if (k=="write_vtk") {
                if (v!="true"&&v!="false") throw std::invalid_argument("write_vtk must be true or false");
                c.write_vtk=v=="true";
            } else throw std::invalid_argument("Unknown setting: "+k);
        }
        c.validate(); return c;
    }
    void validate() const {
        if (collision!="bgk" && collision!="mrt") throw std::invalid_argument("collision must be bgk or mrt");
        if (collision=="mrt" && cs!=0) throw std::invalid_argument("MRT requires cs=0 until a compatible SGS closure is validated");
        if (nx<8||ny<8||chord<8||u_inf<=0||u_inf>=0.15||re<=0||cs<0
            ||surface_points<8||surface_points%2!=0||average_steps>steps
            ||sample_every<1||health_every<1||average_steps/sample_every<20
            ||(steps-average_steps)%sample_every!=0||average_steps%sample_every!=0
            ||x_fraction<=0||x_fraction>=1||y_fraction<=0||y_fraction>=1
            ||max_mach<=0||max_density_deviation<=0||stable_abs<=0||stable_rel<0
            ||min_block_times<=0||ac_min>ac_max||angles.empty())
            throw std::invalid_argument("Invalid grid/run settings; averaging needs >=20 samples and aligned sampling intervals");
        for (size_t i=0;i<angles.size();++i)
            if (std::abs(angles[i])>=90 || (i>0 && angles[i]<=angles[i-1]))
                throw std::invalid_argument("Angles must be strictly increasing and between -90 and 90 degrees");
    }
};
