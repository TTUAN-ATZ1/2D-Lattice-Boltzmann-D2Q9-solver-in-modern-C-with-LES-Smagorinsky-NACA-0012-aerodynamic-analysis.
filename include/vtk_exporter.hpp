#pragma once
#include <string>
#include <vector>
#include "Grid.hpp"
#include "geometry_naca.hpp"

namespace lbm
{

    class VTKExporter
    {
    public:
        double u_inf_lb;
        double rho_real;
        double V_real;
        double P_inf_real;
        double q_real;

        VTKExporter(double u_inf, double rho_r = 1.225, double V_r = 30.0, double P_inf_r = 101325.0);

        // Xuất file trường dòng 2D toàn miền (.vtk STRUCTURED_POINTS)
        void export_field_2d(const std::string &filename,
                             const Grid2D &grid,
                             const std::vector<NodeType> &node_type) const;

        // Xuất file viền cánh 120 điểm mốc (.vtk POLYDATA)
        void export_airfoil_surface(const std::string &filename,
                                    const NACAGeometry &geom,
                                    const std::vector<double> &cp_values) const;
    };

}