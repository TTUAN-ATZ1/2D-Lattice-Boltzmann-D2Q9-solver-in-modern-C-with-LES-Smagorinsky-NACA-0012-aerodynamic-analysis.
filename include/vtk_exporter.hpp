#pragma once
#include <string>
#include <vector>
#include "Grid.hpp"
#include "geometry_naca.hpp"
#include "physical_units.hpp"

namespace lbm
{

    class VTKExporter
    {
    public:
        PhysicalUnits units;
        explicit VTKExporter(const PhysicalUnits& physical_units);

        // Xuất file trường dòng 2D toàn miền (.vtk STRUCTURED_POINTS)
        void export_field_2d(const std::string &filename,
                             const Grid2D &grid,
                             const std::vector<NodeType> &node_type) const;
    };

}
