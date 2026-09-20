#pragma once
#include <vector>
#include "D2Q9.hpp"
#include <cstdint>
namespace lbm
{
    enum NodeType : uint8_t
    {
        NODE_FLUID = 0,   // dong chat luu
        NODE_SOLID = 1,   // Ruot canh
        NODE_BOUNDARY = 2 // Nut chat luu tiep xuc vo canh
    };
    struct Point2D
    {
        double x, y;
    };
    class NACAGeometry
    {
    public:
        double chord;
        double alpha_deg;
        double x_pivot, y_pivot;
        int Nb; //

        std::vector<Point2D> surface_points;
        std::vector<Point2D> raw_points;
        std::vector<double> ds;
        std::vector<Point2D> normal;

        NACAGeometry(double chord_len, double alpha, double xp, double yp, int num_points = 120);

        // tao hinh va xoay goc tan
        void build_profile();
        // phan loai nut
        void classify_grid(int nx, int ny, std::vector<NodeType> &node_type, std::vector<double> &wall_dist) const;

    private:
        // tinh khoang cach ngan nhat tu 1 diem den A,B
        static double distance_to_segment(double px, double py, double ax, double ay, double bx, double by);
    };
}