#pragma once
#include <array>

namespace lbm
{
    constexpr int Q = 9;
    constexpr double CS2 = 1.0 / 3.0;
    constexpr double INV_CS2 = 3.0;
    constexpr double INV_2CS4 = 4.5;
    constexpr double INV_2CS2 = 1.5;

    constexpr std::array<int, Q> EX = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    constexpr std::array<int, Q> EY = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    constexpr std::array<double, Q> WEIGHTS = {4.0 / 9.0,
                                               1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                                               1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    constexpr std::array<int, Q> OPP = {0, 3, 4, 1, 2, 7, 8, 5, 6};
}