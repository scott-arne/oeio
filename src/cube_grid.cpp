/// \file cube_grid.cpp
/// \brief Implementation of CUBE geometry and unit helpers.

#include "oeio/cube_grid.h"

#include <cmath>

namespace oeio {
namespace cube {

bool is_axis_aligned_uniform(const CubeAxes& ax, double tol, double& spacing_out) {
    if (!std::isfinite(tol)) return false;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!std::isfinite(ax.vec[i][j])) return false;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            if (std::fabs(ax.vec[i][j]) > tol) return false;  // off-diagonal -> skew
        }
    }
    const double s0 = ax.vec[0][0];
    const double s1 = ax.vec[1][1];
    const double s2 = ax.vec[2][2];
    if (s0 <= tol || s1 <= tol || s2 <= tol) return false;    // reflected/zero
    if (std::fabs(s0 - s1) > tol || std::fabs(s0 - s2) > tol) return false;  // anisotropic
    spacing_out = s0;
    return true;
}

void origin_to_mid(const double origin[3], const int nvox[3], double spacing,
                   double mid_out[3]) {
    for (int i = 0; i < 3; ++i) {
        mid_out[i] = origin[i] + spacing * ((nvox[i] - 1) / 2.0);
    }
}

void mid_to_origin(const double mid[3], const int nvox[3], double spacing,
                   double origin_out[3]) {
    for (int i = 0; i < 3; ++i) {
        origin_out[i] = mid[i] - spacing * ((nvox[i] - 1) / 2.0);
    }
}

}  // namespace cube
}  // namespace oeio
