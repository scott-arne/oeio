/// \file cube_grid.cpp
/// \brief Implementation of CUBE geometry and unit helpers.

#include "oeio/cube_grid.h"

#include "oeio/exceptions.h"

#include <cmath>
#include <string>

namespace oeio {
namespace cube {

void validate_cube_shape(unsigned long natom, std::size_t ngrid,
                         int nx, int ny, int nz) {
    if (ngrid == 0) {
        throw FormatError("oeio: CUBE requires at least one grid");
    }
    if (ngrid > static_cast<std::size_t>(MAX_ORBITAL_COUNT)) {
        throw FormatError("oeio: CUBE: grid count exceeds MO orbital limit (" +
                          std::to_string(MAX_ORBITAL_COUNT) + ")");
    }
    if (natom > static_cast<unsigned long>(MAX_ATOM_COUNT)) {
        throw FormatError("oeio: CUBE: atom count exceeds limit (" +
                          std::to_string(MAX_ATOM_COUNT) + ")");
    }
    // A multi-grid CUBE is encoded as an MO cube via a NEGATIVE atom count. With
    // zero atoms, negating yields 0, which the reader treats as a single-grid
    // cube and then misparses the orbital header as data. There is no valid MO
    // encoding for an atomless molecule, so reject it.
    if (ngrid > 1 && natom == 0) {
        throw FormatError(
            "oeio: CUBE: cannot write a multi-grid (MO) cube for a molecule "
            "with no atoms");
    }
    if (nx < 1 || ny < 1 || nz < 1) {
        throw FormatError("oeio: CUBE: grid dimensions must be positive");
    }
    if (nx > MAX_VOXELS_PER_DIM || ny > MAX_VOXELS_PER_DIM ||
        nz > MAX_VOXELS_PER_DIM) {
        throw FormatError("oeio: CUBE: voxel count per dimension must be <= " +
                          std::to_string(MAX_VOXELS_PER_DIM));
    }
    // Bound the total element count (voxels * ngrid) to the reader's ceiling.
    // The division form is overflow-safe: ngrid >= 1, so no product is formed.
    const std::size_t voxels =
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
        static_cast<std::size_t>(nz);
    if (voxels > MAX_TOTAL_ELEMENTS / ngrid) {
        throw FormatError("oeio: CUBE: total grid element count exceeds limit (" +
                          std::to_string(MAX_TOTAL_ELEMENTS) + " floats)");
    }
}

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
