#pragma once

/// \file cube_grid.h
/// \brief Geometry and unit helpers for Gaussian CUBE grids.

namespace oeio {
namespace cube {

/// CODATA Bohr radius in Angstrom (1 Bohr = 0.52917721067 A).
constexpr double BOHR_TO_ANGSTROM = 0.52917721067;

/// Raw CUBE header geometry (origin, voxel counts, and the three axis vectors),
/// in the unit the header used (converted to Angstrom by the caller afterward).
struct CubeAxes {
    double origin[3];
    int nvox[3];
    double vec[3][3];  // vec[i] is the step vector along grid dimension i
};

/// \brief Test whether the axes describe an axis-aligned, positively-oriented,
/// cubic-voxel grid, and if so report the common spacing.
///
/// Requires each axis vector parallel to its Cartesian axis (off-diagonal
/// magnitude <= tol), positive diagonal (> tol), and equal diagonal magnitudes
/// across all three axes (within tol).
///
/// \param ax The raw axis geometry.
/// \param tol Absolute tolerance for the checks.
/// \param spacing_out Set to the common axis length when the function returns true.
/// \returns true iff the grid is axis-aligned/positive/cubic.
bool is_axis_aligned_uniform(const CubeAxes& ax, double tol, double& spacing_out);

/// \brief Convert a CUBE origin (corner) to an OEScalarGrid midpoint.
///
/// mid = origin + spacing * ((nx-1)/2, (ny-1)/2, (nz-1)/2)
void origin_to_mid(const double origin[3], const int nvox[3], double spacing,
                   double mid_out[3]);

/// \brief Convert an OEScalarGrid midpoint to a CUBE origin (corner).
///
/// origin = mid - spacing * ((nx-1)/2, (ny-1)/2, (nz-1)/2)
void mid_to_origin(const double mid[3], const int nvox[3], double spacing,
                   double origin_out[3]);

}  // namespace cube
}  // namespace oeio
