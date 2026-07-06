#pragma once

/// \file cube_grid.h
/// \brief Geometry and unit helpers for Gaussian CUBE grids.

#include <cstddef>

namespace oeio {
namespace cube {

/// CODATA Bohr radius in Angstrom (1 Bohr = 0.52917721067 A).
constexpr double BOHR_TO_ANGSTROM = 0.52917721067;

// Shared CUBE limits. These are the single source of truth for the bounds the
// reader enforces on an untrusted file AND the bounds the writer refuses to
// exceed, so a molecule/grid set the writer accepts always reads back (no
// round-trip asymmetry). Keeping them here rather than duplicated in each
// translation unit prevents the reader and writer from drifting apart.

/// Absolute tolerance (in the axis' native unit, i.e. Bohr for a Bohr file) for
/// treating a CUBE axis as aligned/cubic. The reader passes this to
/// is_axis_aligned_uniform(); the writer refuses a spacing that would serialize
/// at or below it (and so read back as "not axis-aligned").
constexpr double AXIS_TOL = 1e-6;

/// Maximum voxel count per grid dimension; bounds per-dimension allocation.
constexpr int MAX_VOXELS_PER_DIM = 8192;

/// Maximum total bytes across all grid buffers (bounds the whole allocation by
/// memory footprint). The derived element ceiling stays well below UINT_MAX, so
/// OEScalarGrid::SetValues() length casts never truncate.
constexpr std::size_t MAX_TOTAL_BYTES = 512u * 1024u * 1024u;  // 512 MiB
constexpr std::size_t MAX_TOTAL_ELEMENTS = MAX_TOTAL_BYTES / sizeof(float);

/// Maximum orbital (grid) count for a multi-grid MO cube.
constexpr int MAX_ORBITAL_COUNT = 1024;

/// Highest supported atomic number (Oganesson, Z=118). Atom records outside
/// [1, MAX_ATOMIC_NUMBER] are rejected on both read and write.
constexpr int MAX_ATOMIC_NUMBER = 118;

/// Maximum atom count. Bounds the atom loop so an absurd count fails cleanly.
constexpr long MAX_ATOM_COUNT = 100000000L;

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

/// \brief Validate the scalar shape of a CUBE record against the shared limits.
///
/// Checks the counts and dimensions (grid/orbital count, atom count, atomless
/// multi-grid encoding, per-dimension voxel bounds, and the total element
/// ceiling) and throws oeio::FormatError describing the first violation. This
/// is pure integer logic with no allocation, so the ceilings a live
/// molecule/grid could not reach without exhausting memory (the atom-count and
/// total-element limits) remain unit-testable. Callers that write a file must
/// invoke this BEFORE opening the target so a rejected write never truncates an
/// existing file.
///
/// \param natom Number of atoms in the molecule.
/// \param ngrid Number of grids to serialize (>= 1).
/// \param nx X voxel count of the (shared) grid geometry.
/// \param ny Y voxel count.
/// \param nz Z voxel count.
/// \raises oeio::FormatError When any scalar shape constraint is violated.
void validate_cube_shape(unsigned long natom, std::size_t ngrid,
                         int nx, int ny, int nz);

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
