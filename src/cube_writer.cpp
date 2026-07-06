/// \file cube_writer.cpp
/// \brief CUBE serialization: OEMol + N OEScalarGrid -> file.

#include "oeio/cube_handler.h"

#include "oeio/cube_grid.h"
#include "oeio/exceptions.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace oeio {
namespace builtin {

namespace {
constexpr double GEOM_TOL = 1e-4;

bool same_geometry(const OESystem::OEScalarGrid& a, const OESystem::OEScalarGrid& b) {
    if (a.GetXDim() != b.GetXDim() || a.GetYDim() != b.GetYDim() ||
        a.GetZDim() != b.GetZDim()) return false;
    if (std::fabs(a.GetSpacing() - b.GetSpacing()) > GEOM_TOL) return false;
    return std::fabs(a.GetXMin() - b.GetXMin()) <= GEOM_TOL &&
           std::fabs(a.GetYMin() - b.GetYMin()) <= GEOM_TOL &&
           std::fabs(a.GetZMin() - b.GetZMin()) <= GEOM_TOL;
}

/// Throw FormatError unless the value is finite (no NaN/Inf).
void require_finite(double value, const char* what) {
    if (!std::isfinite(value)) {
        throw FormatError(std::string("oeio: CUBE: non-finite ") + what);
    }
}
}  // namespace

CubeMolSink::CubeMolSink(const std::string& path) : path_(path) {}

bool CubeMolSink::write(const OEChem::OEMolBase&) {
    throw FormatError("oeio: CUBE requires at least one grid; use "
                      "append(mol, grid) / writer.append(mol, *grids)");
}

bool CubeMolSink::write(const OEChem::OEMolBase& mol,
                       const std::vector<const OESystem::OEScalarGrid*>& grids) {
    if (grids.empty()) {
        throw FormatError("oeio: CUBE requires at least one grid");
    }
    for (const auto* g : grids) {
        if (!g) throw FormatError("oeio: CUBE: null grid");
        if (!same_geometry(*grids[0], *g)) {
            throw FormatError(
                "oeio: CUBE: all grids must share geometry (dimensions, "
                "origin, and spacing)");
        }
    }

    // A CUBE file holds exactly one record. Reject a second write before opening
    // the file so a repeated append through the Writer API cannot truncate and
    // overwrite the first record (silent data loss). This mirrors the reader's
    // single-record contract.
    if (written_) {
        throw FormatError("oeio: CUBE holds a single record; cannot append again");
    }

    const int ngrid = static_cast<int>(grids.size());
    const int natom = mol.NumAtoms();

    // Bound the grid and atom counts to the same ceilings the reader enforces,
    // so a molecule/grid set the writer accepts always reads back. A multi-grid
    // set larger than MAX_ORBITAL_COUNT would emit an orbital header the reader
    // rejects; an atom count above MAX_ATOM_COUNT would exceed the reader's
    // limit (and NumAtoms() is tested unsigned to avoid the int truncation the
    // -natom marker below would otherwise suffer).
    if (ngrid > cube::MAX_ORBITAL_COUNT) {
        throw FormatError("oeio: CUBE: grid count exceeds MO orbital limit (" +
                          std::to_string(cube::MAX_ORBITAL_COUNT) + ")");
    }
    if (mol.NumAtoms() > static_cast<unsigned>(cube::MAX_ATOM_COUNT)) {
        throw FormatError("oeio: CUBE: atom count exceeds limit (" +
                          std::to_string(cube::MAX_ATOM_COUNT) + ")");
    }

    // A multi-grid CUBE is encoded as an MO cube: a NEGATIVE atom count signals
    // the reader to consume the orbital header line. With zero atoms, negating
    // the count yields 0, which the reader treats as a single-grid cube and then
    // misparses the orbital header as volumetric data. There is no valid MO
    // encoding for an atomless molecule, so reject it rather than emit a file
    // that will not round-trip. All input validation happens BEFORE the output
    // file is opened so a rejected write never truncates an existing file.
    if (ngrid > 1 && natom == 0) {
        throw FormatError(
            "oeio: CUBE: cannot write a multi-grid (MO) cube for a molecule "
            "with no atoms");
    }

    // Validate the ENTIRE serializable payload before opening the target, so an
    // invalid grid (zero/negative dimensions, non-finite spacing/coordinates/
    // values) is rejected without truncating an existing file and without
    // emitting a CUBE that CubeMolSource would refuse to read back.
    const OESystem::OEScalarGrid& g0 = *grids[0];
    const int nx = static_cast<int>(g0.GetXDim());
    const int ny = static_cast<int>(g0.GetYDim());
    const int nz = static_cast<int>(g0.GetZDim());
    if (nx < 1 || ny < 1 || nz < 1) {
        throw FormatError("oeio: CUBE: grid dimensions must be positive");
    }
    if (nx > cube::MAX_VOXELS_PER_DIM || ny > cube::MAX_VOXELS_PER_DIM ||
        nz > cube::MAX_VOXELS_PER_DIM) {
        throw FormatError("oeio: CUBE: voxel count per dimension must be <= " +
                          std::to_string(cube::MAX_VOXELS_PER_DIM));
    }
    const float spacing = g0.GetSpacing();
    // Reject a spacing that would serialize (in Bohr) at or below the reader's
    // axis tolerance, since is_axis_aligned_uniform() would then treat the axis
    // as degenerate and refuse the file. For an Angstrom-origin grid the file
    // stores spacing in Bohr, so the Bohr magnitude is what must clear the tol.
    const double spacing_bohr_check =
        static_cast<double>(spacing) / cube::BOHR_TO_ANGSTROM;
    if (!std::isfinite(spacing) || spacing_bohr_check <= cube::AXIS_TOL) {
        throw FormatError("oeio: CUBE: grid spacing must be finite and positive");
    }
    float midx, midy, midz;
    g0.GetMid(midx, midy, midz);
    require_finite(midx, "grid midpoint");
    require_finite(midy, "grid midpoint");
    require_finite(midz, "grid midpoint");
    for (OESystem::OEIter<const OEChem::OEAtomBase> ai = mol.GetAtoms(); ai; ++ai) {
        // The reader rejects an atomic number outside [1, MAX_ATOMIC_NUMBER]; a
        // dummy/wildcard atom (Z == 0, common in OEChem) or an out-of-range Z
        // would therefore write a record the reader refuses. Guard here so the
        // round-trip is symmetric.
        const int z = ai->GetAtomicNum();
        if (z < 1 || z > cube::MAX_ATOMIC_NUMBER) {
            throw FormatError("oeio: CUBE: atomic number out of range [1, " +
                              std::to_string(cube::MAX_ATOMIC_NUMBER) + "]");
        }
        float xyz[3];
        mol.GetCoords(&*ai, xyz);
        require_finite(xyz[0], "atom coordinate");
        require_finite(xyz[1], "atom coordinate");
        require_finite(xyz[2], "atom coordinate");
    }
    const std::size_t voxels = static_cast<std::size_t>(nx) * ny * nz;
    // Bound the total element count (voxels * ngrid) to the reader's ceiling.
    // The division form is overflow-safe: ngrid >= 1, so no product is formed.
    if (voxels > cube::MAX_TOTAL_ELEMENTS / static_cast<std::size_t>(ngrid)) {
        throw FormatError("oeio: CUBE: total grid element count exceeds limit (" +
                          std::to_string(cube::MAX_TOTAL_ELEMENTS) + " floats)");
    }
    for (int g = 0; g < ngrid; ++g) {
        const float* vals = grids[g]->GetValues();
        if (!vals) throw FormatError("oeio: CUBE: grid has no values");
        for (std::size_t vx = 0; vx < voxels; ++vx) {
            require_finite(vals[vx], "volumetric value");
        }
    }

    std::ofstream out(path_);
    if (!out) throw FileError("oeio: unable to open '" + path_ + "' for writing");

    // Emit enough significant digits to round-trip IEEE-754 single precision
    // exactly; the default (~6) would silently truncate non-round coordinates
    // and grid values. max_digits10 for float is 9.
    out << std::setprecision(std::numeric_limits<float>::max_digits10);

    // Convert Angstrom geometry back to Bohr for the file.
    const double inv = 1.0 / cube::BOHR_TO_ANGSTROM;
    const double spacing_bohr = spacing * inv;
    const double mid_ang[3] = { midx, midy, midz };
    const int nvox[3] = { nx, ny, nz };
    double origin_ang[3];
    cube::mid_to_origin(mid_ang, nvox, spacing, origin_ang);
    const double origin_bohr[3] = { origin_ang[0] * inv, origin_ang[1] * inv,
                                    origin_ang[2] * inv };

    out << "CUBE file written by oeio\n";
    out << "OEScalarGrid volumetric data\n";

    // Atom count: negative if multi-grid (MO cube).
    const int atom_count_field = (ngrid > 1) ? -natom : natom;
    out << atom_count_field << " " << origin_bohr[0] << " " << origin_bohr[1]
        << " " << origin_bohr[2] << "\n";

    // Axis vectors: positive diagonal = spacing (Bohr), zero off-diagonal.
    out << nx << " " << spacing_bohr << " 0.0 0.0\n";
    out << ny << " 0.0 " << spacing_bohr << " 0.0\n";
    out << nz << " 0.0 0.0 " << spacing_bohr << "\n";

    // Atoms: atomic number, charge (use atomic number as nuclear charge), x,y,z (Bohr).
    for (OESystem::OEIter<const OEChem::OEAtomBase> ai = mol.GetAtoms(); ai; ++ai) {
        float xyz[3];
        mol.GetCoords(&*ai, xyz);
        const int z = ai->GetAtomicNum();
        out << z << " " << static_cast<double>(z) << " "
            << xyz[0] * inv << " " << xyz[1] * inv << " " << xyz[2] * inv << "\n";
    }

    // MO orbital header line for multi-grid.
    if (ngrid > 1) {
        out << ngrid;
        for (int g = 1; g <= ngrid; ++g) out << " " << g;
        out << "\n";
    }

    // Volumetric block: orbital-fastest then z,y,x. `voxels` was computed during
    // pre-open validation above.
    std::vector<const float*> data(ngrid);
    for (int g = 0; g < ngrid; ++g) data[g] = grids[g]->GetValues();
    for (std::size_t vx = 0; vx < voxels; ++vx) {
        for (int g = 0; g < ngrid; ++g) {
            out << " " << data[g][vx];
        }
        out << "\n";
    }

    // Close explicitly and re-check the stream state so a failure that only
    // surfaces during the final flush/close (disk full, quota, network
    // filesystem) is observed here rather than being swallowed by the
    // destructor after success was already reported. Mark the single record
    // consumed only after this close-checked success, so a failed write leaves
    // the sink retryable.
    out.close();
    const bool ok = static_cast<bool>(out);
    if (ok) written_ = true;
    return ok;
}

void CubeMolSink::close() {}

}  // namespace builtin
}  // namespace oeio
