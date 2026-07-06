/// \file cube_writer.cpp
/// \brief CUBE serialization: OEMol + N OEScalarGrid -> file.

#include "oeio/cube_handler.h"

#include "oeio/cube_grid.h"
#include "oeio/exceptions.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
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

    std::ofstream out(path_);
    if (!out) throw FileError("oeio: unable to open '" + path_ + "' for writing");

    // Emit enough significant digits to round-trip IEEE-754 single precision
    // exactly; the default (~6) would silently truncate non-round coordinates
    // and grid values. max_digits10 for float is 9.
    out << std::setprecision(std::numeric_limits<float>::max_digits10);

    const OESystem::OEScalarGrid& g0 = *grids[0];
    const int nx = static_cast<int>(g0.GetXDim());
    const int ny = static_cast<int>(g0.GetYDim());
    const int nz = static_cast<int>(g0.GetZDim());
    const int ngrid = static_cast<int>(grids.size());

    // Convert Angstrom geometry back to Bohr for the file.
    const double inv = 1.0 / cube::BOHR_TO_ANGSTROM;
    const double spacing_bohr = g0.GetSpacing() * inv;
    float midx, midy, midz;
    g0.GetMid(midx, midy, midz);
    const double mid_ang[3] = { midx, midy, midz };
    const int nvox[3] = { nx, ny, nz };
    double origin_ang[3];
    cube::mid_to_origin(mid_ang, nvox, g0.GetSpacing(), origin_ang);
    const double origin_bohr[3] = { origin_ang[0] * inv, origin_ang[1] * inv,
                                    origin_ang[2] * inv };

    out << "CUBE file written by oeio\n";
    out << "OEScalarGrid volumetric data\n";

    // Atom count: negative if multi-grid (MO cube).
    const int natom = mol.NumAtoms();
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

    // Volumetric block: orbital-fastest then z,y,x.
    const std::size_t voxels = static_cast<std::size_t>(nx) * ny * nz;
    std::vector<const float*> data(ngrid);
    for (int g = 0; g < ngrid; ++g) data[g] = grids[g]->GetValues();
    for (std::size_t vx = 0; vx < voxels; ++vx) {
        for (int g = 0; g < ngrid; ++g) {
            out << " " << data[g][vx];
        }
        out << "\n";
    }

    return static_cast<bool>(out);
}

void CubeMolSink::close() {}

}  // namespace builtin
}  // namespace oeio
