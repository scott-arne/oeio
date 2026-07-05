/// \file cube_reader.cpp
/// \brief CUBE parsing: file -> OEMol + N OEScalarGrid.

#include "oeio/cube_handler.h"

#include "oeio/cube_grid.h"
#include "oeio/exceptions.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

namespace oeio {
namespace builtin {

namespace {
constexpr double AXIS_TOL = 1e-6;

/// Maximum allowed voxel count per dimension; prevents absurd allocations.
constexpr int MAX_VOXELS_PER_DIM = 8192;

/// Maximum total voxel count (nx*ny*nz); must fit in unsigned int for SetValues.
constexpr std::size_t MAX_TOTAL_VOXELS = 536870912u;  // 512 MB @ 4 bytes/float

/// Maximum orbital count for MO cubes; prevents huge grid arrays.
constexpr int MAX_ORBITAL_COUNT = 1024;

/// Read the next whitespace-delimited token stream line, throwing on failure.
std::vector<double> parse_doubles(std::istream& in, int count, const char* what) {
    std::vector<double> out;
    out.reserve(count);
    double v;
    for (int i = 0; i < count; ++i) {
        if (!(in >> v)) {
            throw FormatError(std::string("oeio: CUBE: malformed ") + what);
        }
        out.push_back(v);
    }
    return out;
}
}  // namespace

CubeMolSource::CubeMolSource(const std::string& path) : path_(path) {}

bool CubeMolSource::next(OEChem::OEGraphMol& mol) {
    std::vector<OESystem::OEScalarGrid> grids;  // discarded
    return read_record(mol, grids);
}

bool CubeMolSource::next(OEChem::OEMolBase& mol,
                        const std::vector<OESystem::OEScalarGrid*>& grids,
                        int* num_grids) {
    std::vector<OESystem::OEScalarGrid> owned;
    if (!read_record(mol, owned)) return false;
    if (num_grids) *num_grids = static_cast<int>(owned.size());
    const std::size_t k = std::min(grids.size(), owned.size());
    for (std::size_t i = 0; i < k; ++i) {
        if (grids[i]) *grids[i] = owned[i];  // copy-assign into caller's grid
    }
    return true;
}

bool CubeMolSource::next(OEChem::OEMolBase& mol,
                        std::vector<OESystem::OEScalarGrid>& grids) {
    return read_record(mol, grids);
}

bool CubeMolSource::read_record(OEChem::OEMolBase& mol,
                               std::vector<OESystem::OEScalarGrid>& grids) {
    mol.Clear();
    grids.clear();
    if (consumed_) return false;
    consumed_ = true;

    std::ifstream in(path_);
    if (!in) {
        throw FileError("oeio: unable to open '" + path_ + "' for reading");
    }

    std::string comment1, comment2;
    std::getline(in, comment1);
    std::getline(in, comment2);

    // Line 3: atom count + origin.
    long natom_raw = 0;
    if (!(in >> natom_raw)) throw FormatError("oeio: CUBE: missing atom count");
    const bool mo_cube = natom_raw < 0;
    const long natom = std::labs(natom_raw);
    auto origin = parse_doubles(in, 3, "origin");

    // Lines 4-6: voxel counts + axis vectors.
    cube::CubeAxes ax{};
    ax.origin[0] = origin[0]; ax.origin[1] = origin[1]; ax.origin[2] = origin[2];
    bool angstrom = false;
    for (int i = 0; i < 3; ++i) {
        long nv = 0;
        if (!(in >> nv)) throw FormatError("oeio: CUBE: missing voxel count");
        if (nv < 0) angstrom = true;   // negative voxel count -> Angstrom
        const long nv_abs = std::labs(nv);
        if (nv_abs < 1 || nv_abs > MAX_VOXELS_PER_DIM) {
            throw FormatError("oeio: CUBE: voxel count per dimension must be in [1, " +
                              std::to_string(MAX_VOXELS_PER_DIM) + "]");
        }
        ax.nvox[i] = static_cast<int>(nv_abs);
        auto row = parse_doubles(in, 3, "axis vector");
        ax.vec[i][0] = row[0]; ax.vec[i][1] = row[1]; ax.vec[i][2] = row[2];
    }

    double spacing = 0.0;
    if (!cube::is_axis_aligned_uniform(ax, AXIS_TOL, spacing)) {
        throw FormatError(
            "oeio: CUBE: only axis-aligned, positively-oriented, cubic-voxel "
            "grids are supported (v1); rotated/skewed/reflected/anisotropic "
            "axes rejected");
    }

    // Unit conversion factor to Angstrom (positive voxel count -> Bohr).
    const double to_ang = angstrom ? 1.0 : cube::BOHR_TO_ANGSTROM;

    // Atoms: atomic number, charge, x, y, z. Build the molecule (coords in A).
    for (long a = 0; a < natom; ++a) {
        auto atom_line = parse_doubles(in, 5, "atom line");
        OEChem::OEAtomBase* atom = mol.NewAtom(static_cast<unsigned int>(atom_line[0]));
        const float coords[3] = {
            static_cast<float>(atom_line[2] * to_ang),
            static_cast<float>(atom_line[3] * to_ang),
            static_cast<float>(atom_line[4] * to_ang)
        };
        mol.SetCoords(atom, coords);
    }

    // MO cube: extra header line "M id1 id2 ... idM".
    int ngrid = 1;
    if (mo_cube) {
        long m = 0;
        if (!(in >> m)) throw FormatError("oeio: CUBE: missing orbital count");
        if (m < 1 || m > MAX_ORBITAL_COUNT) {
            throw FormatError("oeio: CUBE: orbital count must be in [1, " +
                              std::to_string(MAX_ORBITAL_COUNT) + "]");
        }
        ngrid = static_cast<int>(m);
        for (long i = 0; i < m; ++i) {
            long id;
            if (!(in >> id)) throw FormatError("oeio: CUBE: malformed orbital ids");
        }
    }

    // Grid geometry in Angstrom.
    const double spacing_ang = spacing * to_ang;
    double origin_ang[3] = { ax.origin[0] * to_ang, ax.origin[1] * to_ang,
                             ax.origin[2] * to_ang };
    double mid[3];
    cube::origin_to_mid(origin_ang, ax.nvox, spacing_ang, mid);

    const std::size_t voxels =
        static_cast<std::size_t>(ax.nvox[0]) * ax.nvox[1] * ax.nvox[2];

    // Validate total voxel count to prevent overflow and absurd allocations.
    if (voxels > MAX_TOTAL_VOXELS) {
        throw FormatError("oeio: CUBE: total voxel count (" + std::to_string(voxels) +
                          ") exceeds limit (" + std::to_string(MAX_TOTAL_VOXELS) + ")");
    }
    // Check that total element count (voxels * ngrid) won't overflow.
    if (ngrid > 1 && voxels > MAX_TOTAL_VOXELS / static_cast<std::size_t>(ngrid)) {
        throw FormatError("oeio: CUBE: total grid element count would overflow");
    }

    // Allocate N grids, each nx*ny*nz.
    grids.reserve(ngrid);
    std::vector<std::vector<float>> buffers(ngrid, std::vector<float>(voxels));

    // Volumetric block: voxels*ngrid values, orbital-fastest then z,y,x.
    for (std::size_t vx = 0; vx < voxels; ++vx) {
        for (int g = 0; g < ngrid; ++g) {
            double val;
            if (!(in >> val)) {
                throw FormatError("oeio: CUBE: truncated volumetric block");
            }
            buffers[g][vx] = static_cast<float>(val);
        }
    }

    for (int g = 0; g < ngrid; ++g) {
        OESystem::OEScalarGrid grid(
            ax.nvox[0], ax.nvox[1], ax.nvox[2],
            static_cast<float>(mid[0]), static_cast<float>(mid[1]),
            static_cast<float>(mid[2]), static_cast<float>(spacing_ang));
        grid.SetValues(buffers[g].data(), static_cast<unsigned int>(voxels));
        grids.push_back(grid);
    }

    return true;
}

}  // namespace builtin
}  // namespace oeio
