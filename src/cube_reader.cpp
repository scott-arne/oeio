/// \file cube_reader.cpp
/// \brief CUBE parsing: file -> OEMol + N OEScalarGrid.

#include "oeio/cube_handler.h"

#include "oeio/cube_grid.h"
#include "oeio/exceptions.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace oeio {
namespace builtin {

namespace {
constexpr double AXIS_TOL = 1e-6;

/// Maximum allowed voxel count per dimension; prevents absurd allocations.
constexpr int MAX_VOXELS_PER_DIM = 8192;

/// Maximum total bytes across all grid buffers. Bounds the whole allocation by
/// memory footprint so a small malformed header cannot force a multi-gigabyte
/// reservation before any volumetric data is read. The derived element ceiling
/// stays well below UINT_MAX, so SetValues() lengths never truncate.
constexpr std::size_t MAX_TOTAL_BYTES = 512u * 1024u * 1024u;  // 512 MiB
constexpr std::size_t MAX_TOTAL_ELEMENTS = MAX_TOTAL_BYTES / sizeof(float);

/// Maximum orbital count for MO cubes; prevents huge grid arrays.
constexpr int MAX_ORBITAL_COUNT = 1024;

/// Highest supported atomic number (Oganesson, Z=118). Atom records outside
/// [1, MAX_ATOMIC_NUMBER] are rejected before reaching OpenEye.
constexpr int MAX_ATOMIC_NUMBER = 118;

/// Maximum atom count. Bounds the atom-parsing loop so a malformed header
/// advertising an absurd atom count fails cleanly rather than spinning.
constexpr long MAX_ATOM_COUNT = 100000000L;

/// Read the next line of a CUBE file into a stringstream for record-anchored
/// parsing. Fixed CUBE records occupy exactly one line each; reading a whole
/// line first prevents a short record from silently consuming tokens that
/// belong to the following line. Throws FormatError at end-of-file.
std::istringstream next_record_line(std::istream& in, const char* what) {
    std::string line;
    if (!std::getline(in, line)) {
        throw FormatError(std::string("oeio: CUBE: truncated file, expected ") + what);
    }
    return std::istringstream(line);
}

/// Parse exactly `count` doubles from a single record line, rejecting a line
/// that has too few fields or trailing non-whitespace content. The trailing
/// check keeps a malformed fixed record from being silently accepted.
std::vector<double> parse_doubles(std::istringstream& line, int count, const char* what) {
    std::vector<double> out;
    out.reserve(count);
    double v;
    for (int i = 0; i < count; ++i) {
        if (!(line >> v)) {
            throw FormatError(std::string("oeio: CUBE: malformed ") + what);
        }
        out.push_back(v);
    }
    std::string leftover;
    if (line >> leftover) {
        throw FormatError(std::string("oeio: CUBE: unexpected extra data on ") + what);
    }
    return out;
}

/// Absolute value of an untrusted signed count, guarding the one input whose
/// magnitude is not representable. std::labs(LONG_MIN) is undefined behavior
/// because -LONG_MIN overflows; reject that value (and anything that streamed
/// as it) as malformed before taking the magnitude.
long checked_abs(long value, const char* what) {
    if (value == std::numeric_limits<long>::min()) {
        throw FormatError(std::string("oeio: CUBE: ") + what + " out of range");
    }
    return std::labs(value);
}

/// Extract a whole integer token from a record line as a long. Reading an int
/// count with `operator>>` accepts an integer prefix and leaves any suffix
/// behind, so a token like "1.25" would parse as 1 and let ".25" splice into
/// the following double field. Pull the next whitespace-delimited token and
/// require it to be an optional sign followed by digits only, then range-check
/// the conversion, so a malformed count is rejected rather than silently split.
long parse_int_token(std::istream& line, const char* what) {
    std::string tok;
    if (!(line >> tok)) {
        throw FormatError(std::string("oeio: CUBE: missing ") + what);
    }
    std::size_t i = 0;
    if (i < tok.size() && (tok[i] == '+' || tok[i] == '-')) ++i;
    if (i == tok.size()) {  // sign with no digits
        throw FormatError(std::string("oeio: CUBE: malformed ") + what);
    }
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') {
            throw FormatError(std::string("oeio: CUBE: malformed ") + what);
        }
    }
    try {
        std::size_t consumed = 0;
        const long value = std::stol(tok, &consumed);
        if (consumed != tok.size()) {
            throw FormatError(std::string("oeio: CUBE: malformed ") + what);
        }
        return value;
    } catch (const std::out_of_range&) {
        throw FormatError(std::string("oeio: CUBE: ") + what + " out of range");
    } catch (const std::invalid_argument&) {
        throw FormatError(std::string("oeio: CUBE: malformed ") + what);
    }
}

/// Narrow an untrusted double to float, rejecting values that are non-finite or
/// whose magnitude exceeds the float range. A finite double such as 1e308 passes
/// an isfinite() check yet becomes +/-inf once cast to float; validating after
/// unit conversion at the single point where the value is narrowed prevents a
/// non-finite coordinate/geometry value from reaching OpenEye or OEScalarGrid.
float to_finite_float(double value, const char* what) {
    if (!std::isfinite(value) ||
        std::fabs(value) > static_cast<double>(std::numeric_limits<float>::max())) {
        throw FormatError(std::string("oeio: CUBE: ") + what +
                          " is non-finite or out of range");
    }
    return static_cast<float>(value);
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
    if (!std::getline(in, comment1) || !std::getline(in, comment2)) {
        throw FormatError("oeio: CUBE: truncated header (missing comment lines)");
    }

    // Line 3: atom count + origin (a single fixed record line).
    std::istringstream atom_count_line = next_record_line(in, "atom count line");
    const long natom_raw = parse_int_token(atom_count_line, "atom count");
    const bool mo_cube = natom_raw < 0;
    const long natom = checked_abs(natom_raw, "atom count");
    if (natom > MAX_ATOM_COUNT) {
        throw FormatError("oeio: CUBE: atom count exceeds limit (" +
                          std::to_string(MAX_ATOM_COUNT) + ")");
    }
    auto origin = parse_doubles(atom_count_line, 3, "origin");
    for (double o : origin) {
        if (!std::isfinite(o)) throw FormatError("oeio: CUBE: non-finite origin");
    }

    // Lines 4-6: voxel counts + axis vectors.
    cube::CubeAxes ax{};
    ax.origin[0] = origin[0]; ax.origin[1] = origin[1]; ax.origin[2] = origin[2];
    bool angstrom = false;
    for (int i = 0; i < 3; ++i) {
        // Each voxel count + axis vector is a single fixed record line.
        std::istringstream axis_line = next_record_line(in, "voxel count line");
        const long nv = parse_int_token(axis_line, "voxel count");
        if (nv < 0) angstrom = true;   // negative voxel count -> Angstrom
        const long nv_abs = checked_abs(nv, "voxel count");
        if (nv_abs < 1 || nv_abs > MAX_VOXELS_PER_DIM) {
            throw FormatError("oeio: CUBE: voxel count per dimension must be in [1, " +
                              std::to_string(MAX_VOXELS_PER_DIM) + "]");
        }
        ax.nvox[i] = static_cast<int>(nv_abs);
        auto row = parse_doubles(axis_line, 3, "axis vector");
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
        // Each atom is a single fixed record line; reading the whole line first
        // prevents a short atom record from consuming a volumetric value.
        std::istringstream atom_record = next_record_line(in, "atom line");
        auto atom_line = parse_doubles(atom_record, 5, "atom line");
        // The atomic number streams in as a double; validate it is a finite,
        // integral value in the supported element range before the narrowing
        // cast. A negative/NaN/Inf/out-of-range value would otherwise invoke
        // undefined float-to-unsigned conversion or hand OpenEye a bogus element.
        const double z = atom_line[0];
        if (!std::isfinite(z) || z < 1.0 || z > static_cast<double>(MAX_ATOMIC_NUMBER) ||
            z != std::floor(z)) {
            throw FormatError("oeio: CUBE: atomic number out of range [1, " +
                              std::to_string(MAX_ATOMIC_NUMBER) + "]");
        }
        // Convert to Angstrom and narrow to float, rejecting values that are
        // non-finite or overflow the float range (a finite double like 1e308
        // becomes +/-inf once cast). Validating the converted value guards the
        // coordinate that actually reaches the molecule.
        const float coords[3] = {
            to_finite_float(atom_line[2] * to_ang, "atom coordinate"),
            to_finite_float(atom_line[3] * to_ang, "atom coordinate"),
            to_finite_float(atom_line[4] * to_ang, "atom coordinate")
        };
        OEChem::OEAtomBase* atom = mol.NewAtom(static_cast<unsigned int>(z));
        if (!atom) {
            throw FormatError("oeio: CUBE: failed to create atom");
        }
        mol.SetCoords(atom, coords);
    }

    // MO cube: extra header line "M id1 id2 ... idM".
    int ngrid = 1;
    if (mo_cube) {
        const long m = parse_int_token(in, "orbital count");
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

    // Bound the total element count (voxels * ngrid) that a well-formed file may
    // request. The division form is overflow-safe: ngrid is already >= 1, so the
    // divisor is nonzero and no multiplication is evaluated. The comparison is
    // strict so a legitimate maximum grid (e.g. 512^3 == MAX_TOTAL_ELEMENTS) is
    // accepted. This ceiling also guarantees voxels <= MAX_TOTAL_ELEMENTS <
    // UINT_MAX, so the SetValues() length cast below cannot truncate.
    if (voxels > MAX_TOTAL_ELEMENTS / static_cast<std::size_t>(ngrid)) {
        throw FormatError("oeio: CUBE: total grid element count (" +
                          std::to_string(voxels) + " x " + std::to_string(ngrid) +
                          ") exceeds limit (" + std::to_string(MAX_TOTAL_ELEMENTS) +
                          " floats)");
    }

    // Read the volumetric block, growing per-grid buffers as values are actually
    // consumed. We deliberately do NOT pre-size to the header-declared voxel
    // count: a malformed or truncated file that advertises large dimensions but
    // supplies little or no data must fail fast, allocating only in proportion to
    // the bytes it genuinely contains rather than to the untrusted length prefix.
    // The MAX_TOTAL_ELEMENTS ceiling above still bounds a well-formed oversized
    // file. Values are interleaved orbital-fastest, then z, y, x.
    std::vector<std::vector<float>> buffers(ngrid);  // ngrid empty buffers
    for (std::size_t vx = 0; vx < voxels; ++vx) {
        for (int g = 0; g < ngrid; ++g) {
            double val;
            if (!(in >> val)) {
                throw FormatError("oeio: CUBE: truncated volumetric block");
            }
            // Physical density/orbital values are finite; reject NaN/Inf (which
            // stream as valid doubles) rather than storing corrupt grid data.
            buffers[g].push_back(to_finite_float(val, "volumetric value"));
        }
    }

    // Materialize each grid, releasing its source buffer immediately afterwards so
    // peak memory stays near a single grid's worth rather than all buffers plus
    // all grids at once.
    // Narrow the geometry to float once, rejecting values that overflow the
    // float range after unit conversion, so the grid never receives a non-finite
    // midpoint or spacing.
    const float mid_f[3] = {
        to_finite_float(mid[0], "grid midpoint"),
        to_finite_float(mid[1], "grid midpoint"),
        to_finite_float(mid[2], "grid midpoint")
    };
    const float spacing_f = to_finite_float(spacing_ang, "grid spacing");
    grids.reserve(ngrid);
    for (int g = 0; g < ngrid; ++g) {
        OESystem::OEScalarGrid grid(
            ax.nvox[0], ax.nvox[1], ax.nvox[2],
            mid_f[0], mid_f[1], mid_f[2], spacing_f);
        grid.SetValues(buffers[g].data(), static_cast<unsigned int>(voxels));
        grids.push_back(grid);
        std::vector<float>().swap(buffers[g]);  // release source buffer
    }

    return true;
}

}  // namespace builtin
}  // namespace oeio
