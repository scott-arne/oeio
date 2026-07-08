/// \file fchk_reader.cpp
/// \brief FCHK parsing: file -> OEMol + typed QM scalars.

#include "oeio/fchk_handler.h"

#include "oeio/detail/parse_util.h"
#include "oeio/exceptions.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace oeio {
namespace builtin {

namespace {
using detail::BOHR_TO_ANGSTROM;
using detail::MAX_ATOM_COUNT;
using detail::MAX_ATOMIC_NUMBER;
using detail::parse_int_token;
using detail::to_finite_float;

/// Values-per-line stride for each FCHK array type code (formchk.txt write
/// formats: 6I12 / 5E16.8 / 5A12 / 72L1 / 9A8). Drives the line-based skip.
int values_per_line(char type) {
    switch (type) {
        case 'I': return 6;
        case 'R': return 5;
        case 'C': return 5;
        case 'L': return 72;
        case 'H': return 9;
        default: return 0;
    }
}

/// Right-trim trailing whitespace.
std::string rtrim(const std::string& s) {
    std::size_t end = s.find_last_not_of(" \t\r\n");
    return end == std::string::npos ? std::string() : s.substr(0, end + 1);
}

/// Trim leading and trailing whitespace. Used for the fixed-width method/basis
/// header fields, which are space-padded on BOTH sides in real files (e.g. the
/// example fixture's basis field is leading-space padded before "6-311+G(2d,p)").
std::string trim(const std::string& s) {
    std::size_t begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return std::string();
    std::size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

/// Upper bound on an array's declared N=. A real FCHK array is at most a few
/// million values; this ceiling (a) rejects absurd/hostile counts cleanly and
/// (b) keeps the payload-line ceiling arithmetic (n + per - 1) well clear of
/// signed-long overflow. It is far above any legitimate array (e.g. a 10k-basis
/// job's MO coefficients are ~1e8, still under this bound only if genuinely
/// present — and such a payload would be bounded by the file's own line count).
constexpr long MAX_ARRAY_COUNT = 1000000000L;  // 1e9

/// A parsed FCHK record header. is_array indicates an "N=" array declaration.
struct Header {
    std::string label;  // columns 1-40, right-trimmed
    char type = ' ';    // column 44 (A40,3X,A1)
    bool is_array = false;
    long count = 0;     // N= for arrays
};

/// Parse a record header line by fixed columns. Throws on an unrecognized type
/// code, a missing/blank label, or an out-of-range array count.
Header parse_header(const std::string& line) {
    if (line.size() < 44) {
        throw FormatError("oeio: FCHK: malformed record header (line too short)");
    }
    Header h;
    h.label = rtrim(line.substr(0, 40));
    if (h.label.empty()) {
        throw FormatError("oeio: FCHK: missing record label");
    }
    h.type = line[43];
    if (h.type != 'I' && h.type != 'R' && h.type != 'C' &&
        h.type != 'L' && h.type != 'H') {
        throw FormatError("oeio: FCHK: unrecognized type code in record '" +
                          h.label + "'");
    }
    std::size_t npos = line.find("N=", 44);
    if (npos != std::string::npos) {
        h.is_array = true;
        std::istringstream cnt(line.substr(npos + 2));
        h.count = parse_int_token(cnt, ("FCHK array count for '" + h.label + "'").c_str());
        // parse_int_token already range-checks against long overflow; bound the
        // magnitude here so the payload-line ceiling arithmetic cannot overflow
        // and an absurd count fails cleanly rather than misbehaving.
        if (h.count < 0 || h.count > MAX_ARRAY_COUNT) {
            throw FormatError("oeio: FCHK: array count out of range [0, " +
                              std::to_string(MAX_ARRAY_COUNT) + "] for '" +
                              h.label + "'");
        }
    }
    return h;
}

/// True if a line is a well-formed record header (used by the resync invariant).
bool looks_like_header(const std::string& line) {
    if (line.size() < 44) return false;
    if (rtrim(line.substr(0, 40)).empty()) return false;
    char t = line[43];
    return t == 'I' || t == 'R' || t == 'C' || t == 'L' || t == 'H';
}

/// Number of continuation lines an array of N values of the given type occupies.
long payload_lines(long n, char type) {
    int per = values_per_line(type);
    if (per <= 0 || n <= 0) return 0;
    return (n + per - 1) / per;  // ceil(n / per)
}

/// Read exactly `n` whitespace-delimited tokens from `nlines` lines of `lines`
/// starting at index `first`, erroring on too few tokens. The caller has already
/// bounds-checked that `first + nlines <= lines.size()`.
std::vector<std::string> read_tokens(const std::vector<std::string>& lines,
                                     std::size_t first, long n, long nlines,
                                     const std::string& label) {
    std::vector<std::string> out;
    out.reserve(static_cast<std::size_t>(n));
    for (long i = 0; i < nlines; ++i) {
        std::istringstream ss(lines[first + static_cast<std::size_t>(i)]);
        std::string tok;
        while (ss >> tok) out.push_back(tok);
    }
    if (static_cast<long>(out.size()) < n) {
        throw FormatError("oeio: FCHK: too few values for '" + label + "'");
    }
    out.resize(static_cast<std::size_t>(n));
    return out;
}

// Message convention (Task 1 option b): the shared parse_util helpers emit
// "oeio: <what>", and every caller passes a `what` that begins with the format
// tag, e.g. "FCHK atomic number". These FCHK wrappers therefore also emit
// neutral "oeio: ..." text (their `what` already carries "FCHK"), so no message
// double-tags the format and none ever says "CUBE" for an FCHK failure.

/// Parse a strict whole-integer token (reusing the shared parse_int_token rule).
long parse_strict_int(const std::string& tok, const char* what) {
    std::istringstream ss(tok);
    long v = parse_int_token(ss, what);
    std::string leftover;
    if (ss >> leftover) {
        throw FormatError(std::string("oeio: malformed ") + what);
    }
    return v;
}

/// Parse a finite double token.
double parse_finite_double(const std::string& tok, const char* what) {
    try {
        std::size_t consumed = 0;
        double v = std::stod(tok, &consumed);
        if (consumed != tok.size() || !std::isfinite(v)) {
            throw FormatError(std::string("oeio: malformed ") + what);
        }
        return v;
    } catch (const std::invalid_argument&) {
        throw FormatError(std::string("oeio: malformed ") + what);
    } catch (const std::out_of_range&) {
        throw FormatError(std::string("oeio: ") + what + " out of range");
    }
}
}  // namespace

FchkMolSource::FchkMolSource(const std::string& path) : path_(path) {}

bool FchkMolSource::next(OEChem::OEGraphMol& mol) {
    return read_record(mol);
}

bool FchkMolSource::read_record(OEChem::OEMolBase& mol) {
    mol.Clear();
    if (consumed_) return false;
    consumed_ = true;

    std::ifstream in(path_);
    if (!in) {
        throw FileError("oeio: unable to open '" + path_ + "' for reading");
    }

    // Read the whole file into a line vector. FCHK files are small (the real
    // fixture is ~230 KB); indexing lines is simpler and safer than interleaving
    // getline() with tellg()/seekg() to implement the resync peek.
    std::vector<std::string> lines;
    std::string raw;
    while (std::getline(in, raw)) {
        if (!raw.empty() && raw.back() == '\r') raw.pop_back();  // tolerate CRLF
        lines.push_back(raw);
    }
    if (lines.empty()) {
        throw FormatError("oeio: FCHK: empty file (missing title line)");
    }
    if (lines.size() < 2) {
        throw FormatError("oeio: FCHK: truncated header (missing calc line)");
    }

    // Header line 1: job title.
    mol.SetTitle(rtrim(lines[0]).c_str());

    // Header line 2: A10 calc type, A30 method, A30 basis.
    const std::string& calc_line = lines[1];
    auto col = [&](std::size_t a, std::size_t len) -> std::string {
        if (a >= calc_line.size()) return std::string();
        return trim(calc_line.substr(a, len));  // both-sides trim: fields are space-padded
    };
    const std::string method = col(10, 30);
    const std::string basis = col(40, 30);
    if (!method.empty()) mol.SetData("Method", method);
    if (!basis.empty()) mol.SetData("Basis", basis);

    // Record scan (order-independent). i walks the line vector; each record
    // advances i past its header and (for arrays) its payload lines.
    bool have_natom = false;
    long natom = 0;
    std::vector<long> atomic_numbers;
    std::vector<double> coords_bohr;

    std::size_t i = 2;
    while (i < lines.size()) {
        if (rtrim(lines[i]).empty()) {
            // A blank line is allowed only as trailing EOF padding: if anything
            // non-blank follows, a header was expected here -> malformed.
            for (std::size_t j = i + 1; j < lines.size(); ++j) {
                if (!rtrim(lines[j]).empty()) {
                    throw FormatError("oeio: FCHK: unexpected blank line where a "
                                      "record header was expected");
                }
            }
            break;  // only trailing blanks remain
        }
        Header h = parse_header(lines[i]);

        if (!h.is_array) {
            // Scalar: value is inline on the header line, after the type code.
            const std::string value =
                lines[i].size() > 44 ? rtrim(lines[i].substr(44)) : std::string();
            if (h.label == "Number of atoms") {
                natom = parse_strict_int(value, "FCHK atom count");
                have_natom = true;
            } else if (h.label == "Charge") {
                mol.SetData("Charge",
                    static_cast<int>(parse_strict_int(value, "FCHK charge")));
            } else if (h.label == "Multiplicity") {
                mol.SetData("Multiplicity",
                    static_cast<int>(parse_strict_int(value, "FCHK multiplicity")));
            } else if (h.label == "Number of electrons") {
                mol.SetData("Number of electrons",
                    static_cast<int>(parse_strict_int(value, "FCHK electron count")));
            } else if (h.label == "Total Energy") {
                mol.SetData("Total Energy",
                    parse_finite_double(value, "FCHK total energy"));
            }
            // Unrecognized scalar: ignore value.
            ++i;
            continue;
        }

        // Array: the payload occupies `plines` lines after the header.
        const long plines = payload_lines(h.count, h.type);
        const std::size_t first = i + 1;
        if (first + static_cast<std::size_t>(plines) > lines.size()) {
            throw FormatError("oeio: FCHK: truncated payload for '" + h.label + "'");
        }
        if (h.label == "Atomic numbers" && h.type == 'I') {
            auto toks = read_tokens(lines, first, h.count, plines, h.label);
            atomic_numbers.reserve(toks.size());
            for (const auto& t : toks) {
                atomic_numbers.push_back(parse_strict_int(t, "FCHK atomic number"));
            }
        } else if (h.label == "Current cartesian coordinates" && h.type == 'R') {
            auto toks = read_tokens(lines, first, h.count, plines, h.label);
            coords_bohr.reserve(toks.size());
            for (const auto& t : toks) {
                coords_bohr.push_back(parse_finite_double(t, "FCHK coordinate"));
            }
        }
        // else: skipped array — advance past its payload without reading.

        // Advance past header + payload.
        i = first + static_cast<std::size_t>(plines);

        // Resync invariant: the next non-blank line (if any) must be a header. A
        // wrong stride for a rare type would otherwise silently misread data.
        std::size_t k = i;
        while (k < lines.size() && rtrim(lines[k]).empty()) ++k;
        if (k < lines.size() && !looks_like_header(lines[k])) {
            throw FormatError("oeio: FCHK: lost record sync after '" + h.label + "'");
        }
    }

    // Validation.
    if (!have_natom) {
        throw FormatError("oeio: FCHK: missing 'Number of atoms'");
    }
    if (natom < 1 || natom > MAX_ATOM_COUNT) {
        throw FormatError("oeio: FCHK: atom count out of range [1, " +
                          std::to_string(MAX_ATOM_COUNT) + "]");
    }
    if (atomic_numbers.empty()) {
        throw FormatError("oeio: FCHK: missing 'Atomic numbers'");
    }
    if (static_cast<long>(atomic_numbers.size()) != natom) {
        throw FormatError("oeio: FCHK: atomic-number count does not match atom count");
    }
    if (coords_bohr.empty()) {
        throw FormatError("oeio: FCHK: missing 'Current cartesian coordinates'");
    }
    if (static_cast<long>(coords_bohr.size()) != 3 * natom) {
        throw FormatError("oeio: FCHK: coordinate count does not match 3 * atom count");
    }

    // Build the molecule.
    for (long a = 0; a < natom; ++a) {
        const long z = atomic_numbers[static_cast<std::size_t>(a)];
        if (z < 1 || z > MAX_ATOMIC_NUMBER) {
            throw FormatError("oeio: FCHK: atomic number out of range [1, " +
                              std::to_string(MAX_ATOMIC_NUMBER) + "]");
        }
        const float xyz[3] = {
            to_finite_float(coords_bohr[static_cast<std::size_t>(3 * a + 0)] * BOHR_TO_ANGSTROM,
                            "FCHK atom coordinate"),
            to_finite_float(coords_bohr[static_cast<std::size_t>(3 * a + 1)] * BOHR_TO_ANGSTROM,
                            "FCHK atom coordinate"),
            to_finite_float(coords_bohr[static_cast<std::size_t>(3 * a + 2)] * BOHR_TO_ANGSTROM,
                            "FCHK atom coordinate")
        };
        OEChem::OEAtomBase* atom = mol.NewAtom(static_cast<unsigned int>(z));
        if (!atom) {
            throw FormatError("oeio: FCHK: failed to create atom");
        }
        mol.SetCoords(atom, xyz);
    }
    return true;
}

}  // namespace builtin
}  // namespace oeio
