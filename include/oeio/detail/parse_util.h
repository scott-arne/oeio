#pragma once

/// \file parse_util.h
/// \brief Shared parse/validation primitives for text-based molecular formats
/// (CUBE, FCHK). Diagnostics are format-neutral ("oeio: " + what); each caller
/// passes a context-carrying `what` (e.g. "CUBE atom coordinate",
/// "FCHK atomic number") so the message names the right format.

#include <cmath>
#include <cstdlib>   // std::labs
#include <istream>
#include <limits>
#include <stdexcept>
#include <string>

#include "oeio/exceptions.h"

namespace oeio {
namespace detail {

/// CODATA Bohr radius in Angstrom (1 Bohr = 0.52917721067 A). Single numeric
/// source of truth for every reader/writer that converts atomic-unit lengths.
constexpr double BOHR_TO_ANGSTROM = 0.52917721067;

/// Highest supported atomic number (Oganesson, Z=118). Atom records outside
/// [1, MAX_ATOMIC_NUMBER] are rejected.
constexpr int MAX_ATOMIC_NUMBER = 118;

/// Maximum atom count. Bounds atom loops so an absurd count fails cleanly.
constexpr long MAX_ATOM_COUNT = 100000000L;

/// Absolute value of an untrusted signed count, guarding the one input whose
/// magnitude is not representable. std::labs(LONG_MIN) is undefined behavior
/// because -LONG_MIN overflows; reject that value before taking the magnitude.
inline long checked_abs(long value, const char* what) {
    if (value == std::numeric_limits<long>::min()) {
        throw FormatError(std::string("oeio: ") + what + " out of range");
    }
    return std::labs(value);
}

/// Extract a whole integer token from a record line as a long. operator>> would
/// accept an integer prefix and leave a suffix behind (so "1.25" parses as 1 and
/// splices ".25" into the next field). Pull the next whitespace-delimited token,
/// require an optional sign then digits only, and range-check the conversion.
inline long parse_int_token(std::istream& line, const char* what) {
    std::string tok;
    if (!(line >> tok)) {
        throw FormatError(std::string("oeio: missing ") + what);
    }
    std::size_t i = 0;
    if (i < tok.size() && (tok[i] == '+' || tok[i] == '-')) ++i;
    if (i == tok.size()) {  // sign with no digits
        throw FormatError(std::string("oeio: malformed ") + what);
    }
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') {
            throw FormatError(std::string("oeio: malformed ") + what);
        }
    }
    try {
        std::size_t consumed = 0;
        const long value = std::stol(tok, &consumed);
        if (consumed != tok.size()) {
            throw FormatError(std::string("oeio: malformed ") + what);
        }
        return value;
    } catch (const std::out_of_range&) {
        throw FormatError(std::string("oeio: ") + what + " out of range");
    } catch (const std::invalid_argument&) {
        throw FormatError(std::string("oeio: malformed ") + what);
    }
}

/// Narrow an untrusted double to float, rejecting values that are non-finite or
/// whose magnitude exceeds the float range. A finite double such as 1e308 passes
/// isfinite() yet becomes +/-inf once cast; validate after any unit conversion.
inline float to_finite_float(double value, const char* what) {
    if (!std::isfinite(value) ||
        std::fabs(value) > static_cast<double>(std::numeric_limits<float>::max())) {
        throw FormatError(std::string("oeio: CUBE: ") + what +
                          " is non-finite or out of range");
    }
    return static_cast<float>(value);
}

}  // namespace detail
}  // namespace oeio
