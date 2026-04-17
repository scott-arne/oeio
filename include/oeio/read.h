#pragma once

/// \file read.h
/// \brief Functions for reading molecules from files and streams.

#include <any>
#include <string>

#include <oechem.h>
#include <oeplatform.h>

#include "oeio/format_registry.h"
#include "oeio/mol_range.h"

namespace oeio {

// ============================================================================
// Tier 1: Dynamic dispatch via FormatRegistry
// ============================================================================

/// \brief Read molecules from a file path with auto-detected format.
///
/// \param path The file path to read from.
/// \returns A MolRange for iterating over molecules in the file.
MolRange read(const std::string& path);

/// \brief Read molecules from a file path with format-specific configuration.
///
/// \param path The file path to read from.
/// \param config Format-specific configuration (type-erased via std::any).
/// \returns A MolRange for iterating over molecules in the file.
MolRange read(const std::string& path, std::any config);

/// \brief Read molecules from an input stream.
///
/// \param stream The input stream to read from.
/// \param format_hint Optional format hint (file extension or format name).
/// \returns A MolRange for iterating over molecules in the stream.
MolRange read(OEPlatform::oeifstream& stream, const std::string& format_hint = "");

// ============================================================================
// Tier 2: Compile-time dispatch for type-safe configuration
// ============================================================================

/// \brief Read molecules from a file path with compile-time format selection.
///
/// This template function provides type-safe configuration by resolving the
/// format handler at compile time.
///
/// \tparam Format The format handler type (must have a static handler() method).
/// \param path The file path to read from.
/// \param config Format-specific configuration (defaults to default-constructed).
/// \returns A MolRange for iterating over molecules in the file.
template <typename Format>
MolRange read(const std::string& path,
              const typename Format::ReaderConfig& config = {}) {
    auto source = Format::handler()->make_reader(path, std::any(config));
    return MolRange(std::move(source));
}

}  // namespace oeio
