#pragma once

/// \file write.h
/// \brief Functions and classes for writing molecules to files.

#include <any>
#include <memory>
#include <string>

#include <oechem.h>

#include "oeio/format_handler.h"
#include "oeio/format_registry.h"

namespace oeio {

/// \brief RAII wrapper for writing molecules to a file.
///
/// Writer manages the lifetime of a MolSink and provides a convenient
/// interface for writing molecules one at a time.
class Writer {
public:
    /// \brief Construct a Writer from a MolSink.
    ///
    /// \param sink Unique pointer to the MolSink (ownership transferred).
    explicit Writer(std::unique_ptr<MolSink> sink);

    /// \brief Destructor automatically closes the sink.
    ~Writer();

    // Move constructor and assignment.
    Writer(Writer&&) noexcept;
    Writer& operator=(Writer&&) noexcept;

    // Delete copy constructor and assignment.
    Writer(const Writer&) = delete;
    Writer& operator=(const Writer&) = delete;

    /// \brief Write a molecule to the file.
    ///
    /// \param mol The molecule to write.
    /// \returns true if the molecule was successfully written, false on error.
    bool add(const OEChem::OEMolBase& mol);

    /// \brief Close the writer and flush buffered data.
    void close();

private:
    std::unique_ptr<MolSink> sink_;
};

// ============================================================================
// Tier 1: Dynamic dispatch via FormatRegistry
// ============================================================================

/// \brief Create a Writer for a file path with auto-detected format.
///
/// \param path The file path to write to.
/// \returns A Writer for writing molecules to the file.
Writer write(const std::string& path);

/// \brief Create a Writer for a file path with format-specific configuration.
///
/// \param path The file path to write to.
/// \param config Format-specific configuration (type-erased via std::any).
/// \returns A Writer for writing molecules to the file.
Writer write(const std::string& path, std::any config);

// ============================================================================
// Tier 2: Compile-time dispatch for type-safe configuration
// ============================================================================

/// \brief Create a Writer with compile-time format selection.
///
/// This template function provides type-safe configuration by resolving the
/// format handler at compile time.
///
/// \tparam Format The format handler type (must have a static handler() method).
/// \param path The file path to write to.
/// \param config Format-specific configuration (defaults to default-constructed).
/// \returns A Writer for writing molecules to the file.
template <typename Format>
Writer write(const std::string& path,
             const typename Format::WriterConfig& config = {}) {
    auto sink = Format::handler()->make_writer(path, std::any(config));
    return Writer(std::move(sink));
}

}  // namespace oeio
