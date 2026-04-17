#pragma once

/// \file format_handler.h
/// \brief Abstract base classes for molecular I/O format handling.

#include <any>
#include <memory>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeplatform.h>

namespace oeio {

/// \brief Metadata describing a molecular file format.
struct FormatInfo {
    /// The human-readable name of the format (e.g., "Maestro", "SDF").
    std::string name;

    /// File extensions supported by this format (including leading dot, e.g., {".mae.gz", ".mae"}).
    std::vector<std::string> extensions;

    /// Human-readable description of the format.
    std::string description;

    /// Whether this format supports reading molecules.
    bool supports_read = true;

    /// Whether this format supports writing molecules.
    bool supports_write = true;

    /// Whether this format supports multi-threaded reading.
    bool supports_threaded_read = false;

    /// Whether this format supports multi-threaded writing.
    bool supports_threaded_write = false;
};

/// \brief Abstract base class for reading molecules from a source.
///
/// MolSource provides a streaming interface for reading molecules one at a time.
class MolSource {
public:
    virtual ~MolSource() = default;

    /// \brief Read the next molecule from the source.
    ///
    /// \param mol The molecule to populate with data from the source.
    /// \returns true if a molecule was successfully read, false if end-of-stream or error.
    virtual bool next(OEChem::OEGraphMol& mol) = 0;
};

/// \brief Abstract base class for writing molecules to a destination.
///
/// MolSink provides a streaming interface for writing molecules one at a time.
class MolSink {
public:
    virtual ~MolSink() = default;

    /// \brief Write a molecule to the sink.
    ///
    /// \param mol The molecule to write.
    /// \returns true if the molecule was successfully written, false on error.
    virtual bool write(const OEChem::OEMolBase& mol) = 0;

    /// \brief Close the sink and flush any buffered data.
    virtual void close() = 0;
};

/// \brief Abstract base class for format-specific I/O handlers.
///
/// FormatHandler implementations provide factory methods for creating readers
/// and writers for a specific molecular file format.
class FormatHandler {
public:
    virtual ~FormatHandler() = default;

    /// \brief Return metadata describing this format.
    ///
    /// \returns A FormatInfo structure with format metadata.
    virtual FormatInfo info() const = 0;

    /// \brief Create a reader for a file path.
    ///
    /// \param path The file path to read from.
    /// \param config Format-specific configuration (type-erased via std::any).
    /// \returns A unique_ptr to a MolSource for reading molecules.
    virtual std::unique_ptr<MolSource> make_reader(const std::string& path,
                                                    const std::any& config) const = 0;

    /// \brief Create a reader for an input stream.
    ///
    /// Default implementation throws an exception. Override if stream-based reading is supported.
    ///
    /// \param stream The input stream to read from.
    /// \param config Format-specific configuration (type-erased via std::any).
    /// \returns A unique_ptr to a MolSource for reading molecules.
    virtual std::unique_ptr<MolSource> make_reader(OEPlatform::oeifstream& stream,
                                                    const std::any& config) const;

    /// \brief Create a writer for a file path.
    ///
    /// \param path The file path to write to.
    /// \param config Format-specific configuration (type-erased via std::any).
    /// \returns A unique_ptr to a MolSink for writing molecules.
    virtual std::unique_ptr<MolSink> make_writer(const std::string& path,
                                                  const std::any& config) const = 0;
};

}  // namespace oeio
