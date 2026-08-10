#pragma once

/// \file format_handler.h
/// \brief Abstract base classes for molecular I/O format handling.

#include <any>
#include <memory>
#include <string>
#include <vector>

#include <oechem.h>
#include <oegrid.h>
#include <oeplatform.h>

#include <oeio/read_status.h>

#if defined(__GNUC__) || defined(__clang__)
#define OEIO_HOT [[gnu::hot]]
#define OEIO_FLATTEN __attribute__((flatten))
#else
#define OEIO_HOT
#define OEIO_FLATTEN
#endif

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
    OEIO_HOT virtual bool next(OEChem::OEGraphMol& mol) = 0;

    /// \brief Read the next molecule into an OEMolBase reference.
    ///
    /// Default implementation reads into a temp OEGraphMol and copies.
    /// Override for zero-copy when the handler can read into OEMolBase directly.
    ///
    /// \param mol The molecule base reference to populate.
    /// \returns true if a molecule was successfully read, false if end-of-stream or error.
    OEIO_HOT virtual bool next(OEChem::OEMolBase& mol) {
        OEChem::OEGraphMol temp;
        if (!next(temp)) return false;
        OEChem::OECopyMol(mol, temp);
        return true;
    }

    /// \brief Read the next molecule and up to grids.size() scalar grids
    /// (caller-owned, by-reference).
    ///
    /// \param mol The molecule to populate.
    /// \param grids Caller-owned grids to fill; min(grids.size(), N) are filled.
    /// \param num_grids If non-null, receives N (the record's grid count).
    /// \returns true if a record was read, false at end-of-stream.
    ///
    /// Default reads the molecule only and reports N=0, so non-grid handlers
    /// ignore grids for free.
    virtual bool next(OEChem::OEMolBase& mol,
                      const std::vector<OESystem::OEScalarGrid*>& grids,
                      int* num_grids = nullptr) {
        (void)grids;
        if (!next(mol)) return false;   // EOF: return false, do NOT touch *num_grids
        if (num_grids) *num_grids = 0;  // success: molecule only, N=0
        return true;
    }

    /// \brief Read the next molecule and all N scalar grids into an
    /// owned vector (resized to N).
    ///
    /// Used by the grid iterator view, which does not know N in advance.
    ///
    /// \param mol The molecule to populate.
    /// \param grids Resized to N and filled with the record's grids.
    /// \returns true if a record was read, false at end-of-stream.
    ///
    /// Default reads the molecule only and clears grids (N=0).
    virtual bool next(OEChem::OEMolBase& mol,
                      std::vector<OESystem::OEScalarGrid>& grids) {
        grids.clear();
        return next(mol);
    }

    /// \returns True if this format can skip a malformed record and continue at
    ///          the next record boundary. Record-delimited text formats can;
    ///          block-structured binary formats generally cannot.
    ///
    /// The default is false so that a format that has not been audited is
    /// reported as unrecoverable rather than silently assumed safe.
    virtual bool can_resynchronize() const { return false; }

    /// Reads one record, distinguishing end-of-stream from a record failure.
    ///
    /// The default implementation adapts the boolean next(): it cannot tell the
    /// two apart, so it reports EndOfStream. Formats that can tell must override
    /// this; oeio's own OEChem-backed source does.
    ///
    /// On RecordError, mol is left in a cleared state.
    virtual ReadResult try_next(OEChem::OEMolBase& mol) {
        return next(mol) ? read_ok() : read_end();
    }
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

    /// \brief Write a molecule plus grids.
    ///
    /// Default rejects a non-empty grid list (format has no grid support) and
    /// delegates to the molecule-only write when grids is empty.
    ///
    /// \param mol The molecule to write.
    /// \param grids The grids to write alongside the molecule.
    /// \returns true on success, false otherwise.
    virtual bool write(const OEChem::OEMolBase& mol,
                       const std::vector<const OESystem::OEScalarGrid*>& grids) {
        return grids.empty() ? write(mol) : false;
    }

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
