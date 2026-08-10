/// \file oechem_handler.cpp
/// \brief Built-in OEChem format handler for native OpenEye file formats.

#include "oeio/exceptions.h"
#include "oeio/format_handler.h"
#include "oeio/format_registry.h"
#include "oeio/oechem_config.h"
#include "oeio/read_status.h"

#include <oechem.h>
#include <oesystem.h>

#include <cstddef>
#include <string>

namespace oeio {
namespace builtin {

namespace {

/// Reads one record using the OEReadMolecule overload that matches the
/// molecule's dynamic type.
///
/// OEReadMolecule's overloads resolve by static type, and the OEMolBase&
/// overload flattens conformers. Dispatching on the dynamic type is what lets an
/// OEMol (OEMCMolBase) assemble a multi-conformer record and an OEQMol
/// (OEQMolBase) read as a query.
bool read_dispatched(OEChem::oemolistream& ifs, OEChem::OEMolBase& mol) {
    if (auto* mc = dynamic_cast<OEChem::OEMCMolBase*>(&mol)) {
        return OEChem::OEReadMolecule(ifs, *mc);
    }
    if (auto* q = dynamic_cast<OEChem::OEQMolBase*>(&mol)) {
        return OEChem::OEReadMolecule(ifs, *q);
    }
    return OEChem::OEReadMolecule(ifs, mol);
}

}  // namespace

// ============================================================================
// OEChemMolSource — reads molecules via oemolistream
// ============================================================================

/// \brief MolSource implementation that reads molecules using OEChem's oemolistream.
class OEChemMolSource : public MolSource {
public:
    /// \brief Construct a reader from a file path.
    ///
    /// \param path The file path to read from.
    /// \param cfg Reader configuration.
    OEChemMolSource(const std::string& path, const oechem::ReaderConfig& cfg) {
        if (!ifs_.open(path)) {
            throw FileError("oeio: unable to open '" + path + "' for reading");
        }
        if (cfg.format != 0) {
            ifs_.SetFormat(cfg.format);
        }
        if (cfg.iflavor != 0) {
            unsigned int fmt = cfg.iflavor_format ? cfg.iflavor_format
                                                  : ifs_.GetFormat();
            ifs_.SetFlavor(fmt, cfg.iflavor);
        }
    }

    bool next(OEChem::OEGraphMol& mol) override {
        mol.Clear();
        return OEChem::OEReadMolecule(ifs_, mol);
    }

    bool next(OEChem::OEMolBase& mol) override {
        mol.Clear();
        return read_dispatched(ifs_, mol);
    }

    bool can_resynchronize() const override {
        // OEChem provides no skip-to-next-record primitive. Claiming resynchronized=true
        // would send callers into a retry loop on a stuck stream (the ReadResult contract
        // says resynchronized=true means positioned at the next boundary). A handler that
        // owns its tokenizer (e.g., oemaestro when implemented) can implement true recovery;
        // OEChemMolSource cannot. Return false to prevent infinite loops.
        return false;
    }

    oeio::ReadResult try_next(OEChem::OEMolBase& mol) override {
        // Terminal state: once a non-EOF failure occurs, all subsequent calls return EndOfStream.
        if (failed_) {
            mol.Clear();
            return oeio::read_end();
        }

        mol.Clear();

        // Same dispatch as next(), and for the same reason. It matters more here:
        // oeviz classifies an imported record as a trajectory by NumConfs(), so
        // reading through the OEMolBase& overload would flatten a multi-conformer
        // record and silently reclassify it as a single static structure.
        const bool read = read_dispatched(ifs_, mol);
        if (read) {
            ++records_read_;
            return oeio::read_ok();
        }

        // This is the distinction the boolean API cannot make: the stream is only at
        // its end if oemolistream says so. Otherwise the record failed.
        //
        // IMPORTANT: Whether a malformed record surfaces as RecordError depends on the
        // underlying reader's strictness. OEChem's SDF reader is very lenient — it
        // returns true and creates a partial molecule even for invalid data (emitting
        // only a stderr warning). Binary formats like OEB are stricter. This eof()
        // check is correct for readers that DO fail; for lenient readers, RecordError
        // is simply unreachable for most corruption.
        if (ifs_.eof()) {
            mol.Clear();
            return oeio::read_end();
        }

        // Non-EOF failure without resynchronization support: transition to terminal state.
        // A failed read without resynchronization support ends the stream; subsequent
        // calls return EndOfStream to prevent retry loops.
        failed_ = true;

        mol.Clear();
        std::string message = "record " + std::to_string(records_read_ + 1) + " of this " +
                              std::string(OEChem::OEGetFormatString(ifs_.GetFormat())) +
                              " stream could not be parsed";
        return oeio::read_error(std::move(message), false);
    }

private:
    mutable OEChem::oemolistream ifs_;
    std::size_t records_read_ = 0;
    bool failed_ = false;  // Terminal state: true after first non-EOF failure
};

// ============================================================================
// OEChemMolSink — writes molecules via oemolostream
// ============================================================================

/// \brief MolSink implementation that writes molecules using OEChem's oemolostream.
class OEChemMolSink : public MolSink {
public:
    /// \brief Construct a writer from a file path.
    ///
    /// \param path The file path to write to.
    /// \param cfg Writer configuration.
    OEChemMolSink(const std::string& path, const oechem::WriterConfig& cfg) {
        if (!ofs_.open(path)) {
            throw FileError("oeio: unable to open '" + path + "' for writing");
        }
        if (cfg.format != 0) {
            ofs_.SetFormat(cfg.format);
        }
        if (cfg.oflavor != 0) {
            unsigned int fmt = cfg.oflavor_format ? cfg.oflavor_format
                                                  : ofs_.GetFormat();
            ofs_.SetFlavor(fmt, cfg.oflavor);
        }
    }

    bool write(const OEChem::OEMolBase& mol) override {
        // OEWriteMolecule overloads resolve by static type. The OEMCMolBase&
        // overload writes every conformer as its own record using each
        // conformer's own title/data, whereas the OEMolBase& overload writes a
        // single record using the molecule-level title/data. Route genuinely
        // multi-conformer molecules through the OEMCMolBase& overload so their
        // conformers survive; keep single-conformer molecules on the
        // OEMolBase& overload so molecule-level title/SD data (as set by
        // transforms) is what gets written. OEWriteConstMolecule preserves
        // const-correctness (no const_cast).
        //
        // Tradeoff: for a genuinely multi-conformer molecule, per-record output
        // reflects each conformer's own title/data, so a molecule-level-only
        // title change may not surface — the necessary cost of not collapsing
        // conformers. No conformer data is lost.
        if (auto* mc = dynamic_cast<const OEChem::OEMCMolBase*>(&mol)) {
            if (mc->NumConfs() > 1) {
                return OEChem::OEWriteConstMolecule(ofs_, *mc) != 0;
            }
        }
        if (auto* q = dynamic_cast<const OEChem::OEQMolBase*>(&mol)) {
            return OEChem::OEWriteConstMolecule(ofs_, *q) != 0;
        }
        return OEChem::OEWriteConstMolecule(ofs_, mol) != 0;
    }

    void close() override {
        ofs_.close();
    }

private:
    OEChem::oemolostream ofs_;
};

// ============================================================================
// OEChemHandler — FormatHandler for OEChem native formats
// ============================================================================

/// \brief FormatHandler implementation for all OEChem-supported file formats.
class OEChemHandler : public FormatHandler {
public:
    FormatInfo info() const override {
        return {
            "OEChem",
            {
                // Compound (gzipped) extensions first — registry sorts by length,
                // but having them pre-sorted is cleaner for readability.
                ".oeb.gz", ".sdf.gz", ".mol2.gz", ".pdb.gz", ".smi.gz",
                ".csv.gz", ".xyz.gz", ".ent.gz",
                // Simple extensions.
                ".sdf", ".mol", ".mol2", ".pdb", ".ent",
                ".oeb", ".oez",
                ".smi", ".ism", ".can",
                ".csv", ".xyz",
                ".fasta", ".cif", ".mmcif",
                ".mopac", ".cdx"
            },
            "OpenEye OEChem native formats",
            true,   // supports_read
            true,   // supports_write
            false,  // supports_threaded_read  (oemolithread is separate)
            false   // supports_threaded_write (oemolothread is separate)
        };
    }

    std::unique_ptr<MolSource> make_reader(
        const std::string& path, const std::any& config) const override
    {
        oechem::ReaderConfig cfg;
        if (config.has_value()) {
            try {
                cfg = std::any_cast<oechem::ReaderConfig>(config);
            } catch (const std::bad_any_cast&) {
                OESystem::OEThrow.Warning(
                    "oeio: OEChem reader received unexpected config type; "
                    "using defaults");
            }
        }
        return std::make_unique<OEChemMolSource>(path, cfg);
    }

    std::unique_ptr<MolSink> make_writer(
        const std::string& path, const std::any& config) const override
    {
        oechem::WriterConfig cfg;
        if (config.has_value()) {
            try {
                cfg = std::any_cast<oechem::WriterConfig>(config);
            } catch (const std::bad_any_cast&) {
                OESystem::OEThrow.Warning(
                    "oeio: OEChem writer received unexpected config type; "
                    "using defaults");
            }
        }
        return std::make_unique<OEChemMolSink>(path, cfg);
    }
};

}  // namespace builtin

/// Force-link function to ensure OEChem handler gets registered
/// when linking statically.
void oeio_force_link_oechem_handler() {
    // This function exists solely to be referenced from tests
    // to prevent the linker from dropping oechem_handler.o
}

}  // namespace oeio

OEIO_REGISTER_FORMAT(oeio::builtin::OEChemHandler)
