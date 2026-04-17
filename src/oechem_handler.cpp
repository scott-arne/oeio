/// \file oechem_handler.cpp
/// \brief Built-in OEChem format handler for native OpenEye file formats.

#include "oeio/format_handler.h"
#include "oeio/format_registry.h"
#include "oeio/oechem_config.h"

#include <oechem.h>
#include <oesystem.h>

namespace oeio {
namespace builtin {

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
            OESystem::OEThrow.Fatal(
                "oeio: unable to open '%s' for reading", path.c_str());
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
        return OEChem::OEReadMolecule(ifs_, mol);
    }

private:
    OEChem::oemolistream ifs_;
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
            OESystem::OEThrow.Fatal(
                "oeio: unable to open '%s' for writing", path.c_str());
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
        // OEWriteMolecule takes a non-const reference; the OEChem API does not
        // modify the molecule but the signature is historical.  The const_cast
        // is safe here.
        auto& mutable_mol = const_cast<OEChem::OEMolBase&>(mol);
        return OEChem::OEWriteMolecule(ofs_, mutable_mol) != 0;
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
