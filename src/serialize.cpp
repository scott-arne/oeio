/// \file serialize.cpp
/// \brief Implementation of in-memory molecule serialization.

#include "oeio/serialize.h"

#include <vector>

#include "oeio/exceptions.h"

namespace oeio {

namespace {

/// Resolve flavor 0 to the per-format default write flavor.
unsigned out_flavor(unsigned fmt, unsigned flavor) {
    return flavor != 0 ? flavor : OEChem::OEGetDefaultOFlavor(fmt);
}

/// Resolve flavor 0 to the per-format default read flavor.
unsigned in_flavor(unsigned fmt, unsigned flavor) {
    return flavor != 0 ? flavor : OEChem::OEGetDefaultIFlavor(fmt);
}

/// Guarded write dispatch mirroring OEChemMolSink::write: only genuinely
/// multi-conformer molecules go through the OEMCMolBase overload (per-conformer
/// records); single-conformer molecules use the OEMolBase overload so
/// molecule-level title/SD data survives; queries use the OEQMolBase overload.
/// NOTE: the installed OpenEye 2025.2.3 flavor overload REQUIRES the gzip arg
/// (OEWriteMolToString(fmt, flavor, gzip, mol)); there is no (fmt, flavor, mol)
/// overload — so gzip is threaded through from the start (text callers pass
/// false).
std::string write_dispatch(const OEChem::OEMolBase& mol, unsigned fmt,
                           unsigned flavor, bool gzip) {
    const unsigned fl = out_flavor(fmt, flavor);
    if (auto* mc = dynamic_cast<const OEChem::OEMCMolBase*>(&mol)) {
        if (mc->NumConfs() > 1) {
            return OEChem::OEWriteMolToString(fmt, fl, gzip, *mc);
        }
    }
    // OEQMolBase and single-conformer molecules both route here; OEQMol is an
    // OEMolBase and serializes via this overload.
    return OEChem::OEWriteMolToString(fmt, fl, gzip, mol);
}

/// Guarded read dispatch mirroring OEChemMolSource::next(OEMolBase&).
bool read_dispatch(OEChem::OEMolBase& mol, const std::string& data,
                   unsigned fmt, unsigned flavor, bool gzip) {
    const unsigned fl = in_flavor(fmt, flavor);
    if (auto* mc = dynamic_cast<OEChem::OEMCMolBase*>(&mol)) {
        return OEChem::OEReadMolFromString(*mc, fmt, fl, gzip, data);
    }
    return OEChem::OEReadMolFromString(mol, fmt, fl, gzip, data);
}

}  // namespace

unsigned resolve_format(const std::string& fmt) {
    std::string ext = fmt;
    if (!ext.empty() && ext.front() == '.') ext.erase(0, 1);
    const unsigned code = OEChem::OEGetFileType(ext.c_str());
    if (code == OEChem::OEFormat::UNDEFINED) {
        throw FormatError("oeio: unknown or unsupported format '" + fmt +
                          "' (not an OEChem string/bytes format)");
    }
    return code;
}

std::string mol_to_string(const OEChem::OEMolBase& mol, unsigned fmt,
                          unsigned flavor) {
    if (!OEChem::OEIsWriteable(fmt)) {
        throw FormatError("oeio: format is not writeable as a string");
    }
    return write_dispatch(mol, fmt, flavor, /*gzip=*/false);
}

std::string mol_to_string(const OEChem::OEMolBase& mol, const std::string& fmt,
                          unsigned flavor) {
    return mol_to_string(mol, resolve_format(fmt), flavor);
}

bool mol_from_string(OEChem::OEMolBase& mol, const std::string& data,
                     unsigned fmt, unsigned flavor) {
    if (!OEChem::OEIsReadable(fmt)) {
        throw FormatError("oeio: format is not readable from a string");
    }
    return read_dispatch(mol, data, fmt, flavor, /*gzip=*/false);
}

bool mol_from_string(OEChem::OEMolBase& mol, const std::string& data,
                     const std::string& fmt, unsigned flavor) {
    return mol_from_string(mol, data, resolve_format(fmt), flavor);
}

std::string mol_to_bytes(const OEChem::OEMolBase& mol, unsigned fmt,
                         unsigned flavor, bool gzip) {
    if (!OEChem::OEIsWriteable(fmt)) {
        throw FormatError("oeio: format is not writeable as bytes");
    }
    return write_dispatch(mol, fmt, flavor, gzip);
}

std::string mol_to_bytes(const OEChem::OEMolBase& mol, const std::string& fmt,
                         unsigned flavor, bool gzip) {
    return mol_to_bytes(mol, resolve_format(fmt), flavor, gzip);
}

bool mol_from_bytes(OEChem::OEMolBase& mol, const std::string& data,
                    unsigned fmt, unsigned flavor, bool gzip) {
    if (!OEChem::OEIsReadable(fmt)) {
        throw FormatError("oeio: format is not readable from bytes");
    }
    return read_dispatch(mol, data, fmt, flavor, gzip);
}

bool mol_from_bytes(OEChem::OEMolBase& mol, const std::string& data,
                    const std::string& fmt, unsigned flavor, bool gzip) {
    return mol_from_bytes(mol, data, resolve_format(fmt), flavor, gzip);
}

}  // namespace oeio
