#pragma once

/// \file serialize.h
/// \brief In-memory single-molecule serialization to/from bytes and strings.
///
/// Thin wrappers over OpenEye's OEReadMolFromString / OEWriteMolToString for
/// OEChem-native formats. "bytes" is a std::string holding binary data (e.g.
/// OEB). These read/write exactly one record; multi-record buffers are out of
/// scope (see the streaming capability). CUBE/FCHK are not OEChem string
/// formats and raise FormatError.

#include <string>

#include <oechem.h>

namespace oeio {

/// \brief Resolve a format token ("sdf", ".oeb", "smi") to an OEFormat code.
/// \param fmt Extension/format token, with or without a leading dot.
/// \returns The OEChem::OEFormat code.
/// \raises FormatError if the token is not a recognized OEChem format.
unsigned resolve_format(const std::string& fmt);

/// \brief Serialize a molecule to a text string in the given format.
/// \param mol The molecule (single- or multi-conformer; dispatch is internal).
/// \param fmt An OEFormat code.
/// \param flavor Output flavor; 0 uses the per-format default.
/// \returns The serialized text.
/// \raises FormatError if the format is not writeable.
std::string mol_to_string(const OEChem::OEMolBase& mol, unsigned fmt,
                          unsigned flavor = 0);

/// \brief String-format overload of mol_to_string.
std::string mol_to_string(const OEChem::OEMolBase& mol, const std::string& fmt,
                          unsigned flavor = 0);

/// \brief Parse a molecule from a text string in the given format.
/// \param mol The molecule to populate (any OEMolBase subclass).
/// \param data The text to parse (first record only).
/// \param fmt An OEFormat code.
/// \param flavor Input flavor; 0 uses the per-format default.
/// \returns true if a molecule was read, false otherwise.
/// \raises FormatError if the format is not readable.
bool mol_from_string(OEChem::OEMolBase& mol, const std::string& data,
                     unsigned fmt, unsigned flavor = 0);

/// \brief String-format overload of mol_from_string.
bool mol_from_string(OEChem::OEMolBase& mol, const std::string& data,
                     const std::string& fmt, unsigned flavor = 0);

}  // namespace oeio
