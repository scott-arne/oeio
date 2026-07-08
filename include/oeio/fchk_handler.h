#pragma once

/// \file fchk_handler.h
/// \brief Built-in Gaussian formatted-checkpoint (FCHK) format handler (read-only).

#include <any>
#include <memory>
#include <string>

#include <oechem.h>

#include "oeio/format_handler.h"

namespace oeio {
namespace builtin {

/// \brief MolSource that reads a Gaussian formatted checkpoint file into one
/// OEMol (atoms + Bohr->Angstrom coordinates) plus QM scalars as typed data.
///
/// An FCHK file holds exactly one molecule; the second next() returns false.
class FchkMolSource : public MolSource {
public:
    explicit FchkMolSource(const std::string& path);

    bool next(OEChem::OEGraphMol& mol) override;

private:
    /// Parse the single FCHK record from path_ into mol. Returns false once the
    /// record has already been consumed.
    bool read_record(OEChem::OEMolBase& mol);

    std::string path_;
    bool consumed_ = false;  // FCHK holds exactly one record
};

/// \brief FormatHandler for Gaussian FCHK files (read-only; no writer).
class FchkHandler : public FormatHandler {
public:
    FormatInfo info() const override;
    std::unique_ptr<MolSource> make_reader(const std::string& path,
                                           const std::any& config) const override;
    /// Writing FCHK is unsupported (a faithful file needs basis/MO data the
    /// reader does not retain); this raises FormatError.
    std::unique_ptr<MolSink> make_writer(const std::string& path,
                                         const std::any& config) const override;
};

}  // namespace builtin
}  // namespace oeio
