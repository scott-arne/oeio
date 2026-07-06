#pragma once

/// \file cube_handler.h
/// \brief Built-in Gaussian CUBE format handler.

#include <memory>
#include <string>
#include <vector>

#include <oechem.h>
#include <oegrid.h>

#include "oeio/format_handler.h"

namespace oeio {
namespace builtin {

/// \brief MolSource that reads Gaussian CUBE files (molecule + N scalar grids).
class CubeMolSource : public MolSource {
public:
    explicit CubeMolSource(const std::string& path);

    bool next(OEChem::OEGraphMol& mol) override;
    bool next(OEChem::OEMolBase& mol,
              const std::vector<OESystem::OEScalarGrid*>& grids,
              int* num_grids = nullptr) override;
    bool next(OEChem::OEMolBase& mol,
              std::vector<OESystem::OEScalarGrid>& grids) override;

private:
    /// Reads one CUBE record from the file into mol + owned grids.
    /// Returns false if the stream is exhausted (already consumed the single record).
    bool read_record(OEChem::OEMolBase& mol,
                     std::vector<OESystem::OEScalarGrid>& grids);

    std::string path_;
    bool consumed_ = false;  // CUBE holds exactly one record
};

/// \brief MolSink that writes Gaussian CUBE files (molecule + N scalar grids).
class CubeMolSink : public MolSink {
public:
    explicit CubeMolSink(const std::string& path);

    /// Writing a CUBE requires at least one grid; this raises FormatError.
    bool write(const OEChem::OEMolBase& mol) override;

    bool write(const OEChem::OEMolBase& mol,
               const std::vector<const OESystem::OEScalarGrid*>& grids) override;

    void close() override;

private:
    std::string path_;
};

}  // namespace builtin
}  // namespace oeio
