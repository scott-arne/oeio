#pragma once

/// \file cube_handler.h
/// \brief Built-in Gaussian CUBE format handler.

#include <any>
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
    bool written_ = false;  // CUBE holds exactly one record; reject a second write
};

/// \brief FormatHandler for Gaussian CUBE files.
class CubeHandler : public FormatHandler {
public:
    FormatInfo info() const override;
    std::unique_ptr<MolSource> make_reader(const std::string& path,
                                           const std::any& config) const override;
    std::unique_ptr<MolSink> make_writer(const std::string& path,
                                         const std::any& config) const override;
};

}  // namespace builtin
}  // namespace oeio
