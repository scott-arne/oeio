/// \file cube_handler.cpp
/// \brief CubeHandler: registration and factory for the CUBE format.

#include "oeio/cube_handler.h"

#include "oeio/format_registry.h"

namespace oeio {
namespace builtin {

FormatInfo CubeHandler::info() const {
    return {
        "CUBE",
        { ".cube", ".cub" },
        "Gaussian CUBE volumetric format (molecule + scalar grids)",
        true,   // supports_read
        true,   // supports_write
        false,  // supports_threaded_read
        false   // supports_threaded_write
    };
}

std::unique_ptr<MolSource> CubeHandler::make_reader(
    const std::string& path, const std::any&) const {
    return std::make_unique<CubeMolSource>(path);
}

std::unique_ptr<MolSink> CubeHandler::make_writer(
    const std::string& path, const std::any&) const {
    return std::make_unique<CubeMolSink>(path);
}

}  // namespace builtin

/// Force-link function so the CUBE handler's static registration survives
/// static linking (mirrors oeio_force_link_oechem_handler).
void oeio_force_link_cube_handler() {
    // This function exists solely to be referenced from the SWIG module and
    // tests so the linker does not drop cube_handler.o (and its
    // OEIO_REGISTER_FORMAT static initializer) from the static library.
}

}  // namespace oeio

OEIO_REGISTER_FORMAT(oeio::builtin::CubeHandler)
