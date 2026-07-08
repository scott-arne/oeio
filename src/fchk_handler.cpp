/// \file fchk_handler.cpp
/// \brief FchkHandler: registration and factory for the FCHK format (read-only).

#include "oeio/fchk_handler.h"

#include "oeio/exceptions.h"
#include "oeio/format_registry.h"

namespace oeio {
namespace builtin {

FormatInfo FchkHandler::info() const {
    return {
        "FCHK",
        { ".fchk", ".fch" },
        "Gaussian formatted checkpoint (molecule + QM scalars; read-only)",
        true,   // supports_read
        false,  // supports_write
        false,  // supports_threaded_read
        false   // supports_threaded_write
    };
}

std::unique_ptr<MolSource> FchkHandler::make_reader(
    const std::string& path, const std::any&) const {
    return std::make_unique<FchkMolSource>(path);
}

std::unique_ptr<MolSink> FchkHandler::make_writer(
    const std::string&, const std::any&) const {
    throw FormatError("oeio: FCHK is read-only; writing is not supported");
}

}  // namespace builtin

/// Force-link function so the FCHK handler's static registration survives static
/// linking (mirrors oeio_force_link_cube_handler).
void oeio_force_link_fchk_handler() {
    // Exists solely to be referenced from the SWIG module and tests so the
    // linker does not drop fchk_handler.o (and its OEIO_REGISTER_FORMAT static
    // initializer) from the static library.
}

}  // namespace oeio

OEIO_REGISTER_FORMAT(oeio::builtin::FchkHandler)
