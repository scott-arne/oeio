/// \file read.cpp
/// \brief Implementation of Tier 1 read() overloads.

#include "oeio/read.h"
#include "oeio/format_registry.h"

#include <oesystem.h>

namespace oeio {

MolRange read(const std::string& path) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        OESystem::OEThrow.Fatal(
            "oeio: unrecognized file extension for '%s'", path.c_str());
    }
    auto source = handler->make_reader(path, std::any{});
    if (!source) {
        OESystem::OEThrow.Fatal(
            "oeio: failed to create reader for '%s'", path.c_str());
    }
    return MolRange(std::move(source));
}

MolRange read(const std::string& path, std::any config) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        OESystem::OEThrow.Fatal(
            "oeio: unrecognized file extension for '%s'", path.c_str());
    }
    auto source = handler->make_reader(path, config);
    if (!source) {
        OESystem::OEThrow.Fatal(
            "oeio: failed to create reader for '%s'", path.c_str());
    }
    return MolRange(std::move(source));
}

MolRange read(OEPlatform::oeifstream& stream, const std::string& format_hint) {
    if (format_hint.empty()) {
        OESystem::OEThrow.Fatal(
            "oeio: stream-based read requires a format hint");
    }

    // Construct the extension with a leading dot if not already present.
    std::string ext = format_hint;
    if (ext.front() != '.') {
        ext = "." + ext;
    }

    auto* handler = FormatRegistry::instance().lookup_ext(ext);
    if (!handler) {
        OESystem::OEThrow.Fatal(
            "oeio: unrecognized format hint '%s'", format_hint.c_str());
    }

    auto source = handler->make_reader(stream, std::any{});
    if (!source) {
        OESystem::OEThrow.Fatal(
            "oeio: failed to create stream reader for format '%s'",
            format_hint.c_str());
    }
    return MolRange(std::move(source));
}

}  // namespace oeio
