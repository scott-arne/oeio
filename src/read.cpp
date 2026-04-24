/// \file read.cpp
/// \brief Implementation of Tier 1 read() overloads.

#include "oeio/read.h"
#include "oeio/exceptions.h"
#include "oeio/format_registry.h"

namespace oeio {

MolRange read(const std::string& path) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto source = handler->make_reader(path, std::any{});
    if (!source) {
        throw FileError(
            "oeio: failed to create reader for '" + path + "'");
    }
    return MolRange(std::move(source));
}

MolRange read(const std::string& path, std::any config) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto source = handler->make_reader(path, config);
    if (!source) {
        throw FileError(
            "oeio: failed to create reader for '" + path + "'");
    }
    return MolRange(std::move(source));
}

MolRange read(OEPlatform::oeifstream& stream, const std::string& format_hint) {
    if (format_hint.empty()) {
        throw FormatError(
            "oeio: stream-based read requires a format hint");
    }

    // Construct the extension with a leading dot if not already present.
    std::string ext = format_hint;
    if (ext.front() != '.') {
        ext = "." + ext;
    }

    auto* handler = FormatRegistry::instance().lookup_ext(ext);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized format hint '" + format_hint + "'");
    }

    auto source = handler->make_reader(stream, std::any{});
    if (!source) {
        throw FileError(
            "oeio: failed to create stream reader for format '" + format_hint + "'");
    }
    return MolRange(std::move(source));
}

}  // namespace oeio
