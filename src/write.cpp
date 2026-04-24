/// \file write.cpp
/// \brief Implementation of Writer class and Tier 1 write() overloads.

#include "oeio/write.h"
#include "oeio/exceptions.h"
#include "oeio/format_registry.h"

namespace oeio {

// ============================================================================
// Writer
// ============================================================================

Writer::Writer(std::unique_ptr<MolSink> sink)
    : sink_(std::move(sink)) {}

Writer::~Writer() {
    close();
}

Writer::Writer(Writer&&) noexcept = default;
Writer& Writer::operator=(Writer&&) noexcept = default;

bool Writer::append(const OEChem::OEMolBase& mol) {
    return sink_ ? sink_->write(mol) : false;
}

void Writer::close() {
    if (sink_) {
        sink_->close();
        sink_.reset();
    }
}

// ============================================================================
// Free functions
// ============================================================================

Writer write(const std::string& path) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto sink = handler->make_writer(path, std::any{});
    if (!sink) {
        throw FileError(
            "oeio: failed to create writer for '" + path + "'");
    }
    return Writer(std::move(sink));
}

Writer write(const std::string& path, std::any config) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto sink = handler->make_writer(path, config);
    if (!sink) {
        throw FileError(
            "oeio: failed to create writer for '" + path + "'");
    }
    return Writer(std::move(sink));
}

}  // namespace oeio
