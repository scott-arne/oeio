/// \file write.cpp
/// \brief Implementation of Writer class and Tier 1 write() overloads.

#include "oeio/write.h"
#include "oeio/format_registry.h"

#include <oesystem.h>

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
        OESystem::OEThrow.Fatal(
            "oeio: unrecognized file extension for '%s'", path.c_str());
    }
    auto sink = handler->make_writer(path, std::any{});
    if (!sink) {
        OESystem::OEThrow.Fatal(
            "oeio: failed to create writer for '%s'", path.c_str());
    }
    return Writer(std::move(sink));
}

Writer write(const std::string& path, std::any config) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        OESystem::OEThrow.Fatal(
            "oeio: unrecognized file extension for '%s'", path.c_str());
    }
    auto sink = handler->make_writer(path, config);
    if (!sink) {
        OESystem::OEThrow.Fatal(
            "oeio: failed to create writer for '%s'", path.c_str());
    }
    return Writer(std::move(sink));
}

}  // namespace oeio
