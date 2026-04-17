/// \file format_registry.cpp
/// \brief Implementation of FormatRegistry singleton and FormatHandler defaults.

#include "oeio/format_registry.h"

#include <algorithm>

#include <oesystem.h>

namespace oeio {

// ============================================================================
// FormatRegistry
// ============================================================================

FormatRegistry& FormatRegistry::instance() {
    static FormatRegistry registry;
    return registry;
}

void FormatRegistry::register_handler(std::unique_ptr<FormatHandler> handler) {
    std::unique_lock lock(mutex_);

    auto info = handler->info();
    handlers_.push_back(std::move(handler));
    std::size_t index = handlers_.size() - 1;

    for (const auto& ext : info.extensions) {
        bool found = false;
        for (auto& entry : ext_map_) {
            if (entry.extension == ext) {
                OESystem::OEThrow.Warning(
                    "oeio: extension '%s' re-registered (was %s, now %s)",
                    ext.c_str(),
                    handlers_[entry.handler_index]->info().name.c_str(),
                    info.name.c_str());
                entry.handler_index = index;
                found = true;
                break;
            }
        }
        if (!found) {
            ext_map_.push_back({ext, index});
        }
    }

    // Sort by extension length descending for longest-match-first lookup.
    std::sort(ext_map_.begin(), ext_map_.end(),
              [](const ExtEntry& a, const ExtEntry& b) {
                  return a.extension.size() > b.extension.size();
              });
}

const FormatHandler* FormatRegistry::lookup(const std::string& path) const {
    std::shared_lock lock(mutex_);

    for (const auto& entry : ext_map_) {
        if (path.size() >= entry.extension.size() &&
            path.compare(path.size() - entry.extension.size(),
                         entry.extension.size(), entry.extension) == 0) {
            return handlers_[entry.handler_index].get();
        }
    }
    return nullptr;
}

const FormatHandler* FormatRegistry::lookup_ext(const std::string& ext) const {
    std::shared_lock lock(mutex_);

    for (const auto& entry : ext_map_) {
        if (entry.extension == ext) {
            return handlers_[entry.handler_index].get();
        }
    }
    return nullptr;
}

std::vector<FormatInfo> FormatRegistry::formats() const {
    std::shared_lock lock(mutex_);

    std::vector<FormatInfo> result;
    result.reserve(handlers_.size());
    for (const auto& handler : handlers_) {
        result.push_back(handler->info());
    }
    return result;
}

// ============================================================================
// FormatHandler default stream-based reader
// ============================================================================

std::unique_ptr<MolSource> FormatHandler::make_reader(
    OEPlatform::oeifstream& /*stream*/, const std::any& /*config*/) const {
    OESystem::OEThrow.Fatal(
        "oeio: stream-based reading not supported for this format");
    return nullptr;
}

}  // namespace oeio
