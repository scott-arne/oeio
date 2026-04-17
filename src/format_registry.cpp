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
        // Count dots to classify simple vs compound
        int dot_count = 0;
        for (char c : ext) {
            if (c == '.') ++dot_count;
        }

        if (dot_count > 1) {
            // Compound extension (e.g., ".sdf.gz") → sorted vector
            bool found = false;
            for (auto& entry : compound_ext_map_) {
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
                compound_ext_map_.push_back({ext, index});
            }
        } else {
            // Simple extension (e.g., ".sdf") → hash map
            auto it = ext_hash_.find(ext);
            if (it != ext_hash_.end()) {
                OESystem::OEThrow.Warning(
                    "oeio: extension '%s' re-registered (was %s, now %s)",
                    ext.c_str(),
                    handlers_[it->second]->info().name.c_str(),
                    info.name.c_str());
                it->second = index;
            } else {
                ext_hash_.emplace(ext, index);
            }
        }
    }

    // Sort compound extensions by length descending (longest match wins)
    std::sort(compound_ext_map_.begin(), compound_ext_map_.end(),
              [](const ExtEntry& a, const ExtEntry& b) {
                  return a.extension.size() > b.extension.size();
              });

    // Invalidate cache
    cached_ext_.clear();
    cached_handler_ = nullptr;
}

const FormatHandler* FormatRegistry::lookup(const std::string& path) const {
    // Uses unique_lock (not shared_lock) because cache writes are not
    // thread-safe under concurrent shared access. Lookup latency is
    // nanoseconds, so exclusive locking has negligible contention impact.
    std::unique_lock lock(mutex_);

    // Fast path: single-entry cache
    if (!cached_ext_.empty() &&
        path.size() >= cached_ext_.size() &&
        path.compare(path.size() - cached_ext_.size(),
                     cached_ext_.size(), cached_ext_) == 0) {
        return cached_handler_;
    }

    // Tier 1: compound extensions (longest match wins)
    for (const auto& entry : compound_ext_map_) {
        if (path.size() >= entry.extension.size() &&
            path.compare(path.size() - entry.extension.size(),
                         entry.extension.size(), entry.extension) == 0) {
            cached_ext_ = entry.extension;
            cached_handler_ = handlers_[entry.handler_index].get();
            return cached_handler_;
        }
    }

    // Tier 2: simple extension via hash map
    auto dot_pos = path.rfind('.');
    if (dot_pos != std::string::npos) {
        std::string ext = path.substr(dot_pos);
        auto it = ext_hash_.find(ext);
        if (it != ext_hash_.end()) {
            cached_ext_ = ext;
            cached_handler_ = handlers_[it->second].get();
            return cached_handler_;
        }
    }

    return nullptr;
}

const FormatHandler* FormatRegistry::lookup_ext(const std::string& ext) const {
    std::unique_lock lock(mutex_);

    // Fast path: single-entry cache
    if (ext == cached_ext_) {
        return cached_handler_;
    }

    // Try hash map first (handles simple extensions)
    auto it = ext_hash_.find(ext);
    if (it != ext_hash_.end()) {
        cached_ext_ = ext;
        cached_handler_ = handlers_[it->second].get();
        return cached_handler_;
    }

    // Fall back to compound vector for multi-dot extensions
    for (const auto& entry : compound_ext_map_) {
        if (entry.extension == ext) {
            cached_ext_ = ext;
            cached_handler_ = handlers_[entry.handler_index].get();
            return cached_handler_;
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
