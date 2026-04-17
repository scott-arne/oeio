#pragma once

/// \file format_registry.h
/// \brief Global registry for molecular file format handlers.

#include <memory>
#include <shared_mutex>
#include <string>
#include <vector>

#include "oeio/format_handler.h"

namespace oeio {

/// \brief Global registry for format handlers.
///
/// FormatRegistry is a thread-safe singleton that manages all registered
/// format handlers and provides lookup by file path or extension.
class FormatRegistry {
public:
    /// \brief Access the global FormatRegistry instance.
    ///
    /// \returns Reference to the singleton FormatRegistry.
    static FormatRegistry& instance();

    // Delete copy and move constructors and assignment operators.
    FormatRegistry(const FormatRegistry&) = delete;
    FormatRegistry& operator=(const FormatRegistry&) = delete;
    FormatRegistry(FormatRegistry&&) = delete;
    FormatRegistry& operator=(FormatRegistry&&) = delete;

    /// \brief Register a new format handler.
    ///
    /// \param handler Unique pointer to the handler to register (ownership transferred).
    void register_handler(std::unique_ptr<FormatHandler> handler);

    /// \brief Look up a format handler by file path.
    ///
    /// Matches against file extensions registered by each handler.
    ///
    /// \param path The file path to match.
    /// \returns Pointer to the matching FormatHandler, or nullptr if not found.
    const FormatHandler* lookup(const std::string& path) const;

    /// \brief Look up a format handler by file extension.
    ///
    /// \param ext The file extension to match (including leading dot, e.g., ".sdf").
    /// \returns Pointer to the matching FormatHandler, or nullptr if not found.
    const FormatHandler* lookup_ext(const std::string& ext) const;

    /// \brief Get metadata for all registered formats.
    ///
    /// \returns A vector of FormatInfo structures for all registered handlers.
    std::vector<FormatInfo> formats() const;

private:
    FormatRegistry() = default;

    mutable std::shared_mutex mutex_;
    std::vector<std::unique_ptr<FormatHandler>> handlers_;

    /// \brief Entry mapping an extension to a handler index.
    struct ExtEntry {
        std::string extension;
        std::size_t handler_index;
    };

    /// \brief Extension map sorted by extension length (descending) for longest-match lookup.
    std::vector<ExtEntry> ext_map_;
};

}  // namespace oeio

/// \brief Register a format handler at static initialization time.
///
/// Place this macro at file scope (outside the oeio namespace) to automatically
/// register a format handler when the shared library is loaded.
///
/// \param HandlerClass The fully-qualified name of the FormatHandler subclass
///        to instantiate and register.
#define OEIO_REGISTER_FORMAT_IMPL(HandlerClass, uid)              \
    namespace {                                                    \
    static const bool oeio_reg_##uid = [] {                        \
        ::oeio::FormatRegistry::instance().register_handler(       \
            std::make_unique<HandlerClass>());                     \
        return true;                                               \
    }();                                                           \
    }

#define OEIO_REGISTER_FORMAT(HandlerClass)                        \
    OEIO_REGISTER_FORMAT_IMPL(HandlerClass, __COUNTER__)
