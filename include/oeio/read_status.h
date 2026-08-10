#pragma once

/// \file
/// Distinguishes a clean end of stream from a record that failed to parse.
///
/// The boolean next() cannot express the difference, so a corrupt record in the
/// middle of a file is indistinguishable from the end of the file and the read
/// silently truncates. ReadResult is the additive fix: existing callers keep the
/// boolean API, and callers that care get the reason and whether the reader
/// managed to skip past the bad record.

#include <string>

namespace oeio {

enum class ReadStatus {
    Ok,           ///< A molecule was read.
    EndOfStream,  ///< No more records; the stream ended cleanly.
    RecordError,  ///< One record failed. See resynchronized for whether to continue.
};

struct ReadResult {
    ReadStatus status = ReadStatus::EndOfStream;

    /// Human-readable reason, populated only for RecordError.
    std::string message;

    /// True when the reader positioned itself at the next record boundary, so
    /// the caller may keep reading. False means the stream must be abandoned.
    bool resynchronized = false;

    /// \returns True only for Ok, so `if (result)` reads as "got a molecule".
    explicit operator bool() const { return status == ReadStatus::Ok; }
};

inline ReadResult read_ok() { return ReadResult{ReadStatus::Ok, {}, false}; }
inline ReadResult read_end() { return ReadResult{ReadStatus::EndOfStream, {}, false}; }
inline ReadResult read_error(std::string message, bool resynchronized) {
    return ReadResult{ReadStatus::RecordError, std::move(message), resynchronized};
}

}  // namespace oeio
