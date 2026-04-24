#pragma once

/// \file exceptions.h
/// \brief Exception hierarchy for oeio.

#include <stdexcept>
#include <string>

namespace oeio {

/// \brief Base class for all oeio exceptions.
class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// \brief Raised when a file extension or format hint is not registered.
class FormatError : public Error {
public:
    using Error::Error;
};

/// \brief Raised when a file cannot be opened or a reader/writer cannot be created.
class FileError : public Error {
public:
    using Error::Error;
};

}  // namespace oeio
