#pragma once

/// \file oechem_config.h
/// \brief Configuration structures for OEChem-based I/O.

namespace oeio {
namespace oechem {

/// \brief Configuration for OEChem-based molecule readers.
struct ReaderConfig {
    /// The file format constant (e.g., OEFormat::SDF). 0 means auto-detect.
    unsigned int format = 0;

    /// Input flavor format specifier.
    unsigned int iflavor_format = 0;

    /// Input flavor flags.
    unsigned int iflavor = 0;

    /// Number of threads to use for reading (0 means single-threaded).
    unsigned int num_threads = 0;
};

/// \brief Configuration for OEChem-based molecule writers.
struct WriterConfig {
    /// The file format constant (e.g., OEFormat::SDF). 0 means auto-detect.
    unsigned int format = 0;

    /// Output flavor format specifier.
    unsigned int oflavor_format = 0;

    /// Output flavor flags.
    unsigned int oflavor = 0;

    /// Number of threads to use for writing (0 means single-threaded).
    unsigned int num_threads = 0;
};

}  // namespace oechem
}  // namespace oeio
