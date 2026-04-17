#pragma once

/// \file pipeline.h
/// \brief Pipeline operators for functional-style molecule processing.

#include <functional>

#include <oechem.h>

#include "oeio/mol_range.h"
#include "oeio/write.h"

namespace oeio {

// ============================================================================
// Pipe operator: MolRange | Writer
// ============================================================================

/// \brief Pipe a MolRange into a Writer (rvalue Writer).
///
/// Reads all molecules from the range and writes them to the writer.
///
/// \param range The molecule range to read from (rvalue).
/// \param writer The writer to write to (rvalue).
OEIO_FLATTEN void operator|(MolRange&& range, Writer&& writer);

/// \brief Pipe a MolRange into a Writer (lvalue Writer).
///
/// Reads all molecules from the range and writes them to the writer.
///
/// \param range The molecule range to read from (rvalue).
/// \param writer The writer to write to (lvalue reference).
OEIO_FLATTEN void operator|(MolRange&& range, Writer& writer);

// ============================================================================
// Functional pipeline operations
// ============================================================================

/// \brief Filter molecules by a predicate.
///
/// Returns a new MolRange that yields only molecules for which the predicate
/// returns true.
///
/// \param range The input molecule range (rvalue).
/// \param pred The predicate function (returns true to keep the molecule).
/// \returns A new MolRange with filtered molecules.
OEIO_HOT MolRange filter(MolRange&& range, std::function<bool(const OEChem::OEMolBase&)> pred);

/// \brief Transform molecules by applying a function.
///
/// Returns a new MolRange that applies the transformation function to each
/// molecule as it is read.
///
/// \param range The input molecule range (rvalue).
/// \param fn The transformation function (modifies the molecule in-place).
/// \returns A new MolRange with transformed molecules.
OEIO_HOT MolRange transform(MolRange&& range, std::function<void(OEChem::OEGraphMol&)> fn);

}  // namespace oeio
