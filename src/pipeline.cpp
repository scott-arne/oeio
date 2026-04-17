/// \file pipeline.cpp
/// \brief Implementation of pipeline operators and functional combinators.

#include "oeio/pipeline.h"

namespace oeio {
namespace detail {

// ============================================================================
// FilteredMolSource — skips molecules that don't match a predicate
// ============================================================================

/// \brief MolSource adapter that filters molecules by a predicate.
class FilteredMolSource : public MolSource {
public:
    /// \brief Construct a filtered source.
    ///
    /// \param upstream The upstream source to read from.
    /// \param pred Predicate that returns true for molecules to keep.
    FilteredMolSource(std::unique_ptr<MolSource> upstream,
                      std::function<bool(const OEChem::OEMolBase&)> pred)
        : upstream_(std::move(upstream)), pred_(std::move(pred)) {}

    bool next(OEChem::OEGraphMol& mol) override {
        while (upstream_->next(mol)) {
            if (pred_(mol))
                return true;
        }
        return false;
    }

private:
    std::unique_ptr<MolSource> upstream_;
    std::function<bool(const OEChem::OEMolBase&)> pred_;
};

// ============================================================================
// TransformedMolSource — applies a function to each molecule
// ============================================================================

/// \brief MolSource adapter that transforms molecules in-place.
class TransformedMolSource : public MolSource {
public:
    /// \brief Construct a transformed source.
    ///
    /// \param upstream The upstream source to read from.
    /// \param fn Transformation function applied to each molecule.
    TransformedMolSource(std::unique_ptr<MolSource> upstream,
                         std::function<void(OEChem::OEGraphMol&)> fn)
        : upstream_(std::move(upstream)), fn_(std::move(fn)) {}

    bool next(OEChem::OEGraphMol& mol) override {
        if (!upstream_->next(mol))
            return false;
        fn_(mol);
        return true;
    }

private:
    std::unique_ptr<MolSource> upstream_;
    std::function<void(OEChem::OEGraphMol&)> fn_;
};

}  // namespace detail

// ============================================================================
// Pipe operators
// ============================================================================

void operator|(MolRange&& range, Writer&& writer) {
    for (auto& mol : range)
        writer.add(mol);
    writer.close();
}

void operator|(MolRange&& range, Writer& writer) {
    for (auto& mol : range)
        writer.add(mol);
    // Do not close — caller owns the writer.
}

// ============================================================================
// Functional combinators
// ============================================================================

MolRange filter(MolRange&& range, std::function<bool(const OEChem::OEMolBase&)> pred) {
    auto source = range.release_source();
    return MolRange(std::make_unique<detail::FilteredMolSource>(
        std::move(source), std::move(pred)));
}

MolRange transform(MolRange&& range, std::function<void(OEChem::OEGraphMol&)> fn) {
    auto source = range.release_source();
    return MolRange(std::make_unique<detail::TransformedMolSource>(
        std::move(source), std::move(fn)));
}

}  // namespace oeio
