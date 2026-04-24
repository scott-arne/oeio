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

    OEIO_HOT bool next(OEChem::OEGraphMol& mol) override {
        while (upstream_->next(mol)) {
            if (pred_(mol))
                return true;
        }
        return false;
    }

    OEIO_HOT bool next(OEChem::OEMolBase& mol) override {
        while (upstream_->next(mol)) {
            if (pred_(mol))
                return true;
        }
        return false;
    }

    /// Collapse a second predicate into this filter, returning a new
    /// FilteredMolSource that combines both predicates.
    std::unique_ptr<MolSource> collapse_with(
            std::function<bool(const OEChem::OEMolBase&)> other_pred) {
        auto combined = [p1 = std::move(pred_), p2 = std::move(other_pred)]
                        (const OEChem::OEMolBase& mol) {
            return p1(mol) && p2(mol);
        };
        return std::make_unique<FilteredMolSource>(
            std::move(upstream_), std::move(combined));
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

    OEIO_HOT bool next(OEChem::OEGraphMol& mol) override {
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

OEIO_FLATTEN void operator|(MolRange&& range, Writer&& writer) {
    for (auto& mol : range)
        writer.append(mol);
    writer.close();
}

OEIO_FLATTEN void operator|(MolRange&& range, Writer& writer) {
    for (auto& mol : range)
        writer.append(mol);
    // Do not close — caller owns the writer.
}

// ============================================================================
// Functional combinators
// ============================================================================

MolRange filter(MolRange&& range, std::function<bool(const OEChem::OEMolBase&)> pred) {
    auto source = range.release_source();

    // Filter chain collapsing: if the upstream is already a FilteredMolSource,
    // combine the predicates into a single filter to eliminate one virtual
    // dispatch level per molecule.
    auto* upstream_filter = dynamic_cast<detail::FilteredMolSource*>(source.get());
    if (upstream_filter) {
        auto combined_source = upstream_filter->collapse_with(std::move(pred));
        source.reset();
        return MolRange(std::move(combined_source));
    }

    return MolRange(std::make_unique<detail::FilteredMolSource>(
        std::move(source), std::move(pred)));
}

MolRange transform(MolRange&& range, std::function<void(OEChem::OEGraphMol&)> fn) {
    auto source = range.release_source();
    return MolRange(std::make_unique<detail::TransformedMolSource>(
        std::move(source), std::move(fn)));
}

}  // namespace oeio
