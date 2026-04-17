#pragma once

/// \file mol_range.h
/// \brief Range-based interface for streaming molecules.

#include <functional>
#include <memory>

#include <oechem.h>

#include "oeio/format_handler.h"

namespace oeio {

/// \brief A range representing a stream of molecules.
///
/// MolRange provides a C++ range interface (begin/end) over a MolSource,
/// enabling range-based for loops and functional pipeline operations.
class MolRange {
public:
    /// \brief Construct a MolRange from a MolSource.
    ///
    /// \param source Unique pointer to the MolSource (ownership transferred).
    explicit MolRange(std::unique_ptr<MolSource> source);

    // Move constructor and assignment.
    MolRange(MolRange&&) noexcept;
    MolRange& operator=(MolRange&&) noexcept;

    // Delete copy constructor and assignment.
    MolRange(const MolRange&) = delete;
    MolRange& operator=(const MolRange&) = delete;

    /// \brief Sentinel type for range-based iteration.
    struct Sentinel {};

    /// \brief Iterator for streaming molecules from a MolSource.
    ///
    /// This is an input iterator that reads molecules one at a time.
    class Iterator {
    public:
        using value_type = OEChem::OEGraphMol;
        using reference = OEChem::OEGraphMol&;
        using pointer = OEChem::OEGraphMol*;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        Iterator();

        reference operator*();
        pointer operator->();
        Iterator& operator++();

        bool operator==(const Sentinel&) const;
        bool operator!=(const Sentinel& s) const;

    private:
        friend class MolRange;

        explicit Iterator(MolSource* source);

        MolSource* source_ = nullptr;
        OEChem::OEGraphMol mol_;
        bool done_ = true;
    };

    /// \brief Return an iterator to the first molecule.
    ///
    /// \returns Iterator positioned at the first molecule.
    Iterator begin();

    /// \brief Return a sentinel representing the end of the range.
    ///
    /// \returns Sentinel value for end-of-range comparison.
    Sentinel end() const;

private:
    std::unique_ptr<MolSource> release_source();

    friend MolRange filter(MolRange&&, std::function<bool(const OEChem::OEMolBase&)>);
    friend MolRange transform(MolRange&&, std::function<void(OEChem::OEGraphMol&)>);

    std::unique_ptr<MolSource> source_;
};

}  // namespace oeio
