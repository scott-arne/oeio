#pragma once

/// \file mol_range.h
/// \brief Range-based interface for streaming molecules.

#include <functional>
#include <iterator>
#include <memory>
#include <type_traits>

#include <oechem.h>

#include "oeio/format_handler.h"

namespace oeio {

template <typename T>
class TypedMolRange;

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
        using value_type = OEChem::OEMol;
        using reference = OEChem::OEMol&;
        using pointer = OEChem::OEMol*;
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
        OEChem::OEMol mol_;
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

    /// \brief Read the next molecule into a caller-owned molecule (zero-copy).
    ///
    /// Forwards the caller's reference to the underlying source, so the
    /// dynamic type selects the correct read behavior (an OEMol reads
    /// multi-conformer, an OEGraphMol single-conformer).
    ///
    /// \param mol The molecule to populate.
    /// \returns true if a molecule was read, false if the source is exhausted or null.
    bool read_into(OEChem::OEMolBase& mol);

    /// \brief Return a typed view that yields molecules of type T.
    ///
    /// Moves the underlying MolSource into the returned view, leaving this
    /// MolRange empty (safe for the read("x").as<OEMol>() temporary idiom).
    /// T must be convertible to OEChem::OEMolBase& (e.g. OEMol, OEGraphMol,
    /// OEQMol).
    ///
    /// \tparam T The molecule type to read into.
    /// \returns A TypedMolRange<T> over the moved source.
    template <typename T>
    TypedMolRange<T> as();

private:
    std::unique_ptr<MolSource> release_source();

    friend MolRange filter(MolRange&&, std::function<bool(const OEChem::OEMolBase&)>);
    friend MolRange transform(MolRange&&, std::function<void(OEChem::OEMolBase&)>);

    std::unique_ptr<MolSource> source_;
};

/// \brief A typed view over a MolSource that yields molecules of type T.
///
/// T must be convertible to OEChem::OEMolBase& (verified via static_assert).
/// The iterator reuses a single T buffer, mirroring MolRange::Iterator.
template <typename T>
class TypedMolRange {
    static_assert(std::is_convertible<T&, OEChem::OEMolBase&>::value,
                  "TypedMolRange<T>: T must be convertible to OEChem::OEMolBase&");
    static_assert(std::is_default_constructible<T>::value,
                  "TypedMolRange<T>: T must be default-constructible");

public:
    explicit TypedMolRange(std::unique_ptr<MolSource> source)
        : source_(std::move(source)) {}

    TypedMolRange(TypedMolRange&&) noexcept = default;
    TypedMolRange& operator=(TypedMolRange&&) noexcept = default;
    TypedMolRange(const TypedMolRange&) = delete;
    TypedMolRange& operator=(const TypedMolRange&) = delete;

    struct Sentinel {};

    class Iterator {
    public:
        using value_type = T;
        using reference = T&;
        using pointer = T*;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        Iterator() = default;
        explicit Iterator(MolSource* source) : source_(source), done_(false) {
            ++(*this);
        }

        reference operator*() { return mol_; }
        pointer operator->() { return &mol_; }
        Iterator& operator++() {
            OEChem::OEMolBase& ref = mol_;  // rely on operator OEMolBase&()
            done_ = !source_ || !source_->next(ref);
            return *this;
        }
        bool operator==(const Sentinel&) const { return done_; }
        bool operator!=(const Sentinel& s) const { return !(*this == s); }

    private:
        MolSource* source_ = nullptr;
        T mol_;
        bool done_ = true;
    };

    Iterator begin() { return Iterator(source_.get()); }
    Sentinel end() const { return Sentinel{}; }

private:
    std::unique_ptr<MolSource> source_;
};

template <typename T>
TypedMolRange<T> MolRange::as() {
    return TypedMolRange<T>(release_source());
}

}  // namespace oeio
