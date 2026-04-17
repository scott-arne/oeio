/// \file mol_range.cpp
/// \brief Implementation of MolRange and its Iterator.

#include "oeio/mol_range.h"

namespace oeio {

// ============================================================================
// MolRange
// ============================================================================

MolRange::MolRange(std::unique_ptr<MolSource> source)
    : source_(std::move(source)) {}

MolRange::MolRange(MolRange&&) noexcept = default;
MolRange& MolRange::operator=(MolRange&&) noexcept = default;

MolRange::Iterator MolRange::begin() {
    return Iterator(source_.get());
}

MolRange::Sentinel MolRange::end() const {
    return Sentinel{};
}

std::unique_ptr<MolSource> MolRange::release_source() {
    return std::move(source_);
}

// ============================================================================
// MolRange::Iterator
// ============================================================================

MolRange::Iterator::Iterator() = default;

MolRange::Iterator::Iterator(MolSource* source)
    : source_(source), done_(false)
{
    ++(*this);  // Prime the first molecule.
}

MolRange::Iterator::reference MolRange::Iterator::operator*() {
    return mol_;
}

MolRange::Iterator::pointer MolRange::Iterator::operator->() {
    return &mol_;
}

MolRange::Iterator& MolRange::Iterator::operator++() {
    done_ = !source_ || !source_->next(mol_);
    return *this;
}

bool MolRange::Iterator::operator==(const Sentinel&) const {
    return done_;
}

bool MolRange::Iterator::operator!=(const Sentinel& s) const {
    return !(*this == s);
}

}  // namespace oeio
