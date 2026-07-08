/// \file test_fchk_handler.cpp
/// \brief Tests for the Gaussian FCHK format handler.

#include <gtest/gtest.h>

#include <any>
#include <string>
#include <vector>

#include <oechem.h>
#include <oesystem.h>  // OESystem::OEIter for GetAtoms()

#include "oeio/exceptions.h"
#include "oeio/fchk_handler.h"
#include "oeio/format_registry.h"

namespace {

std::string data_path(const std::string& name) {
    return std::string(OEIO_TEST_DATA_DIR) + "/" + name;
}

TEST(FchkReader, MinFixtureReadsMoleculeAndScalars) {
    oeio::builtin::FchkMolSource src(data_path("fchk_min.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 2u);

    std::vector<unsigned int> zs;
    // Repo pattern: GetAtoms() returns OESystem::OEIter<const OEChem::OEAtomBase>
    // (see tests/cpp/test_cube_handler.cpp). There is no OEAtomBaseIter alias.
    for (OESystem::OEIter<const OEChem::OEAtomBase> it = mol.GetAtoms(); it; ++it) {
        zs.push_back(it->GetAtomicNum());
    }
    ASSERT_EQ(zs.size(), 2u);
    EXPECT_EQ(zs[0], 8u);
    EXPECT_EQ(zs[1], 1u);

    EXPECT_STREQ(mol.GetTitle(), "Min FCHK fixture");
    EXPECT_TRUE(mol.HasData("Total Energy"));
    EXPECT_NEAR(mol.GetDoubleData("Total Energy"), -75.0, 1e-9);
    EXPECT_EQ(mol.GetIntData("Charge"), -1);
    EXPECT_EQ(mol.GetIntData("Multiplicity"), 1);
    EXPECT_EQ(mol.GetIntData("Number of electrons"), 10);
    EXPECT_EQ(mol.GetStringData("Method"), "RHF");
    EXPECT_EQ(mol.GetStringData("Basis"), "STO-3G");

    // Second read: EOF.
    EXPECT_FALSE(src.next(mol));
}

TEST(FchkReader, RealExampleReadsHofAndSkipsLargeArrays) {
    // Real Gaussian 16 APFD/6-311+G(2d,p) single point on HOF. Exercises the
    // skip stride over large unconsumed arrays (Alpha MO coefficients N=3600,
    // Total SCF Density N=1830) and the empty-payload C array (Atom Types).
    oeio::builtin::FchkMolSource src(data_path("example.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 3u);

    std::vector<unsigned int> zs;
    for (OESystem::OEIter<const OEChem::OEAtomBase> it = mol.GetAtoms(); it; ++it) {
        zs.push_back(it->GetAtomicNum());
    }
    ASSERT_EQ(zs.size(), 3u);
    EXPECT_EQ(zs[0], 8u);  // O
    EXPECT_EQ(zs[1], 9u);  // F
    EXPECT_EQ(zs[2], 1u);  // H

    EXPECT_STREQ(mol.GetTitle(), "Example");
    EXPECT_NEAR(mol.GetDoubleData("Total Energy"), -175.4626954507162, 1e-9);
    EXPECT_EQ(mol.GetIntData("Charge"), 0);
    EXPECT_EQ(mol.GetIntData("Multiplicity"), 1);
    EXPECT_EQ(mol.GetIntData("Number of electrons"), 18);
    EXPECT_EQ(mol.GetStringData("Basis"), "6-311+G(2d,p)");
    EXPECT_FALSE(src.next(mol));
}

class FchkReject : public ::testing::TestWithParam<const char*> {};

TEST_P(FchkReject, MalformedFixtureThrowsFormatError) {
    oeio::builtin::FchkMolSource src(data_path(GetParam()));
    OEChem::OEGraphMol mol;
    EXPECT_THROW(src.next(mol), oeio::FormatError);
}

INSTANTIATE_TEST_SUITE_P(
    Malformed, FchkReject,
    ::testing::Values(
        "fchk_missing_natom.fchk",
        "fchk_missing_atomicnums.fchk",
        "fchk_missing_coords.fchk",
        "fchk_coord_count_mismatch.fchk",
        "fchk_atomicnum_count_mismatch.fchk",
        "fchk_bad_z.fchk",
        "fchk_decimal_z.fchk",
        "fchk_exp_z.fchk",
        "fchk_suffix_z.fchk",
        "fchk_nonfinite_energy.fchk",
        "fchk_decimal_natom.fchk",
        "fchk_bad_charge.fchk",
        "fchk_decimal_ncount.fchk",
        "fchk_oversized_ncount.fchk",
        "fchk_bad_typecode.fchk",
        "fchk_coord_overflow.fchk",
        "fchk_resync_break.fchk",
        "fchk_interior_blank.fchk",
        "fchk_truncated_skip.fchk"));

TEST(FchkReader, OptionalScalarsAbsentStillReadsGeometry) {
    oeio::builtin::FchkMolSource src(data_path("fchk_no_optional.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 2u);
    EXPECT_FALSE(mol.HasData("Total Energy"));
    EXPECT_FALSE(mol.HasData("Charge"));
}

TEST(FchkReader, UnrecognizedScalarWithOddValueIsIgnored) {
    // An unrecognized scalar is line-skipped without value parsing, so a
    // syntactically odd value must NOT raise; only consumed scalars can.
    oeio::builtin::FchkMolSource src(data_path("fchk_unknown_scalar_bad_value.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 2u);
    EXPECT_FALSE(mol.HasData("Job Status"));
}

TEST(FchkReader, SkipsLogicalArray) {
    oeio::builtin::FchkMolSource src(data_path("fchk_skip_logical.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 2u);
}

TEST(FchkReader, SkipsHollerithArray) {
    oeio::builtin::FchkMolSource src(data_path("fchk_skip_hollerith.fchk"));
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(src.next(mol));
    EXPECT_EQ(mol.NumAtoms(), 2u);
}

}  // namespace
