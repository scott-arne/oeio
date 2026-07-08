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

struct FchkRejectCase {
    const char* fixture;
    const char* expected_message_substring;
};

class FchkReject : public ::testing::TestWithParam<FchkRejectCase> {};

TEST_P(FchkReject, MalformedFixtureThrowsFormatErrorWithExpectedMessage) {
    const auto& param = GetParam();
    oeio::builtin::FchkMolSource src(data_path(param.fixture));
    OEChem::OEGraphMol mol;

    try {
        src.next(mol);
        FAIL() << "Expected FormatError for " << param.fixture;
    } catch (const oeio::FormatError& e) {
        std::string message(e.what());
        EXPECT_NE(message.find(param.expected_message_substring), std::string::npos)
            << "For " << param.fixture << ":\n"
            << "  Expected substring: \"" << param.expected_message_substring << "\"\n"
            << "  Actual message: \"" << message << "\"";
    }
}

INSTANTIATE_TEST_SUITE_P(
    Malformed, FchkReject,
    ::testing::Values(
        FchkRejectCase{"fchk_missing_natom.fchk", "missing 'Number of atoms'"},
        FchkRejectCase{"fchk_missing_atomicnums.fchk", "missing 'Atomic numbers'"},
        FchkRejectCase{"fchk_missing_coords.fchk", "missing 'Current cartesian coordinates'"},
        FchkRejectCase{"fchk_coord_count_mismatch.fchk", "coordinate count does not match"},
        FchkRejectCase{"fchk_atomicnum_count_mismatch.fchk", "atomic-number count does not match"},
        FchkRejectCase{"fchk_bad_z.fchk", "atomic number out of range"},
        FchkRejectCase{"fchk_decimal_z.fchk", "malformed FCHK atomic number"},
        FchkRejectCase{"fchk_exp_z.fchk", "malformed FCHK atomic number"},
        FchkRejectCase{"fchk_suffix_z.fchk", "malformed FCHK atomic number"},
        FchkRejectCase{"fchk_nonfinite_energy.fchk", "malformed FCHK total energy"},
        FchkRejectCase{"fchk_decimal_natom.fchk", "malformed FCHK atom count"},
        FchkRejectCase{"fchk_bad_charge.fchk", "malformed FCHK charge"},
        FchkRejectCase{"fchk_decimal_ncount.fchk", "malformed FCHK array count"},
        FchkRejectCase{"fchk_oversized_ncount.fchk", "array count out of range"},
        FchkRejectCase{"fchk_bad_typecode.fchk", "unrecognized type code"},
        FchkRejectCase{"fchk_coord_overflow.fchk", "non-finite or out of range"},
        FchkRejectCase{"fchk_resync_break.fchk", "lost record sync"},
        FchkRejectCase{"fchk_interior_blank.fchk", "unexpected blank line"},
        FchkRejectCase{"fchk_truncated_skip.fchk", "truncated payload"},
        FchkRejectCase{"fchk_overflow_charge.fchk", "FCHK charge out of range"}));

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
