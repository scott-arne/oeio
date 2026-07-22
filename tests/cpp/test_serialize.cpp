#include <gtest/gtest.h>

#include <oechem.h>

#include "oeio/serialize.h"
#include "oeio/exceptions.h"

using namespace OEChem;

namespace {

OEGraphMol make_ethanol() {
    OEGraphMol mol;
    OESmilesToMol(mol, "CCO");
    mol.SetTitle("ethanol");
    return mol;
}

}  // namespace

TEST(Serialize, StringRoundTripSmiles) {
    OEGraphMol mol = make_ethanol();
    const std::string s = oeio::mol_to_string(mol, "smi");
    ASSERT_FALSE(s.empty());

    OEGraphMol back;
    ASSERT_TRUE(oeio::mol_from_string(back, s, "smi"));
    EXPECT_EQ(back.NumAtoms(), 3u);
    EXPECT_STREQ(back.GetTitle(), "ethanol");
}

TEST(Serialize, StringRoundTripSdfPreservesSdTag) {
    OEGraphMol mol = make_ethanol();
    OESetSDData(mol, "prop", "42");
    const std::string s = oeio::mol_to_string(mol, "sdf");

    OEGraphMol back;
    ASSERT_TRUE(oeio::mol_from_string(back, s, "sdf"));
    EXPECT_TRUE(OEHasSDData(back, "prop"));
    EXPECT_EQ(OEGetSDData(back, "prop"), "42");
}

TEST(Serialize, UnknownFormatThrows) {
    OEGraphMol mol = make_ethanol();
    EXPECT_THROW(oeio::mol_to_string(mol, "cube"), oeio::FormatError);
    EXPECT_THROW(oeio::mol_to_string(mol, "not_a_format"), oeio::FormatError);
}

TEST(Serialize, SingleConfOEMolKeepsMoleculeTitle) {
    // A single-conformer OEMol (OEMCMolBase) must serialize via the OEMolBase
    // path so the molecule-level title survives, not a per-conformer record.
    OEMol mol;
    OESmilesToMol(mol, "CCO");
    mol.SetTitle("mol-level-title");
    const std::string s = oeio::mol_to_string(mol, "sdf");
    OEGraphMol back;
    ASSERT_TRUE(oeio::mol_from_string(back, s, "sdf"));
    EXPECT_STREQ(back.GetTitle(), "mol-level-title");
}
