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

TEST(Serialize, OebBytesRoundTrip) {
    OEGraphMol mol = make_ethanol();
    const std::string b = oeio::mol_to_bytes(mol);  // default OEB
    ASSERT_FALSE(b.empty());
    OEGraphMol back;
    ASSERT_TRUE(oeio::mol_from_bytes(back, b));
    EXPECT_EQ(back.NumAtoms(), 3u);
    EXPECT_STREQ(back.GetTitle(), "ethanol");
}

TEST(Serialize, GzipSdfBytesRoundTrip) {
    OEGraphMol mol = make_ethanol();
    const std::string b = oeio::mol_to_bytes(mol, "sdf", /*flavor=*/0, /*gzip=*/true);
    OEGraphMol back;
    ASSERT_TRUE(oeio::mol_from_bytes(back, b, "sdf", /*flavor=*/0, /*gzip=*/true));
    EXPECT_EQ(back.NumAtoms(), 3u);
}

TEST(Serialize, MultiConformerOebPreservesConformers) {
    OEMol mol;
    OESmilesToMol(mol, "CCO");
    OEAddExplicitHydrogens(mol);
    // Create a second conformer by duplicating the first.
    OEMol src(mol);
    OESystem::OEIter<OEConfBase> iter = src.GetConfs();
    while (iter) {
        mol.NewConf(&*iter);
        ++iter;
    }
    ASSERT_GT(mol.NumConfs(), 1u);

    const std::string b = oeio::mol_to_bytes(mol);  // OEB
    OEMol back;
    ASSERT_TRUE(oeio::mol_from_bytes(back, b));
    EXPECT_EQ(back.NumConfs(), mol.NumConfs());
}
