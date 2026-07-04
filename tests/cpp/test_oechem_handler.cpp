#include <gtest/gtest.h>

#include "oeio/oeio.h"

#include <oechem.h>

#include <filesystem>
#include <string>
#include <vector>

namespace oeio {
namespace test {

class OEChemHandlerTest : public ::testing::Test {
protected:
    void SetUp() override {
        tmpdir_ = std::filesystem::temp_directory_path() /
                  ("oeio_test_" + std::to_string(getpid()));
        std::filesystem::create_directories(tmpdir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(tmpdir_);
    }

    std::string write_temp_smi() {
        std::string path = (tmpdir_ / "test_mols.smi").string();
        OEChem::oemolostream ofs;
        if (!ofs.open(path)) {
            throw std::runtime_error("Failed to open " + path + " for writing");
        }
        ofs.SetFormat(OEChem::OEFormat::SMI);

        OEChem::OEGraphMol mol;
        OEChem::OESmilesToMol(mol, "CCO");
        mol.SetTitle("ethanol");
        OEChem::OEWriteMolecule(ofs, mol);

        mol.Clear();
        OEChem::OESmilesToMol(mol, "CC(=O)O");
        mol.SetTitle("acetic_acid");
        OEChem::OEWriteMolecule(ofs, mol);

        ofs.close();
        return path;
    }

    std::string write_temp_sdf() {
        std::string path = (tmpdir_ / "test_mols.sdf").string();
        OEChem::oemolostream ofs;
        if (!ofs.open(path)) {
            throw std::runtime_error("Failed to open " + path + " for writing");
        }
        ofs.SetFormat(OEChem::OEFormat::SDF);

        OEChem::OEGraphMol mol;
        OEChem::OESmilesToMol(mol, "CCO");
        mol.SetTitle("ethanol");
        OEChem::OEWriteMolecule(ofs, mol);

        mol.Clear();
        OEChem::OESmilesToMol(mol, "CC(=O)O");
        mol.SetTitle("acetic_acid");
        OEChem::OEWriteMolecule(ofs, mol);

        mol.Clear();
        OEChem::OESmilesToMol(mol, "c1ccccc1");
        mol.SetTitle("benzene");
        OEChem::OEWriteMolecule(ofs, mol);

        ofs.close();
        return path;
    }

    std::filesystem::path tmpdir_;
};

TEST_F(OEChemHandlerTest, ReadSMILES) {
    std::string path = write_temp_smi();

    auto range = oeio::read(path);
    int count = 0;
    std::vector<std::string> titles;

    for (auto& mol : range) {
        count++;
        titles.push_back(mol.GetTitle());
    }

    EXPECT_EQ(count, 2);
    ASSERT_EQ(titles.size(), 2);
    EXPECT_EQ(titles[0], "ethanol");
    EXPECT_EQ(titles[1], "acetic_acid");
}

TEST_F(OEChemHandlerTest, ReadSDF) {
    std::string path = write_temp_sdf();

    auto range = oeio::read(path);
    int count = 0;

    for (auto& mol : range) {
        count++;
    }

    EXPECT_EQ(count, 3);
}

TEST_F(OEChemHandlerTest, WriteSDF) {
    std::string out_path = (tmpdir_ / "written.sdf").string();

    auto writer = oeio::write(out_path);

    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, "CCO");
    mol.SetTitle("ethanol");
    writer.append(mol);

    mol.Clear();
    OEChem::OESmilesToMol(mol, "CC(=O)O");
    mol.SetTitle("acetic_acid");
    writer.append(mol);

    writer.close();

    // Read back with OEChem to verify
    OEChem::oemolistream ifs;
    ASSERT_TRUE(ifs.open(out_path));

    int count = 0;
    OEChem::OEGraphMol read_mol;
    while (OEChem::OEReadMolecule(ifs, read_mol)) {
        count++;
    }

    EXPECT_EQ(count, 2);
}

TEST_F(OEChemHandlerTest, RoundTrip) {
    std::string in_path = write_temp_sdf();
    std::string out_path = (tmpdir_ / "roundtrip.sdf").string();

    // Read molecules and collect their properties
    std::vector<std::string> original_titles;
    std::vector<unsigned int> original_atom_counts;

    {
        auto range = oeio::read(in_path);
        for (auto& mol : range) {
            original_titles.push_back(mol.GetTitle());
            original_atom_counts.push_back(mol.NumAtoms());
        }
    }

    // Write molecules
    {
        auto range = oeio::read(in_path);
        auto writer = oeio::write(out_path);
        for (auto& mol : range) {
            writer.append(mol);
        }
        writer.close();
    }

    // Read back and verify
    std::vector<std::string> roundtrip_titles;
    std::vector<unsigned int> roundtrip_atom_counts;

    {
        auto range = oeio::read(out_path);
        for (auto& mol : range) {
            roundtrip_titles.push_back(mol.GetTitle());
            roundtrip_atom_counts.push_back(mol.NumAtoms());
        }
    }

    ASSERT_EQ(original_titles.size(), roundtrip_titles.size());
    ASSERT_EQ(original_atom_counts.size(), roundtrip_atom_counts.size());

    for (size_t i = 0; i < original_titles.size(); i++) {
        EXPECT_EQ(original_titles[i], roundtrip_titles[i]);
        EXPECT_EQ(original_atom_counts[i], roundtrip_atom_counts[i]);
    }
}

TEST_F(OEChemHandlerTest, ReadIntoOEMolBase) {
    std::string path = write_temp_smi();

    auto* handler = oeio::FormatRegistry::instance().lookup(path);
    ASSERT_NE(handler, nullptr);
    auto source = handler->make_reader(path, std::any{});
    ASSERT_NE(source, nullptr);

    // Read into OEMolBase& (the new virtual overload)
    OEChem::OEGraphMol container;
    OEChem::OEMolBase& mol_base = container;
    int count = 0;
    while (source->next(mol_base)) {
        EXPECT_GT(mol_base.NumAtoms(), 0);
        ++count;
        container.Clear();
    }
    EXPECT_EQ(count, 2);
}

TEST_F(OEChemHandlerTest, NextOEMolReadsMultiConformer) {
    // Write a genuine 3-conformer molecule to an OEB file.
    std::string path = (tmpdir_ / "multiconf.oeb").string();
    {
        OEChem::OEMol mc;
        OEChem::OESmilesToMol(mc, "CCO");
        mc.SetTitle("ethanol_mc");
        OEChem::OEAddExplicitHydrogens(mc);
        std::vector<float> coords(mc.NumAtoms() * 3);
        mc.GetCoords(coords.data());
        mc.NewConf(coords.data());
        mc.NewConf(coords.data());
        ASSERT_EQ(mc.NumConfs(), 3u);
        OEChem::oemolostream ofs;
        ASSERT_TRUE(ofs.open(path));
        OEChem::OEWriteMolecule(ofs, mc);
        ofs.close();
    }

    auto* handler = oeio::FormatRegistry::instance().lookup(path);
    ASSERT_NE(handler, nullptr);
    auto source = handler->make_reader(path, std::any{});
    ASSERT_NE(source, nullptr);

    // Read into an OEMol bound as OEMolBase& — must preserve all 3 conformers.
    OEChem::OEMol mol;
    OEChem::OEMolBase& ref = mol;
    ASSERT_TRUE(source->next(ref));
    EXPECT_EQ(mol.NumConfs(), 3u);
}

TEST_F(OEChemHandlerTest, NextOEGraphMolFlattensMultiConformer) {
    // Same 3-conformer OEB as above: through an OEGraphMol path each conformer
    // surfaces as its own single-conformer record (3 records), confirming the
    // dispatch does NOT route OEGraphMol through the multi-conformer overload.
    std::string path = (tmpdir_ / "multiconf.oeb").string();
    {
        OEChem::OEMol mc;
        OEChem::OESmilesToMol(mc, "CCO");
        mc.SetTitle("ethanol_mc");
        OEChem::OEAddExplicitHydrogens(mc);
        std::vector<float> coords(mc.NumAtoms() * 3);
        mc.GetCoords(coords.data());
        mc.NewConf(coords.data());
        mc.NewConf(coords.data());
        ASSERT_EQ(mc.NumConfs(), 3u);
        OEChem::oemolostream ofs;
        ASSERT_TRUE(ofs.open(path));
        OEChem::OEWriteMolecule(ofs, mc);
        ofs.close();
    }
    auto* handler = oeio::FormatRegistry::instance().lookup(path);
    ASSERT_NE(handler, nullptr);
    auto source = handler->make_reader(path, std::any{});
    ASSERT_NE(source, nullptr);

    OEChem::OEGraphMol mol;
    OEChem::OEMolBase& ref = mol;
    int count = 0;
    while (source->next(ref)) { ++count; }
    EXPECT_EQ(count, 3);  // one single-conformer record per conformer
}

}  // namespace test
}  // namespace oeio
