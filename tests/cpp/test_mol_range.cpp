#include <gtest/gtest.h>

#include "oeio/oeio.h"

#include <oechem.h>

#include <filesystem>
#include <string>

namespace oeio {
namespace test {

class MolRangeTest : public ::testing::Test {
protected:
    void SetUp() override {
        tmpdir_ = std::filesystem::temp_directory_path() /
                  ("oeio_test_" + std::to_string(getpid()));
        std::filesystem::create_directories(tmpdir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(tmpdir_);
    }

    std::string write_temp_file(int num_mols) {
        std::string path = (tmpdir_ / "test_mols.smi").string();
        OEChem::oemolostream ofs;
        if (!ofs.open(path)) {
            throw std::runtime_error("Failed to open " + path + " for writing");
        }
        ofs.SetFormat(OEChem::OEFormat::SMI);

        for (int i = 0; i < num_mols; i++) {
            OEChem::OEGraphMol mol;
            OEChem::OESmilesToMol(mol, "CCO");
            mol.SetTitle("mol_" + std::to_string(i));
            OEChem::OEWriteMolecule(ofs, mol);
        }

        ofs.close();
        return path;
    }

    std::filesystem::path tmpdir_;
};

TEST_F(MolRangeTest, IterateProducesCorrectCount) {
    std::string path = write_temp_file(2);

    auto range = oeio::read(path);
    int count = 0;

    for (auto& mol : range) {
        count++;
    }

    EXPECT_EQ(count, 2);
}

TEST_F(MolRangeTest, EmptyFile) {
    std::string path = write_temp_file(0);

    auto range = oeio::read(path);
    int count = 0;

    for (auto& mol : range) {
        count++;
    }

    EXPECT_EQ(count, 0);
}

TEST_F(MolRangeTest, MoveSemantics) {
    std::string path = write_temp_file(2);

    // Test move construction
    MolRange range = oeio::read(path);
    int count = 0;

    for (auto& mol : range) {
        count++;
    }

    EXPECT_EQ(count, 2);
}

TEST_F(MolRangeTest, MoleculeDataPreserved) {
    std::string path = (tmpdir_ / "test_data.smi").string();

    // Write molecules with specific titles and structures
    OEChem::oemolostream ofs;
    ASSERT_TRUE(ofs.open(path));
    ofs.SetFormat(OEChem::OEFormat::SMI);

    OEChem::OEGraphMol mol1;
    OEChem::OESmilesToMol(mol1, "CCO");
    mol1.SetTitle("ethanol");
    OEChem::OEWriteMolecule(ofs, mol1);

    OEChem::OEGraphMol mol2;
    OEChem::OESmilesToMol(mol2, "c1ccccc1");
    mol2.SetTitle("benzene");
    OEChem::OEWriteMolecule(ofs, mol2);

    ofs.close();

    // Read back and verify
    auto range = oeio::read(path);
    std::vector<std::string> titles;
    std::vector<unsigned int> atom_counts;

    for (auto& mol : range) {
        titles.push_back(mol.GetTitle());
        atom_counts.push_back(mol.NumAtoms());
    }

    ASSERT_EQ(titles.size(), 2);
    EXPECT_EQ(titles[0], "ethanol");
    EXPECT_EQ(titles[1], "benzene");

    ASSERT_EQ(atom_counts.size(), 2);
    EXPECT_EQ(atom_counts[0], 3);  // C-C-O
    EXPECT_EQ(atom_counts[1], 6);  // benzene ring
}

}  // namespace test
}  // namespace oeio
