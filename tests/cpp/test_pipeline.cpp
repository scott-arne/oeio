#include <gtest/gtest.h>

#include "oeio/oeio.h"

#include <oechem.h>

#include <filesystem>
#include <string>
#include <unistd.h>

namespace oeio {
namespace test {

class PipelineTest : public ::testing::Test {
protected:
    void SetUp() override {
        tmpdir_ = std::filesystem::temp_directory_path() /
                  ("oeio_test_" + std::to_string(getpid()));
        std::filesystem::create_directories(tmpdir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(tmpdir_);
    }

    std::string write_test_molecules() {
        std::string path = (tmpdir_ / "input.smi").string();
        OEChem::oemolostream ofs;
        if (!ofs.open(path)) {
            throw std::runtime_error("Failed to open " + path + " for writing");
        }
        ofs.SetFormat(OEChem::OEFormat::SMI);

        // Molecule with 3 heavy atoms
        OEChem::OEGraphMol mol;
        OEChem::OESmilesToMol(mol, "CCO");
        mol.SetTitle("ethanol");
        OEChem::OEWriteMolecule(ofs, mol);

        // Molecule with 6 heavy atoms
        mol.Clear();
        OEChem::OESmilesToMol(mol, "c1ccccc1");
        mol.SetTitle("benzene");
        OEChem::OEWriteMolecule(ofs, mol);

        // Molecule with 1 heavy atom
        mol.Clear();
        OEChem::OESmilesToMol(mol, "C");
        mol.SetTitle("methane");
        OEChem::OEWriteMolecule(ofs, mol);

        // Molecule with 4 heavy atoms
        mol.Clear();
        OEChem::OESmilesToMol(mol, "CC(=O)O");
        mol.SetTitle("acetic_acid");
        OEChem::OEWriteMolecule(ofs, mol);

        ofs.close();
        return path;
    }

    int count_molecules(const std::string& path) {
        auto range = oeio::read(path);
        int count = 0;
        for (auto& mol : range) {
            count++;
        }
        return count;
    }

    std::filesystem::path tmpdir_;
};

TEST_F(PipelineTest, PipeReadToWrite) {
    std::string in_path = write_test_molecules();
    std::string out_path = (tmpdir_ / "output.sdf").string();

    oeio::read(in_path) | oeio::write(out_path);

    int in_count = count_molecules(in_path);
    int out_count = count_molecules(out_path);

    EXPECT_EQ(in_count, out_count);
    EXPECT_EQ(out_count, 4);
}

TEST_F(PipelineTest, FilterReducesCount) {
    std::string in_path = write_test_molecules();
    std::string out_path = (tmpdir_ / "filtered.sdf").string();

    // Filter to molecules with > 2 heavy atoms
    oeio::filter(oeio::read(in_path), [](const OEChem::OEMolBase& mol) {
        return mol.NumAtoms() > 2;
    }) | oeio::write(out_path);

    int out_count = count_molecules(out_path);

    // Should have 3 molecules: ethanol (3), benzene (6), acetic_acid (4)
    // Should exclude: methane (1)
    EXPECT_EQ(out_count, 3);
}

TEST_F(PipelineTest, TransformModifiesMolecules) {
    std::string in_path = write_test_molecules();
    std::string out_path = (tmpdir_ / "transformed.sdf").string();

    // Transform that sets all titles to "modified"
    oeio::transform(oeio::read(in_path), [](OEChem::OEGraphMol& mol) {
        mol.SetTitle("modified");
    }) | oeio::write(out_path);

    // Read back and verify titles
    auto range = oeio::read(out_path);
    std::vector<std::string> titles;
    for (auto& mol : range) {
        titles.push_back(mol.GetTitle());
    }

    ASSERT_EQ(titles.size(), 4);
    for (const auto& title : titles) {
        EXPECT_EQ(title, "modified");
    }
}

TEST_F(PipelineTest, ChainedFilterTransformWrite) {
    std::string in_path = write_test_molecules();
    std::string out_path = (tmpdir_ / "chained.sdf").string();

    // Filter to molecules with > 2 heavy atoms, then add prefix to titles
    oeio::transform(
        oeio::filter(oeio::read(in_path), [](const OEChem::OEMolBase& mol) {
            return mol.NumAtoms() > 2;
        }),
        [](OEChem::OEGraphMol& mol) {
            std::string old_title = mol.GetTitle();
            mol.SetTitle("filtered_" + old_title);
        }
    ) | oeio::write(out_path);

    // Read back and verify
    auto range = oeio::read(out_path);
    std::vector<std::string> titles;
    int count = 0;

    for (auto& mol : range) {
        count++;
        titles.push_back(mol.GetTitle());
    }

    EXPECT_EQ(count, 3);
    ASSERT_EQ(titles.size(), 3);

    // Check that titles have the prefix
    for (const auto& title : titles) {
        EXPECT_TRUE(title.substr(0, 9) == "filtered_");
    }
}

TEST_F(PipelineTest, FilteredReadToOEMolBase) {
    std::string in_path = write_test_molecules();
    std::string out_path = (tmpdir_ / "molbase_filtered.sdf").string();

    // Filter + write pipeline (exercises FilteredMolSource)
    oeio::filter(oeio::read(in_path), [](const OEChem::OEMolBase& mol) {
        return mol.NumAtoms() > 2;
    }) | oeio::write(out_path);

    EXPECT_EQ(count_molecules(out_path), 3);
}

}  // namespace test
}  // namespace oeio
