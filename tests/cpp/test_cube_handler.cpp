#include <gtest/gtest.h>

#include "oeio/cube_grid.h"
#include "oeio/format_handler.h"
#include "oeio/mol_range.h"
#include "oeio/write.h"

#include <oechem.h>
#include <oegrid.h>

#include <memory>
#include <vector>

namespace oeio {
namespace test {

// A minimal MolSource that yields exactly one empty molecule and no grids,
// used to prove the defaulted grid overloads delegate to the molecule path.
class OneMolSource : public MolSource {
public:
    // Bring the base grid overloads into scope: declaring next(OEGraphMol&)
    // otherwise HIDES the inherited next(mol, grids, ...) overloads for calls
    // made through the OneMolSource static type (C++ name hiding). This test
    // exercises the DEFAULT base implementations, so we must un-hide them.
    using MolSource::next;

    bool next(OEChem::OEGraphMol& mol) override {
        mol.Clear();
        if (done_) return false;
        done_ = true;
        return true;
    }
private:
    bool done_ = false;
};

TEST(MolSourceGridDefault, GridOverloadDelegatesToMolecule) {
    OneMolSource src;
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid*> grids;  // caller supplies none
    int n = -99;
    EXPECT_TRUE(src.next(mol, grids, &n));   // first record read
    EXPECT_EQ(n, 0);                         // default reports 0 grids
    n = 99;                                  // sentinel: EOF must not touch *num_grids
    EXPECT_FALSE(src.next(mol, grids, &n));  // EOF
    EXPECT_EQ(n, 99);                        // unchanged on EOF (spec contract)
}

TEST(MolSourceGridDefault, OwnedGridOverloadClearsGrids) {
    OneMolSource src;
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids(3);  // pre-populated
    EXPECT_TRUE(src.next(mol, grids));
    EXPECT_TRUE(grids.empty());              // default clears to N=0
    EXPECT_FALSE(src.next(mol, grids));      // EOF
}

// A source yielding one record: one molecule + two 1x1x1 grids (values 10, 20).
class MolWithTwoGridsSource : public MolSource {
public:
    bool next(OEChem::OEGraphMol& mol) override {
        mol.Clear();
        if (done_) return false;
        done_ = true;
        return true;
    }
    bool next(OEChem::OEMolBase& mol, const std::vector<OESystem::OEScalarGrid*>& grids,
              int* num_grids) override {
        mol.Clear();
        if (done_) { return false; }
        done_ = true;
        if (num_grids) *num_grids = 2;
        float v0 = 10.0f, v1 = 20.0f;
        if (grids.size() > 0 && grids[0]) grids[0]->SetValues(&v0, 1);
        if (grids.size() > 1 && grids[1]) grids[1]->SetValues(&v1, 1);
        return true;
    }
    bool next(OEChem::OEMolBase& mol, std::vector<OESystem::OEScalarGrid>& grids) override {
        mol.Clear();
        if (done_) { grids.clear(); return false; }
        done_ = true;
        grids.assign(2, OESystem::OEScalarGrid(1,1,1, 0,0,0, 1.0f));
        float v0 = 10.0f, v1 = 20.0f;
        grids[0].SetValues(&v0, 1);
        grids[1].SetValues(&v1, 1);
        return true;
    }
private:
    bool done_ = false;
};

TEST(MolRangeGrids, ReadIntoFillsMinKN) {
    MolRange range(std::make_unique<MolWithTwoGridsSource>());
    OEChem::OEMol mol;
    OESystem::OEScalarGrid g0;
    std::vector<OESystem::OEScalarGrid*> grids = { &g0 };  // K=1, file N=2
    int n = -1;
    ASSERT_TRUE(range.read_into(mol, grids, &n));
    EXPECT_EQ(n, 2);                          // reports file's N
    EXPECT_FLOAT_EQ(*g0.GetValues(), 10.0f);  // grid 0 filled
    EXPECT_FALSE(range.read_into(mol, grids, &n));  // EOF -> false
}

TEST(MolRangeGrids, WithGridsYieldsAllN) {
    MolRange range(std::make_unique<MolWithTwoGridsSource>());
    int rows = 0;
    for (auto& [mol, grids] : range.with_grids()) {
        (void)mol;
        ++rows;
        ASSERT_EQ(grids.size(), 2u);
        EXPECT_FLOAT_EQ(*grids[0].GetValues(), 10.0f);
        EXPECT_FLOAT_EQ(*grids[1].GetValues(), 20.0f);
    }
    EXPECT_EQ(rows, 1);
}

// A sink that records whether the grid write overload was called and with how
// many grids.
class RecordingSink : public MolSink {
public:
    bool write(const OEChem::OEMolBase&) override { mol_writes++; return true; }
    bool write(const OEChem::OEMolBase&,
               const std::vector<const OESystem::OEScalarGrid*>& grids) override {
        last_grid_count = static_cast<int>(grids.size());
        return true;
    }
    void close() override {}
    int mol_writes = 0;
    int last_grid_count = -1;
};

TEST(WriterGrids, AppendForwardsGrids) {
    auto sink = std::make_unique<RecordingSink>();
    RecordingSink* raw = sink.get();
    Writer writer(std::move(sink));
    OEChem::OEMol mol;
    OESystem::OEScalarGrid g0, g1;
    std::vector<const OESystem::OEScalarGrid*> grids = { &g0, &g1 };
    EXPECT_TRUE(writer.append(mol, grids));
    EXPECT_EQ(raw->last_grid_count, 2);
}

TEST(CubeGrid, AxisAlignedUniformAccepted) {
    oeio::cube::CubeAxes ax{};
    ax.nvox[0] = ax.nvox[1] = ax.nvox[2] = 3;
    ax.vec[0][0] = 0.5; ax.vec[1][1] = 0.5; ax.vec[2][2] = 0.5;  // cubic, positive
    double spacing = 0.0;
    EXPECT_TRUE(oeio::cube::is_axis_aligned_uniform(ax, 1e-6, spacing));
    EXPECT_DOUBLE_EQ(spacing, 0.5);
}

TEST(CubeGrid, SkewedRejected) {
    oeio::cube::CubeAxes ax{};
    ax.nvox[0] = ax.nvox[1] = ax.nvox[2] = 3;
    ax.vec[0][0] = 0.5; ax.vec[0][1] = 0.1;  // off-diagonal -> skew
    ax.vec[1][1] = 0.5; ax.vec[2][2] = 0.5;
    double spacing = 0.0;
    EXPECT_FALSE(oeio::cube::is_axis_aligned_uniform(ax, 1e-6, spacing));
}

TEST(CubeGrid, AnisotropicRejected) {
    oeio::cube::CubeAxes ax{};
    ax.nvox[0] = ax.nvox[1] = ax.nvox[2] = 3;
    ax.vec[0][0] = 0.5; ax.vec[1][1] = 0.6; ax.vec[2][2] = 0.5;  // unequal
    double spacing = 0.0;
    EXPECT_FALSE(oeio::cube::is_axis_aligned_uniform(ax, 1e-6, spacing));
}

TEST(CubeGrid, ReflectedRejected) {
    oeio::cube::CubeAxes ax{};
    ax.nvox[0] = ax.nvox[1] = ax.nvox[2] = 3;
    ax.vec[0][0] = -0.5;  // negative diagonal -> reflected
    ax.vec[1][1] = 0.5; ax.vec[2][2] = 0.5;
    double spacing = 0.0;
    EXPECT_FALSE(oeio::cube::is_axis_aligned_uniform(ax, 1e-6, spacing));
}

TEST(CubeGrid, OriginMidpointRoundTrip) {
    double origin[3] = { 1.0, 2.0, 3.0 };
    int nvox[3] = { 3, 5, 7 };
    double spacing = 0.5;
    double mid[3], origin_back[3];
    oeio::cube::origin_to_mid(origin, nvox, spacing, mid);
    // mid = origin + s*((n-1)/2)
    EXPECT_DOUBLE_EQ(mid[0], 1.0 + 0.5 * 1.0);  // (3-1)/2 = 1
    EXPECT_DOUBLE_EQ(mid[1], 2.0 + 0.5 * 2.0);  // (5-1)/2 = 2
    EXPECT_DOUBLE_EQ(mid[2], 3.0 + 0.5 * 3.0);  // (7-1)/2 = 3
    oeio::cube::mid_to_origin(mid, nvox, spacing, origin_back);
    EXPECT_DOUBLE_EQ(origin_back[0], origin[0]);
    EXPECT_DOUBLE_EQ(origin_back[1], origin[1]);
    EXPECT_DOUBLE_EQ(origin_back[2], origin[2]);
}

}  // namespace test
}  // namespace oeio
