#include <gtest/gtest.h>

#include "oeio/cube_grid.h"
#include "oeio/cube_handler.h"
#include "oeio/exceptions.h"
#include "oeio/format_handler.h"
#include "oeio/mol_range.h"
#include "oeio/write.h"

#include <oechem.h>
#include <oegrid.h>

#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace oeio {
namespace test {

namespace {
std::string data_path(const char* name) {
    return std::string(OEIO_TEST_DATA_DIR) + "/" + name;
}
}

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

TEST(CubeGrid, NonFiniteRejected) {
    oeio::cube::CubeAxes ax{};
    ax.nvox[0] = ax.nvox[1] = ax.nvox[2] = 3;
    ax.vec[0][0] = std::numeric_limits<double>::quiet_NaN();  // NaN diagonal
    ax.vec[1][1] = 0.5; ax.vec[2][2] = 0.5;
    double spacing = 0.0;
    EXPECT_FALSE(oeio::cube::is_axis_aligned_uniform(ax, 1e-6, spacing));

    oeio::cube::CubeAxes ax2{};
    ax2.nvox[0] = ax2.nvox[1] = ax2.nvox[2] = 3;
    ax2.vec[0][0] = 0.5; ax2.vec[1][1] = 0.5; ax2.vec[2][2] = 0.5;
    ax2.vec[0][1] = std::numeric_limits<double>::infinity();  // Inf off-diagonal
    double spacing2 = 0.0;
    EXPECT_FALSE(oeio::cube::is_axis_aligned_uniform(ax2, 1e-6, spacing2));
}

TEST(CubeReader, SingleGridReadsMolAndGrid) {
    oeio::builtin::CubeMolSource src(data_path("single_grid.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    ASSERT_TRUE(src.next(mol, grids));
    EXPECT_EQ(mol.NumAtoms(), 1u);
    ASSERT_EQ(grids.size(), 1u);
    EXPECT_EQ(grids[0].GetXDim(), 2u);
    EXPECT_EQ(grids[0].GetYDim(), 2u);
    EXPECT_EQ(grids[0].GetZDim(), 2u);
    // spacing 0.5 Bohr -> Angstrom
    EXPECT_NEAR(grids[0].GetSpacing(), 0.5 * oeio::cube::BOHR_TO_ANGSTROM, 1e-5);
    EXPECT_FLOAT_EQ(grids[0].GetValues()[0], 0.0f);
    EXPECT_FLOAT_EQ(grids[0].GetValues()[7], 7.0f);
    EXPECT_FALSE(src.next(mol, grids));  // one record only
}

TEST(CubeReader, MoCubeReadsThreeGrids) {
    oeio::builtin::CubeMolSource src(data_path("mo_3grid.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    ASSERT_TRUE(src.next(mol, grids));
    ASSERT_EQ(grids.size(), 3u);
    // orbital-fastest de-interleave: grid0 = {0,10}, grid1 = {1,11}, grid2 = {2,12}
    EXPECT_FLOAT_EQ(grids[0].GetValues()[0], 0.0f);
    EXPECT_FLOAT_EQ(grids[0].GetValues()[1], 10.0f);
    EXPECT_FLOAT_EQ(grids[1].GetValues()[0], 1.0f);
    EXPECT_FLOAT_EQ(grids[2].GetValues()[1], 12.0f);
}

TEST(CubeReader, ReadIntoFillsMinKN) {
    oeio::builtin::CubeMolSource src(data_path("mo_3grid.cube"));
    OEChem::OEMol mol;
    OESystem::OEScalarGrid g0;
    std::vector<OESystem::OEScalarGrid*> grids = { &g0 };  // K=1, N=3
    int n = -1;
    ASSERT_TRUE(src.next(mol, grids, &n));
    EXPECT_EQ(n, 3);
    EXPECT_FLOAT_EQ(g0.GetValues()[0], 0.0f);
}

TEST(CubeReader, SkewedRejected) {
    oeio::builtin::CubeMolSource src(data_path("skewed.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

TEST(CubeReader, AnisotropicRejected) {
    oeio::builtin::CubeMolSource src(data_path("anisotropic.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

TEST(CubeReader, ReflectedRejected) {
    oeio::builtin::CubeMolSource src(data_path("reflected.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

TEST(CubeReader, TruncatedVolumetricBlockRejected) {
    oeio::builtin::CubeMolSource src(data_path("truncated.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

TEST(CubeReader, AngstromUnitsNoBohrConversion) {
    // Negative voxel counts => coordinates already in Angstrom; spacing stays 0.5.
    oeio::builtin::CubeMolSource src(data_path("angstrom.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    ASSERT_TRUE(src.next(mol, grids));
    ASSERT_EQ(grids.size(), 1u);
    EXPECT_NEAR(grids[0].GetSpacing(), 0.5f, 1e-5);  // NOT scaled by Bohr factor
    // Atom at (1,0,0) Angstrom stays at x=1.0 (no Bohr->A scaling).
    float xyz[3];
    OESystem::OEIter<const OEChem::OEAtomBase> ai = mol.GetAtoms();
    ASSERT_TRUE(bool(ai));
    mol.GetCoords(&*ai, xyz);
    EXPECT_NEAR(xyz[0], 1.0f, 1e-4);
}

TEST(CubeReader, ZeroVoxelCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("zero_dim.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

TEST(CubeReader, OversizedVoxelCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("oversized_dim.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// Each dimension is within the per-dimension cap (<= 8192), but their product
// exceeds the total-element ceiling. This exercises the total guard, not the
// per-dimension guard, and must throw before any buffer is allocated. The file
// deliberately has no volumetric data, so reaching allocation would hang/OOM.
TEST(CubeReader, OversizedTotalVoxelCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("oversized_total.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A header at exactly the accepted maximum element count (8192*8192*2 ==
// MAX_TOTAL_ELEMENTS) passes the ceiling, but the file supplies no volumetric
// data. Because buffers grow only as values are consumed, this must fail fast
// with a truncated-block FormatError without pre-allocating the full grid.
TEST(CubeReader, BoundaryHeaderTruncatedFailsFast) {
    oeio::builtin::CubeMolSource src(data_path("boundary_truncated.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A negative (or otherwise out-of-range) atomic number must be rejected with a
// FormatError before the narrowing float-to-unsigned cast reaches OpenEye.
TEST(CubeReader, BadAtomicNumberRejected) {
    oeio::builtin::CubeMolSource src(data_path("bad_atomic_number.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A non-finite atom coordinate must be rejected rather than corrupting geometry.
TEST(CubeReader, NonFiniteCoordRejected) {
    oeio::builtin::CubeMolSource src(data_path("nonfinite_coord.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// LONG_MIN as the atom count must be rejected as FormatError, never passed to
// std::labs (which is undefined behavior for LONG_MIN).
TEST(CubeReader, LongMinAtomCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("longmin_atom_count.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// LONG_MIN as a voxel count must be rejected as FormatError before std::labs.
TEST(CubeReader, LongMinVoxelCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("longmin_voxel_count.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A non-finite origin coordinate must be rejected before it reaches the grid.
TEST(CubeReader, NonFiniteOriginRejected) {
    oeio::builtin::CubeMolSource src(data_path("nonfinite_origin.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A coordinate finite as a double (1e308) but overflowing the float range must
// be rejected rather than silently becoming +/-inf after the narrowing cast.
TEST(CubeReader, OverflowCoordRejected) {
    oeio::builtin::CubeMolSource src(data_path("overflow_coord.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A non-finite volumetric value (inf/nan) must be rejected as corrupt data.
TEST(CubeReader, NonFiniteVolumetricValueRejected) {
    oeio::builtin::CubeMolSource src(data_path("nonfinite_value.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A file that ends before the two required comment lines must fail cleanly.
TEST(CubeReader, TruncatedHeaderRejected) {
    oeio::builtin::CubeMolSource src(data_path("comment_only.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// An atom line with only 4 fields must be rejected, not silently completed by
// splicing a token from the next line. The fixture supplies a compensating
// extra volumetric value so a non-line-anchored parser would parse "success".
TEST(CubeReader, ShortAtomLineRejected) {
    oeio::builtin::CubeMolSource src(data_path("short_atom_line.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// An atom line with an extra trailing field must be rejected as malformed
// rather than silently ignored.
TEST(CubeReader, ExtraAtomFieldRejected) {
    oeio::builtin::CubeMolSource src(data_path("extra_atom_field.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A fractional atom count (e.g. 1.25) must be rejected, not read as int 1 with
// the ".25" suffix spliced into the first origin value.
TEST(CubeReader, FractionalAtomCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("fractional_atom_count.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

// A fractional voxel count (e.g. 2.5) must be rejected, not read as int 2 with
// the ".5" suffix spliced into the first axis-vector value.
TEST(CubeReader, FractionalVoxelCountRejected) {
    oeio::builtin::CubeMolSource src(data_path("fractional_voxel_count.cube"));
    OEChem::OEMol mol;
    std::vector<OESystem::OEScalarGrid> grids;
    EXPECT_THROW(src.next(mol, grids), oeio::FormatError);
}

}  // namespace test
}  // namespace oeio
