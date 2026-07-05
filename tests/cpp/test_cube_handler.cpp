#include <gtest/gtest.h>

#include "oeio/format_handler.h"

#include <oechem.h>
#include <oegrid.h>

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

}  // namespace test
}  // namespace oeio
