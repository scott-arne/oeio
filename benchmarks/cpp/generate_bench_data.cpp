/// \file generate_bench_data.cpp
/// \brief Generate benchmark data files with N molecules.

#include <cstdlib>
#include <iostream>
#include <string>

#include <oechem.h>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <count> <output_path>\n";
        return 1;
    }

    int count = std::atoi(argv[1]);
    std::string path = argv[2];

    OEChem::oemolostream ofs;
    if (!ofs.open(path)) {
        std::cerr << "Failed to open " << path << " for writing\n";
        return 1;
    }

    // Template molecules to cycle through
    const char* smiles[] = {
        "CCO",                      // ethanol (3 heavy atoms)
        "c1ccccc1",                 // benzene (6 heavy atoms)
        "CC(=O)O",                  // acetic acid (4 heavy atoms)
        "CC(=O)OC1=CC=CC=C1C(=O)O" // aspirin (13 heavy atoms)
    };
    constexpr int num_templates = 4;

    OEChem::OEGraphMol mol;
    for (int i = 0; i < count; ++i) {
        mol.Clear();
        OEChem::OESmilesToMol(mol, smiles[i % num_templates]);
        mol.SetTitle(("mol_" + std::to_string(i)).c_str());
        OEChem::OEWriteMolecule(ofs, mol);
    }

    ofs.close();
    std::cout << "Generated " << count << " molecules in " << path << "\n";
    return 0;
}
