"""Shared fixtures and configuration for oeio Python tests."""

import sys
from pathlib import Path

import pytest

# Ensure python/oeio module can be imported
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root / "python"))

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")


@pytest.fixture
def aspirin_mol():
    """Create an aspirin molecule (C9H8O4) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CC(=O)OC1=CC=CC=C1C(=O)O")
    return mol


@pytest.fixture
def ethanol_mol():
    """Create an ethanol molecule (C2H6O) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol


@pytest.fixture
def sdf_file(tmp_path):
    """Write a temp SDF file with 2 molecules and return the path."""
    from openeye import oechem

    path = str(tmp_path / "test.sdf")
    ofs = oechem.oemolostream()
    ofs.open(path)
    ofs.SetFormat(oechem.OEFormat_SDF)

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    mol.SetTitle("ethanol")
    oechem.OEWriteMolecule(ofs, mol)

    mol.Clear()
    oechem.OESmilesToMol(mol, "c1ccccc1")
    mol.SetTitle("benzene")
    oechem.OEWriteMolecule(ofs, mol)

    ofs.close()
    return path


@pytest.fixture
def smi_file(tmp_path):
    """Write a temp SMILES file with 3 molecules and return the path."""
    from openeye import oechem

    path = str(tmp_path / "test.smi")
    ofs = oechem.oemolostream()
    ofs.open(path)
    ofs.SetFormat(oechem.OEFormat_SMI)

    for smi, title in [("CCO", "ethanol"), ("c1ccccc1", "benzene"),
                        ("C", "methane")]:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        mol.SetTitle(title)
        oechem.OEWriteMolecule(ofs, mol)

    ofs.close()
    return path
