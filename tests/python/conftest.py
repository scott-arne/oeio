"""Shared fixtures and configuration for oeio Python tests."""

import pytest

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


@pytest.fixture
def cube_file():
    """Path to the single-grid CUBE test file."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "single_grid.cube")


@pytest.fixture
def mo_cube_file():
    """Path to the 3-orbital MO CUBE test file."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "mo_3grid.cube")


@pytest.fixture
def h2o_dens_cube_file():
    """Path to the real-world water electron-density CUBE file (40^3 grid)."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "h2o-dens.cube")


@pytest.fixture
def multiconf_oeb(tmp_path):
    """Write a temp OEB with one 3-conformer molecule; return the path."""
    from openeye import oechem

    path = str(tmp_path / "multiconf.oeb")
    mc = oechem.OEMol()
    oechem.OESmilesToMol(mc, "CCO")
    mc.SetTitle("ethanol_mc")
    oechem.OEAddExplicitHydrogens(mc)
    coords = oechem.OEFloatArray(mc.NumAtoms() * 3)
    mc.GetCoords(coords)
    mc.NewConf(coords)
    mc.NewConf(coords)
    ofs = oechem.oemolostream()
    ofs.open(path)
    oechem.OEWriteMolecule(ofs, mc)
    ofs.close()
    return path


@pytest.fixture
def fchk_file():
    """Path to the real-world Gaussian FCHK test file (HOF, 3 atoms)."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "example.fchk")


@pytest.fixture
def fchk_min_file():
    """Path to the minimal hand-authored FCHK test file (2 atoms)."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "fchk_min.fchk")


@pytest.fixture
def fchk_bad_z_file():
    """Path to a malformed FCHK file (atomic number out of range)."""
    import pathlib
    return str(pathlib.Path(__file__).parent.parent / "data" / "fchk_bad_z.fchk")
