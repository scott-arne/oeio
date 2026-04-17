"""Tests for oeio Python bindings."""

import pytest
from openeye import oechem


class TestRead:
    """Test oeio.read() molecule reading."""

    def test_read_sdf(self, sdf_file):
        """Read SDF and verify molecule count."""
        import oeio

        mols = list(oeio.read(sdf_file))
        assert len(mols) == 2

    def test_read_smi(self, smi_file):
        """Read SMILES and verify molecule count."""
        import oeio

        mols = list(oeio.read(smi_file))
        assert len(mols) == 3

    def test_read_returns_oegraphmol(self, sdf_file):
        """Returned objects are oechem.OEGraphMol instances."""
        import oeio

        for mol in oeio.read(sdf_file):
            assert isinstance(mol, oechem.OEGraphMol)

    def test_read_preserves_titles(self, sdf_file):
        """Molecule titles are preserved through read."""
        import oeio

        titles = [mol.GetTitle() for mol in oeio.read(sdf_file)]
        assert titles == ["ethanol", "benzene"]

    def test_read_preserves_atoms(self, sdf_file):
        """Atom counts are correct after reading."""
        import oeio

        atom_counts = [mol.NumAtoms() for mol in oeio.read(sdf_file)]
        assert atom_counts == [3, 6]  # ethanol, benzene heavy atoms


class TestWrite:
    """Test oeio.write() molecule writing."""

    def test_write_sdf(self, sdf_file, tmp_path):
        """Write SDF via context manager, read back and verify."""
        import oeio

        out_path = str(tmp_path / "output.sdf")
        with oeio.write(out_path) as writer:
            for mol in oeio.read(sdf_file):
                writer.add(mol)

        # Verify by reading back via OEChem
        ifs = oechem.oemolistream()
        assert ifs.open(out_path)
        count = 0
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            count += 1
        assert count == 2

    def test_write_context_manager(self, sdf_file, tmp_path):
        """Context manager protocol works correctly."""
        import oeio

        out_path = str(tmp_path / "ctx.sdf")
        writer = oeio.write(out_path)
        assert hasattr(writer, "__enter__")
        assert hasattr(writer, "__exit__")

        with writer as w:
            assert w is writer
            for mol in oeio.read(sdf_file):
                w.add(mol)

    def test_roundtrip_preserves_titles(self, sdf_file, tmp_path):
        """Titles survive a read-write-read roundtrip."""
        import oeio

        out_path = str(tmp_path / "roundtrip.sdf")
        with oeio.write(out_path) as writer:
            for mol in oeio.read(sdf_file):
                writer.add(mol)

        titles = [mol.GetTitle() for mol in oeio.read(out_path)]
        assert titles == ["ethanol", "benzene"]


class TestFilter:
    """Test oeio.filter() molecule filtering."""

    def test_filter_reduces_count(self, smi_file):
        """Filter removes molecules that don't match predicate."""
        import oeio

        heavy = list(oeio.filter(
            oeio.read(smi_file),
            lambda mol: mol.NumAtoms() > 1
        ))
        # ethanol (3) and benzene (6) pass; methane (1) fails
        assert len(heavy) == 2

    def test_filter_preserves_data(self, smi_file):
        """Filtered molecules retain their data."""
        import oeio

        titles = [mol.GetTitle() for mol in oeio.filter(
            oeio.read(smi_file),
            lambda mol: mol.NumAtoms() > 1
        )]
        assert titles == ["ethanol", "benzene"]


class TestTransform:
    """Test oeio.transform() molecule transformation."""

    def test_transform_modifies(self, smi_file):
        """Transform applies function to each molecule."""
        import oeio

        def add_prefix(mol):
            mol.SetTitle("mod_" + mol.GetTitle())

        titles = [mol.GetTitle() for mol in oeio.transform(
            oeio.read(smi_file), add_prefix
        )]
        assert titles == ["mod_ethanol", "mod_benzene", "mod_methane"]

    def test_transform_preserves_count(self, smi_file):
        """Transform does not change molecule count."""
        import oeio

        result = list(oeio.transform(
            oeio.read(smi_file),
            lambda mol: None  # no-op
        ))
        assert len(result) == 3


class TestFormats:
    """Test oeio.formats() format listing."""

    def test_formats_nonempty(self):
        """formats() returns at least one registered format."""
        import oeio

        fmts = oeio.formats()
        assert len(fmts) >= 1

    def test_formats_contains_oechem(self):
        """OEChem handler is registered."""
        import oeio

        names = [f.name for f in oeio.formats()]
        assert "OEChem" in names

    def test_format_info_fields(self):
        """FormatInfo objects have expected attributes."""
        import oeio

        fmt = oeio.formats()[0]
        assert hasattr(fmt, "name")
        assert hasattr(fmt, "extensions")
        assert hasattr(fmt, "description")
        assert hasattr(fmt, "supports_read")
        assert hasattr(fmt, "supports_write")


class TestVersion:
    """Test version information."""

    def test_version_string(self):
        """__version__ is accessible and well-formed."""
        import oeio

        assert hasattr(oeio, "__version__")
        assert oeio.__version__ == "0.1.0"

    def test_version_info(self):
        """__version_info__ is a tuple of integers."""
        import oeio

        assert hasattr(oeio, "__version_info__")
        assert oeio.__version_info__ == (0, 1, 0)


class TestBatchRead:
    """Test oeio.batch_read() batch molecule reading."""

    def test_batch_read_returns_all_molecules(self, sdf_file):
        """batch_read yields same molecules as read."""
        import oeio

        mols = list(oeio.batch_read(sdf_file))
        assert len(mols) == 2

    def test_batch_read_preserves_titles(self, sdf_file):
        """Titles preserved through batch reading."""
        import oeio

        titles = [mol.GetTitle() for mol in oeio.batch_read(sdf_file)]
        assert titles == ["ethanol", "benzene"]

    def test_batch_read_returns_oegraphmol(self, sdf_file):
        """batch_read yields oechem.OEGraphMol instances."""
        import oeio

        for mol in oeio.batch_read(sdf_file):
            assert isinstance(mol, oechem.OEGraphMol)

    def test_batch_read_small_batch_size(self, smi_file):
        """batch_read with batch_size=1 still works."""
        import oeio

        mols = list(oeio.batch_read(smi_file, batch_size=1))
        assert len(mols) == 3

    def test_batch_read_large_batch_size(self, sdf_file):
        """batch_read with batch_size larger than file still works."""
        import oeio

        mols = list(oeio.batch_read(sdf_file, batch_size=10000))
        assert len(mols) == 2
