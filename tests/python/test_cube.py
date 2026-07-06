"""Tests for oeio Gaussian CUBE support (Python grid API)."""

import pytest
from openeye import oechem, oegrid


class TestCubeRead:
    def test_next_returns_oemol(self, cube_file):
        import oeio
        with oeio.read(cube_file) as reader:
            mol = next(reader)
        assert isinstance(mol, oechem.OEMol)
        assert mol.NumAtoms() == 1

    def test_next_raises_stopiteration_at_eof(self, cube_file):
        import oeio
        with oeio.read(cube_file) as reader:
            next(reader)  # the single record
            with pytest.raises(StopIteration):
                next(reader)

    def test_iterate_molecule_only(self, cube_file):
        import oeio
        with oeio.read(cube_file) as reader:
            mols = list(reader)
        assert len(mols) == 1

    def test_with_grids_single(self, cube_file):
        import oeio
        with oeio.read(cube_file) as reader:
            rows = list(reader.with_grids())
        assert len(rows) == 1
        mol, grids = rows[0]
        assert isinstance(mol, oechem.OEMol)
        assert len(grids) == 1
        assert isinstance(grids[0], oegrid.OEScalarGrid)
        assert grids[0].GetXDim() == 2
        assert grids[0].GetValues()[7] == 7.0

    def test_with_grids_mo_three(self, mo_cube_file):
        import oeio
        with oeio.read(mo_cube_file) as reader:
            mol, grids = next(iter(reader.with_grids()))
        assert len(grids) == 3
        # Orbital-fastest de-interleave (matches the C++ reader tests).
        assert grids[0].GetValues()[0] == 0.0
        assert grids[0].GetValues()[1] == 10.0
        assert grids[2].GetValues()[1] == 12.0

    def test_read_into_mol_only(self, cube_file):
        import oeio
        mol = oechem.OEMol()
        with oeio.read(cube_file) as reader:
            assert reader.read_into(mol) is True
            assert reader.read_into(mol) is False  # EOF

    def test_read_into_grids_returns_n_and_none_at_eof(self, mo_cube_file):
        import oeio
        mol = oechem.OEMol()
        g0 = oegrid.OEScalarGrid()
        with oeio.read(mo_cube_file) as reader:
            n = reader.read_into(mol, g0)   # K=1, file N=3
            assert n == 3
            assert g0.GetValues()[0] == 0.0
            assert reader.read_into(mol, g0) is None  # EOF -> None

    def test_read_into_grids_fills_min_k_n(self, mo_cube_file):
        import oeio
        # Supply K=2 grids for an N=3 file: both get filled, N reported as 3.
        mol = oechem.OEMol()
        g0 = oegrid.OEScalarGrid()
        g1 = oegrid.OEScalarGrid()
        with oeio.read(mo_cube_file) as reader:
            n = reader.read_into(mol, g0, g1)
            assert n == 3
            assert g0.GetValues()[0] == 0.0   # grid 0
            assert g1.GetValues()[0] == 1.0   # grid 1


class TestCubeWrite:
    def test_write_roundtrip(self, cube_file, tmp_path):
        import oeio
        out = str(tmp_path / "out.cube")
        with oeio.read(cube_file) as reader:
            mol, grids = next(iter(reader.with_grids()))
        with oeio.write(out) as writer:
            writer.append(mol, grids[0])
        with oeio.read(out) as reader:
            mol2, grids2 = next(iter(reader.with_grids()))
        assert len(grids2) == 1
        assert grids2[0].GetValues()[5] == 5.0

    def test_write_multigrid_roundtrip(self, mo_cube_file, tmp_path):
        import oeio
        out = str(tmp_path / "out_mo.cube")
        with oeio.read(mo_cube_file) as reader:
            mol, grids = next(iter(reader.with_grids()))
        assert len(grids) == 3
        with oeio.write(out) as writer:
            writer.append(mol, *grids)
        with oeio.read(out) as reader:
            mol2, grids2 = next(iter(reader.with_grids()))
        assert len(grids2) == 3
        assert grids2[0].GetValues()[1] == 10.0
        assert grids2[2].GetValues()[1] == 12.0

    def test_append_without_grid_raises(self, cube_file, tmp_path):
        import oeio
        out = str(tmp_path / "nogrid.cube")
        with oeio.read(cube_file) as reader:
            mol = next(reader)
        with oeio.write(out) as writer:
            with pytest.raises(oeio.FormatError):
                writer.append(mol)
