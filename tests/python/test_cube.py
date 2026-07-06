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

    def test_append_grids_accepts_lazy_sequence(self, cube_file, tmp_path):
        """A sequence whose __getitem__ returns a fresh wrapper per access must
        round-trip through the raw grid-list typemap. The typemap holds a tuple
        reference for the duration of the C++ call, so the extracted grid
        pointers cannot dangle. Calling the raw ``append_grids`` handle method
        (not the ``append(*grids)`` shim, which materializes a list first)
        exercises the typemap directly — regression for the fixed
        extract-then-DECREF lifetime bug."""
        import oeio
        with oeio.read(cube_file) as reader:
            mol, grids = next(iter(reader.with_grids()))
        base = grids[0]

        class FreshEachGet:
            """Returns a fresh copy of the grid on every __getitem__, so the
            only reference to the returned wrapper is the one the typemap must
            keep alive itself."""
            def __init__(self, src):
                self._src = src

            def __len__(self):
                return 1

            def __getitem__(self, i):
                from openeye import oegrid
                if i != 0:
                    raise IndexError(i)
                return oegrid.OEScalarGrid(self._src)  # deep copy

        out = str(tmp_path / "lazy.cube")
        # write() returns the raw _WriterHandle; call the C++ method directly so
        # the lazy sequence reaches the typemap unmaterialized.
        writer = oeio.write(out)
        assert writer.append_grids(mol, FreshEachGet(base)) is True
        writer.close()
        with oeio.read(out) as reader:
            _, grids2 = next(iter(reader.with_grids()))
        assert grids2[0].GetValues()[5] == 5.0

    def test_append_rejects_non_grid_sequence_element(self, cube_file, tmp_path):
        """A sequence element that is not an OEScalarGrid must raise cleanly
        (TypeError), not crash — the typemap validates each element."""
        import oeio
        out = str(tmp_path / "badseq.cube")
        with oeio.read(cube_file) as reader:
            mol = next(reader)
        with oeio.write(out) as writer:
            with pytest.raises((TypeError, oeio.Error)):
                writer.append(mol, "not a grid")
