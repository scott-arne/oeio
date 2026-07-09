"""Tests for oeio Gaussian FCHK support (read-only, molecule + QM scalars)."""

import pytest
from openeye import oechem


def _openeye_shared_build():
    """Return True when oeio is linked *shared* against OpenEye.

    A shared link means oeio and Python's ``oechem`` load the same
    ``liboesystem`` image and therefore share one generic-data tag registry, so
    a string tag set from C++ (``mol.SetData("Total Energy", v)``) resolves to
    the same integer Python's ``GetData("Total Energy")`` uses. A static wheel
    (the Windows wheel is static because the OpenEye C++ SDK ships static-only
    ``.lib`` there) gives oeio its own registry, so C++-set string tags are not
    visible from Python. On an unknown build we return True so a genuine
    regression on a shared platform is never hidden by this guard.
    """
    try:
        from oeio import _build_info
    except Exception:
        return True
    return getattr(_build_info, "OPENEYE_LIBRARY_TYPE", "STATIC") == "SHARED"


# FCHK typed QM scalars cross into Python only on a shared OpenEye build; skip
# their assertions on static wheels. Molecule geometry, title, and elements are
# unaffected and remain covered on every platform.
requires_shared_openeye = pytest.mark.skipif(
    not _openeye_shared_build(),
    reason="FCHK typed QM scalars require a shared OpenEye link; the static "
           "wheel (e.g. Windows) splits the liboesystem tag registry so a "
           "C++-set string tag is not resolvable from Python's oechem.",
)


class TestFchkRead:
    def test_next_returns_oemol(self, fchk_file):
        import oeio
        with oeio.read(fchk_file) as reader:
            mol = next(reader)
        assert isinstance(mol, oechem.OEMol)
        assert mol.NumAtoms() == 3

    def test_next_raises_stopiteration_at_eof(self, fchk_file):
        import oeio
        with oeio.read(fchk_file) as reader:
            next(reader)  # the single record
            with pytest.raises(StopIteration):
                next(reader)

    def test_iterate_yields_single_molecule(self, fchk_file):
        import oeio
        with oeio.read(fchk_file) as reader:
            mols = list(reader)
        assert len(mols) == 1

    def test_elements_and_title(self, fchk_file):
        import oeio
        with oeio.read(fchk_file) as reader:
            mol = next(reader)
        zs = sorted(a.GetAtomicNum() for a in mol.GetAtoms())
        assert zs == [1, 8, 9]  # H, O, F
        assert mol.GetTitle() == "Example"

    @requires_shared_openeye
    def test_typed_qm_scalars(self, fchk_file):
        import oeio
        with oeio.read(fchk_file) as reader:
            mol = next(reader)
        assert mol.HasData("Total Energy")
        energy = mol.GetData("Total Energy")
        assert isinstance(energy, float)
        assert energy == pytest.approx(-175.4626954507162, abs=1e-9)
        assert mol.GetData("Charge") == 0
        assert isinstance(mol.GetData("Charge"), int)
        assert mol.GetData("Multiplicity") == 1
        assert mol.GetData("Number of electrons") == 18
        assert mol.GetData("Basis") == "6-311+G(2d,p)"

    def test_min_fixture(self, fchk_min_file):
        import oeio
        with oeio.read(fchk_min_file) as reader:
            mol = next(reader)
        assert mol.NumAtoms() == 2

    @requires_shared_openeye
    def test_min_fixture_total_energy(self, fchk_min_file):
        import oeio
        with oeio.read(fchk_min_file) as reader:
            mol = next(reader)
        assert mol.GetData("Total Energy") == pytest.approx(-75.0, abs=1e-9)

    def test_malformed_raises_format_error(self, fchk_bad_z_file):
        import oeio
        with oeio.read(fchk_bad_z_file) as reader:
            with pytest.raises(oeio.FormatError):
                next(reader)
