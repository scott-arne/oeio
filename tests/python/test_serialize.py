"""Tests for oeio in-memory bytes/string molecule serialization (Capability A)."""

import pytest
from openeye import oechem

import oeio


def _ethanol():
    m = oechem.OEGraphMol()
    oechem.OESmilesToMol(m, "CCO")
    m.SetTitle("ethanol")
    return m


class TestBytes:
    def test_to_bytes_returns_bytes(self):
        b = oeio.to_bytes(_ethanol())  # default OEB
        assert isinstance(b, bytes)
        assert len(b) > 0

    def test_oeb_bytes_round_trip(self):
        m = _ethanol()
        b = oeio.to_bytes(m)
        back = oeio.from_bytes(b)
        assert isinstance(back, oechem.OEMol)  # default mol_type
        assert back.NumAtoms() == 3
        assert back.GetTitle() == "ethanol"

    def test_from_bytes_mol_type(self):
        b = oeio.to_bytes(_ethanol())
        back = oeio.from_bytes(b, mol_type=oechem.OEGraphMol)
        assert isinstance(back, oechem.OEGraphMol)

    def test_from_bytes_into(self):
        b = oeio.to_bytes(_ethanol())
        dst = oechem.OEGraphMol()
        assert oeio.from_bytes_into(dst, b) is True
        assert dst.NumAtoms() == 3

    def test_gzip_sdf_bytes_round_trip(self):
        m = _ethanol()
        b = oeio.to_bytes(m, format="sdf", gzip=True)
        assert isinstance(b, bytes)
        back = oeio.from_bytes(b, format="sdf", gzip=True)
        assert back.NumAtoms() == 3

    def test_format_int_accepted(self):
        b = oeio.to_bytes(_ethanol(), format=oechem.OEFormat_OEB)
        assert isinstance(b, bytes)


class TestString:
    def test_smiles_round_trip(self):
        s = oeio.to_string(_ethanol(), "smi")
        assert isinstance(s, str)
        back = oeio.from_string(s, "smi")
        assert back.NumAtoms() == 3

    def test_sdf_sd_tag_round_trip(self):
        m = _ethanol()
        oechem.OESetSDData(m, "prop", "42")
        s = oeio.to_string(m, "sdf")
        back = oeio.from_string(s, "sdf", mol_type=oechem.OEGraphMol)
        assert oechem.OEHasSDData(back, "prop")
        assert oechem.OEGetSDData(back, "prop") == "42"


class TestGenericDataPreservation:
    """SD/generic-data must survive round-trips on both shared and static builds.

    Single-molecule OEB/SDF reads drop molecule-level generic/SD data when the
    target is an ``OEMol`` but preserve it into an ``OEGraphMol`` (inherent
    OpenEye behavior, unrelated to the registry bridge), so the single-molecule
    assertions read back into ``OEGraphMol``. The multi-conformer case uses an
    ``OEMol`` source and target, which preserves molecule-level data.
    """

    def test_oeb_bytes_preserves_generic_data(self):
        m = _ethanol()
        m.SetData("energy", -175.46)
        m.SetData("label", "abc")
        m.SetData("count", 7)
        b = oeio.to_bytes(m)             # OEB
        back = oeio.from_bytes(b, mol_type=oechem.OEGraphMol)
        assert back.GetData("energy") == pytest.approx(-175.46)
        assert back.GetData("label") == "abc"
        assert back.GetData("count") == 7

    def test_sdf_string_preserves_generic_data(self):
        m = _ethanol()
        oechem.OESetSDData(m, "k", "v")
        s = oeio.to_string(m, "sdf")
        back = oeio.from_string(s, "sdf", mol_type=oechem.OEGraphMol)
        assert oechem.OEGetSDData(back, "k") == "v"

    def test_no_duplicate_data_after_from_bytes(self):
        m = _ethanol()
        m.SetData("energy", 1.5)
        back = oeio.from_bytes(oeio.to_bytes(m), mol_type=oechem.OEGraphMol)
        # Exactly one generic-data item named "energy" and no stray items.
        names = [oechem.OEGetTag(dp.GetTag()) for dp in back.GetDataIter()]
        assert names.count("energy") == 1

    def test_multiconformer_with_data_round_trip(self):
        # Exercises the type-preserving bridged temp: conformers AND molecule
        # scalar data must both survive (the serialize bridge must not flatten).
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "CCO")
        oechem.OEAddExplicitHydrogens(mol)
        src = oechem.OEMol(mol)
        for conf in src.GetConfs():
            mol.NewConf(conf)
        assert mol.NumConfs() > 1
        mol.SetData("energy", -1.25)
        back = oeio.from_bytes(oeio.to_bytes(mol))   # OEB
        assert back.NumConfs() == mol.NumConfs()
        assert back.GetData("energy") == pytest.approx(-1.25)


class TestBridgeDirect:
    """Exercise the low-level serialize bridge directly (all builds).

    On a shared build ``_NEEDS_DATA_REATTACH`` is False, so ``to_bytes``/
    ``to_string`` never route through the bridge; nothing else covers the
    type-preserving temp. These call the bridged C++ helpers directly, building
    ``names``/``values`` exactly as ``_scalar_pairs_for_bridge`` would, so the
    bridge is exercised on every build.
    """

    def test_bridged_bytes_preserves_title_and_data(self):
        # Regression: a single-conformer OEB written from an OEMol temp drops the
        # molecule title; the bridge must use an OEGraphMol temp so it survives.
        m = _ethanol()  # titled single-conformer OEGraphMol
        m.SetData("energy", -175.46)
        names, values = oeio._scalar_pairs_for_bridge(m)
        b = oeio._mol_to_bytes_bridged(m, names, values, oechem.OEFormat_OEB, 0, False)
        back = oechem.OEGraphMol()
        assert oeio._mol_from_bytes(back, b, oechem.OEFormat_OEB, 0, False)
        assert back.GetTitle() == "ethanol"
        assert back.GetData("energy") == pytest.approx(-175.46)

    def test_bridged_bytes_multiconformer_preserves_confs_and_data(self):
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, "CCO")
        oechem.OEAddExplicitHydrogens(mol)
        src = oechem.OEMol(mol)
        for conf in src.GetConfs():
            mol.NewConf(conf)
        assert mol.NumConfs() > 1
        mol.SetData("energy", -1.25)
        names, values = oeio._scalar_pairs_for_bridge(mol)
        b = oeio._mol_to_bytes_bridged(mol, names, values, oechem.OEFormat_OEB, 0, False)
        back = oechem.OEMol()
        assert oeio._mol_from_bytes(back, b, oechem.OEFormat_OEB, 0, False)
        assert back.NumConfs() == mol.NumConfs()
        assert back.GetData("energy") == pytest.approx(-1.25)

    def test_bridged_string_sdf_preserves_title(self):
        m = _ethanol()
        m.SetData("energy", -175.46)
        names, values = oeio._scalar_pairs_for_bridge(m)
        s = oeio._mol_to_string_bridged(m, names, values, oechem.OEFormat_SDF, 0)
        assert isinstance(s, str)
        back = oechem.OEGraphMol()
        assert oeio._mol_from_string(back, s, oechem.OEFormat_SDF, 0)
        assert back.GetTitle() == "ethanol"

    def test_bridged_length_mismatch_raises_and_leaves_mol_unchanged(self):
        # Mismatched name/value lengths must fail closed (raise, not crash) and
        # leave the caller's molecule untouched (only the throwaway temp is built).
        m = _ethanol()
        m.SetData("energy", -1.0)
        before_title = m.GetTitle()
        before_energy = m.GetData("energy")
        with pytest.raises(Exception):
            oeio._mol_to_bytes_bridged(m, ["a", "b"], [1], oechem.OEFormat_OEB, 0, False)
        assert m.GetTitle() == before_title
        assert m.GetData("energy") == before_energy
        names = [oechem.OEGetTag(dp.GetTag()) for dp in m.GetDataIter()]
        assert names.count("energy") == 1


class TestErrors:
    def test_unknown_format_raises(self):
        with pytest.raises(oeio.FormatError):
            oeio.to_string(_ethanol(), "not_a_format")

    def test_cube_format_raises(self):
        with pytest.raises(oeio.FormatError):
            oeio.to_string(_ethanol(), "cube")

    def test_bad_mol_type_raises(self):
        b = oeio.to_bytes(_ethanol())
        with pytest.raises(TypeError):
            oeio.from_bytes(b, mol_type=str)

    def test_malformed_string_raises(self):
        with pytest.raises(oeio.FormatError):
            oeio.from_string("this is not an sdf", "sdf")
