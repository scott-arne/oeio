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
    """Molecule data must survive round-trips at full fidelity on all builds.

    Both engines are full-fidelity: on shared builds oeio's C++ is the engine;
    on static builds serialization goes directly through OpenEye's Python
    ``oechem`` (one tag registry). Single-molecule OEB/SDF reads drop
    molecule-level generic/SD data into an ``OEMol`` but preserve it into an
    ``OEGraphMol`` (inherent OpenEye single-mol-into-OEMol behavior), so the
    single-molecule assertions read back into ``OEGraphMol``. The multi-conformer
    case uses an ``OEMol`` source and target, which preserves molecule-level data.
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
        assert oechem.OEHasSDData(back, "k")
        assert oechem.OEGetSDData(back, "k") == "v"

    def test_sdf_string_sd_tag_round_trip(self):
        # Full-fidelity SD tag survives to_string/from_string (sdf).
        m = _ethanol()
        oechem.OESetSDData(m, "prop", "hello")
        back = oeio.from_string(oeio.to_string(m, "sdf"), "sdf",
                                mol_type=oechem.OEGraphMol)
        assert oechem.OEHasSDData(back, "prop")
        assert oechem.OEGetSDData(back, "prop") == "hello"

    def test_oeb_bytes_sd_tag_round_trip(self):
        # SD tags survive to_bytes/from_bytes (oeb) as well.
        m = _ethanol()
        oechem.OESetSDData(m, "prop", "hello")
        back = oeio.from_bytes(oeio.to_bytes(m), mol_type=oechem.OEGraphMol)
        assert oechem.OEHasSDData(back, "prop")
        assert oechem.OEGetSDData(back, "prop") == "hello"

    def test_no_duplicate_data_after_from_bytes(self):
        m = _ethanol()
        m.SetData("energy", 1.5)
        back = oeio.from_bytes(oeio.to_bytes(m), mol_type=oechem.OEGraphMol)
        # Exactly one generic-data item named "energy" and no stray items.
        names = [oechem.OEGetTag(dp.GetTag()) for dp in back.GetDataIter()]
        assert names.count("energy") == 1

    def test_multiconformer_with_data_round_trip(self):
        # Conformers AND molecule scalar data must both survive an OEB round-trip.
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
