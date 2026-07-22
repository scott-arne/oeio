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
