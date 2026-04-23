"""Tests for oeio entry-point plugin discovery."""

import warnings
from unittest.mock import patch, MagicMock

import pytest


class TestLoadPlugins:
    """Test _load_plugins() entry-point discovery."""

    def test_calls_entry_points_with_correct_group(self):
        """_load_plugins queries the 'oeio.plugins' group."""
        from oeio import _load_plugins

        with patch("oeio.entry_points", create=True) as mock_ep:
            mock_ep.return_value = []
            with patch("importlib.metadata.entry_points", mock_ep):
                _load_plugins()
            mock_ep.assert_called_once_with(group="oeio.plugins")

    def test_broken_plugin_warns_but_does_not_raise(self):
        """A plugin that raises on import emits a warning, not an exception."""
        from oeio import _load_plugins

        bad_ep = MagicMock()
        bad_ep.name = "broken_plugin"
        bad_ep.load.side_effect = ImportError("no such module")

        with patch("importlib.metadata.entry_points", return_value=[bad_ep]):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                _load_plugins()
            assert len(w) == 1
            assert "broken_plugin" in str(w[0].message)
            assert issubclass(w[0].category, RuntimeWarning)

    def test_successful_plugin_loaded(self):
        """A valid entry point's load() method is called."""
        from oeio import _load_plugins

        good_ep = MagicMock()
        good_ep.name = "good_plugin"
        good_ep.load.return_value = MagicMock()

        with patch("importlib.metadata.entry_points", return_value=[good_ep]):
            _load_plugins()
        good_ep.load.assert_called_once()

    def test_formats_includes_oechem(self):
        """Sanity check: built-in OEChem handler is always present."""
        import oeio

        names = [f.name for f in oeio.formats()]
        assert "OEChem" in names
