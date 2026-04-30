"""
OEIO - Modular package for deploying new readers and writers for the OpenEye Toolkits
"""

import os
import re
import warnings

__version__ = "0.2.2"
__version_info__ = (0, 2, 2)


def _ensure_library_compat():
    """Create compatibility symlinks when OpenEye library versions differ from build time.

    When this package is built with shared OpenEye libraries, the compiled extension
    records the exact versioned library filenames (e.g., liboechem-4.3.0.1.dylib).
    If the user upgrades openeye-toolkits, these filenames change and the dynamic
    linker fails to load the extension.

    This function detects version mismatches and creates symlinks from the expected
    (build-time) library names to the actual (runtime) library files.
    """
    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return False

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    pkg_dir = os.path.dirname(__file__)
    created_any = False

    for expected_name in expected_libs:
        if os.path.exists(os.path.join(oe_lib_dir, expected_name)):
            continue

        symlink_path = os.path.join(pkg_dir, expected_name)
        if os.path.islink(symlink_path):
            if os.path.exists(symlink_path):
                continue
            try:
                os.unlink(symlink_path)
            except OSError:
                continue
        elif os.path.exists(symlink_path):
            continue

        match = re.match(r'(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$', expected_name)
        if not match:
            continue
        base_name = match.group(1)

        actual_path = None
        for f in os.listdir(oe_lib_dir):
            if f.startswith(base_name + '-') or f.startswith(base_name + '.'):
                actual_path = os.path.join(oe_lib_dir, f)
                break

        if actual_path:
            try:
                os.symlink(actual_path, os.path.join(pkg_dir, expected_name))
                created_any = True
            except OSError:
                pass

    return created_any


def _preload_shared_libs():
    """Preload OpenEye shared libraries so the C extension can find them.

    On Linux, the extension's RUNPATH (set at build time) normally handles
    dependency resolution, but preloading ensures libraries are available
    even if RUNPATH is stripped (e.g. by certain packaging tools).
    On macOS, @rpath references may not resolve without preloading.

    Only the libraries recorded in ``OPENEYE_EXPECTED_LIBS`` are loaded,
    and they are loaded with ``RTLD_GLOBAL`` so that cross-module C++
    symbol references resolve correctly. Loading the entire OpenEye
    library directory (which can contain 70+ unrelated shared objects)
    would pollute the global symbol namespace and cause segfaults in
    unrelated C extensions such as ``_sqlite3``.
    """
    import ctypes
    import sys
    if sys.platform not in ('linux', 'darwin'):
        return

    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return

    if not os.path.isdir(oe_lib_dir):
        return

    pkg_dir = os.path.dirname(__file__)
    for lib_name in expected_libs:
        oe_path = os.path.join(oe_lib_dir, lib_name)
        local_path = os.path.join(pkg_dir, lib_name)
        path = oe_path if os.path.exists(oe_path) else local_path
        if os.path.exists(path) or os.path.islink(path):
            try:
                ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass


def _preload_bundled_libs():
    """Preload libraries bundled by auditwheel from the .libs directory.

    auditwheel repair bundles non-manylinux dependencies (e.g. libraries
    from FetchContent or system packages) into a ``<package>.libs/``
    directory next to the package. The bundled copies have hashed filenames
    and must be loaded before the C extension to satisfy its DT_NEEDED
    entries.

    Libraries may have inter-dependencies, so we do multiple passes
    until no new libraries can be loaded.
    """
    import sys
    if sys.platform != 'linux':
        return

    import ctypes
    pkg_name = __name__
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    site_dir = os.path.dirname(pkg_dir)
    for libs_name in (f'{pkg_name}.libs', f'.{pkg_name}.libs'):
        libs_dir = os.path.join(site_dir, libs_name)
        if not os.path.isdir(libs_dir):
            continue
        remaining = [
            os.path.join(libs_dir, f)
            for f in sorted(os.listdir(libs_dir))
            if '.so' in f
        ]
        # Multi-pass: keep retrying until no progress (handles dep ordering)
        while remaining:
            failed = []
            for lib_path in remaining:
                try:
                    ctypes.CDLL(lib_path)
                except OSError:
                    failed.append(lib_path)
            if len(failed) == len(remaining):
                break  # No progress, stop
            remaining = failed


def _check_openeye_version():
    """Check that the OpenEye version matches what was used at build time."""
    try:
        from . import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    if not build_version:
        return

    try:
        from openeye import oechem
        runtime_version = oechem.OEToolkitsGetRelease()
        if runtime_version and build_version:
            build_parts = build_version.split('.')[:2]
            runtime_parts = runtime_version.split('.')[:2]
            if build_parts != runtime_parts:
                warnings.warn(
                    f"OpenEye version mismatch: oeio was built with "
                    f"OpenEye Toolkits {build_version} but runtime has {runtime_version}. "
                    f"This may cause compatibility issues.",
                    RuntimeWarning
                )
    except ImportError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "Install with: pip install openeye-toolkits",
            ImportWarning
        )


def _load_plugins():
    """Discover and load oeio format handler plugins via entry points.

    Scans for packages that declare an ``oeio.plugins`` entry point group.
    Each entry point should point to a module whose import triggers
    OEIO_REGISTER_FORMAT registration at the C++ level.

    A broken or missing plugin emits a warning but does not prevent oeio
    from loading.
    """
    from importlib.metadata import entry_points

    for ep in entry_points(group="oeio.plugins"):
        try:
            ep.load()
        except Exception as exc:
            warnings.warn(
                f"oeio: failed to load plugin '{ep.name}': {exc}",
                RuntimeWarning,
                stacklevel=2,
            )


_ensure_library_compat()
_preload_shared_libs()
_preload_bundled_libs()
_check_openeye_version()

from .oeio import (
    read as _cpp_read,
    write as _cpp_write,
    formats as _cpp_formats,
    filter,
    transform,
    Reader,
    Error,
    FormatError,
    FileError,
    FormatInfo,
    FormatRegistry,
    ReaderConfig,
    WriterConfig,
    _ReaderHandle,
    _WriterHandle,
    _open_reader,
    _open_writer,
)


# ============================================================================
# Python-level plugin registry
# ============================================================================
# Plugins that are installed as separate packages (e.g. oemaestro) each get
# their own copy of the oeio C++ FormatRegistry when statically linked.
# The Python-level registry bridges this gap: plugins register their
# reader/writer factories here, and the top-level read()/write()/formats()
# functions check this registry alongside the C++ one.

class _PluginHandler:
    """Descriptor for a Python-registered format handler.

    :param name: Human-readable format name (e.g. "Maestro").
    :param extensions: List of file extensions (e.g. [".mae", ".mae.gz"]).
    :param description: Short description of the format.
    :param reader_factory: Callable(path) -> iterator of OEGraphMol, or None.
    :param writer_factory: Callable(path) -> context manager with add(mol), or None.
    """

    def __init__(self, name, extensions, description="",
                 reader_factory=None, writer_factory=None):
        self.name = name
        self.extensions = list(extensions)
        self.description = description
        self.reader_factory = reader_factory
        self.writer_factory = writer_factory

    def matches(self, path):
        """Return True if the path ends with one of this handler's extensions.

        Longest extensions are checked first so that compound extensions
        like ``.mae.gz`` match before ``.gz``.
        """
        lower = path.lower()
        for ext in sorted(self.extensions, key=len, reverse=True):
            if lower.endswith(ext):
                return True
        return False

    def to_format_info(self):
        """Return a FormatInfo-compatible object for this handler."""
        return _PyFormatInfo(
            name=self.name,
            extensions=self.extensions,
            description=self.description,
            supports_read=self.reader_factory is not None,
            supports_write=self.writer_factory is not None,
            supports_threaded_read=False,
            supports_threaded_write=False,
        )


class _PyFormatInfo:
    """Pure-Python FormatInfo compatible with the SWIG FormatInfo struct."""

    def __init__(self, name, extensions, description,
                 supports_read, supports_write,
                 supports_threaded_read, supports_threaded_write):
        self.name = name
        self.extensions = extensions
        self.description = description
        self.supports_read = supports_read
        self.supports_write = supports_write
        self.supports_threaded_read = supports_threaded_read
        self.supports_threaded_write = supports_threaded_write

    def __repr__(self):
        return (f"FormatInfo(name={self.name!r}, "
                f"extensions={self.extensions!r})")


class _PluginRegistry:
    """Python-level format handler registry for cross-package plugins."""

    def __init__(self):
        self._handlers = []

    def register(self, handler):
        """Register a _PluginHandler."""
        self._handlers.append(handler)

    def lookup(self, path):
        """Find a handler matching the given path, or None."""
        for h in self._handlers:
            if h.matches(path):
                return h
        return None

    def formats(self):
        """Return FormatInfo objects for all registered Python handlers."""
        return [h.to_format_info() for h in self._handlers]


_plugin_registry = _PluginRegistry()


def register_handler(name, extensions, description="",
                     reader_factory=None, writer_factory=None):
    """Register a Python-level format handler with oeio.

    This is used by external plugins (e.g. oemaestro) to register their
    format handlers so that ``oeio.read()`` and ``oeio.write()`` can
    dispatch to them.

    :param name: Human-readable format name (e.g. "Maestro").
    :param extensions: List of file extensions (e.g. [".mae", ".mae.gz"]).
    :param description: Short description of the format.
    :param reader_factory: Callable(path) -> iterable of OEGraphMol, or None.
    :param writer_factory: Callable(path) -> context manager with add(mol), or None.
    """
    _plugin_registry.register(_PluginHandler(
        name=name,
        extensions=extensions,
        description=description,
        reader_factory=reader_factory,
        writer_factory=writer_factory,
    ))


def read(path, config=None):
    """Open a molecule reader for ``path``.

    Checks Python-registered plugin handlers first, then falls back to the
    C++ FormatRegistry.

    :param path: Path to a molecular file.
    :param config: Optional handler-specific configuration.
    :returns: A :class:`Reader` or plugin-specific iterable.
    """
    path = str(path)
    handler = _plugin_registry.lookup(path)
    if handler and handler.reader_factory:
        return handler.reader_factory(path)
    return _cpp_read(path, config)


def write(path, config=None):
    """Open a molecule writer for ``path``.

    Checks Python-registered plugin handlers first, then falls back to the
    C++ FormatRegistry.

    :param path: Path to write to.
    :param config: Optional handler-specific configuration.
    :returns: A context manager with an ``add(mol)`` method.
    """
    path = str(path)
    handler = _plugin_registry.lookup(path)
    if handler and handler.writer_factory:
        return handler.writer_factory(path)
    return _cpp_write(path, config)


def formats():
    """List all registered molecular file formats.

    Includes both C++ and Python-registered format handlers.

    :returns: List of FormatInfo objects.
    """
    cpp_fmts = _cpp_formats()
    py_fmts = _plugin_registry.formats()
    return cpp_fmts + py_fmts


_load_plugins()

__all__ = [
    "__version__",
    "__version_info__",
    "read",
    "write",
    "formats",
    "filter",
    "transform",
    "register_handler",
    "Reader",
    "Error",
    "FormatError",
    "FileError",
    "FormatInfo",
    "ReaderConfig",
    "WriterConfig",
]
