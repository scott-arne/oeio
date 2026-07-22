"""
OEIO - Modular package for deploying new readers and writers for the OpenEye Toolkits
"""

import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shutil
import sys
import warnings
from importlib import metadata
from pathlib import Path

__version__ = "0.4.1"
__version_info__ = (0, 4, 1)


_OPENEYE_COMPAT_PRELOAD_PATHS: list[str] = []
_OPENEYE_COMPAT_EXTENSION_DIR: Path | None = None


def _user_cache_root():
    """Return the per-user cache root for OpenEye compatibility aliases."""
    cache_home = os.environ.get("XDG_CACHE_HOME")
    if cache_home:
        return Path(cache_home) / "oeio"
    return Path.home() / ".cache" / "oeio"


def _runtime_openeye_version():
    """Return the installed OpenEye toolkit distribution version if available."""
    try:
        return metadata.version("openeye-toolkits")
    except metadata.PackageNotFoundError:
        return "unknown"


def _cache_key(oe_lib_dir, expected_libs, build_version, runtime_version):
    """Build a stable cache key for one OpenEye runtime library set."""
    key_data = "\n".join(
        [
            os.path.realpath(oe_lib_dir),
            build_version or "unknown",
            runtime_version or "unknown",
            *sorted(expected_libs),
        ]
    )
    return hashlib.sha256(key_data.encode("utf-8")).hexdigest()[:16]


def _runtime_shared_library_names(lib_names):
    """Return filenames that can participate in runtime dynamic loading."""
    return [
        lib_name
        for lib_name in lib_names
        if ".so" in lib_name
        or lib_name.endswith(".dylib")
        or lib_name.endswith(".dll")
    ]


def _is_openeye_runtime_library_name(lib_name):
    """Return whether a dependency belongs to the OpenEye runtime set."""
    return lib_name.startswith("liboe") or lib_name.startswith("libzstd.")


def _find_openeye_runtime_lib_dir(expected_libs=()):
    """Find the OpenEye runtime library directory without importing oechem."""
    search_locations = []
    openeye_module = sys.modules.get("openeye")
    openeye_path = getattr(openeye_module, "__path__", None)
    if openeye_path is not None:
        search_locations.extend(openeye_path)

    if not search_locations:
        try:
            openeye_spec = importlib.util.find_spec("openeye")
        except (ImportError, ValueError):
            openeye_spec = None
        if (
            openeye_spec is not None
            and openeye_spec.submodule_search_locations is not None
        ):
            search_locations.extend(openeye_spec.submodule_search_locations)

    expected_libs = set(_runtime_shared_library_names(expected_libs or ()))
    fallback_dir = None
    for package_root in search_locations:
        libs_root = Path(package_root) / "libs"
        if not libs_root.is_dir():
            continue

        # Importing openeye.libs eagerly imports oechem in some environments.
        # The runtime libraries are shipped below openeye/libs, so filesystem
        # discovery preserves the fresh-import condition.
        for root, _, files in os.walk(libs_root):
            file_set = set(files)
            if expected_libs and expected_libs.intersection(file_set):
                return root
            if fallback_dir is None and any(
                ".dylib" in lib_name or ".so" in lib_name or ".dll" in lib_name
                for lib_name in files
            ):
                fallback_dir = root

    return fallback_dir


def _library_family(lib_name):
    """Return the stable library family name for a versioned shared library."""
    match = re.match(r"(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$", lib_name)
    if match is None:
        return None
    return match.group(1)


def _candidate_runtime_libraries(oe_lib_dir, expected_name):
    """Find runtime libraries with the same family as an expected filename."""
    family = _library_family(expected_name)
    if family is None:
        return []
    candidates = []
    for file_name in os.listdir(oe_lib_dir):
        candidate_path = os.path.join(oe_lib_dir, file_name)
        if not os.path.isfile(candidate_path):
            continue
        if file_name.startswith(f"{family}-") or file_name.startswith(f"{family}."):
            candidates.append(candidate_path)
    return sorted(candidates)


def _compatible_library_path(oe_lib_dir, expected_name):
    """Return a runtime library path and whether it needs an expected-name alias."""
    exact_path = os.path.join(oe_lib_dir, expected_name)
    if os.path.isfile(exact_path):
        return exact_path, False

    candidates = _candidate_runtime_libraries(oe_lib_dir, expected_name)
    if len(candidates) != 1:
        candidate_names = ", ".join(os.path.basename(path) for path in candidates)
        raise ImportError(
            f"Could not find a compatible OpenEye runtime library for "
            f"{expected_name!r} in {oe_lib_dir!r}. "
            f"Candidates: {candidate_names or 'none'}."
        )
    return candidates[0], True


def _extension_runtime_library_names(pkg_dir):
    """Return OpenEye runtime library names recorded by the extension."""
    extension_path = _find_extension_module_path(pkg_dir)
    if extension_path is None:
        return []

    if sys.platform == "darwin":
        return _mach_o_runtime_library_names(extension_path)
    if sys.platform.startswith("linux"):
        return _elf_runtime_library_names(extension_path)
    return []


def _mach_o_runtime_library_names(extension_path):
    """Return OpenEye dylib dependencies recorded in a Mach-O extension."""
    import subprocess

    try:
        result = subprocess.run(
            ["otool", "-L", str(extension_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return []

    dependencies = []
    for line in result.stdout.splitlines()[1:]:
        dependency = line.strip().split(" ", 1)[0]
        lib_name = os.path.basename(dependency)
        if _is_openeye_runtime_library_name(lib_name):
            dependencies.append(lib_name)
    return dependencies


def _elf_runtime_library_names(extension_path):
    """Return OpenEye shared-library dependencies recorded in an ELF extension."""
    import subprocess

    try:
        result = subprocess.run(
            ["readelf", "-d", str(extension_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return []

    dependencies = []
    for match in re.finditer(r"Shared library: \[(?P<name>[^\]]+)\]", result.stdout):
        lib_name = match.group("name")
        if _is_openeye_runtime_library_name(lib_name):
            dependencies.append(lib_name)
    return dependencies


def _ensure_cache_alias(cache_dir, expected_name, target_path):
    """Create or refresh an expected-name symlink in the user cache."""
    alias_path = cache_dir / expected_name
    if alias_path.is_symlink():
        if alias_path.resolve() == Path(target_path).resolve():
            return alias_path
        alias_path.unlink()
    elif alias_path.exists():
        raise ImportError(
            f"Cannot create OpenEye compatibility alias {alias_path}: "
            "a non-symlink file already exists at that path."
        )

    try:
        alias_path.symlink_to(target_path)
    except OSError as exc:
        raise ImportError(
            f"Could not create OpenEye compatibility alias "
            f"{alias_path} -> {target_path}: {exc}"
        ) from exc
    return alias_path


def _ensure_library_compat():
    """Prepare compatibility aliases when OpenEye library filenames drift.

    When oeio is built with shared OpenEye libraries, the compiled extension
    records the exact versioned library filenames (e.g., liboechem-4.3.0.1.dylib).
    If the user upgrades openeye-toolkits, these filenames change and the dynamic
    linker fails to load the extension.

    This function creates expected-name aliases in a user-writable cache instead
    of mutating the installed package directory. When aliases are needed, the
    extension is later loaded from the same cache directory so its $ORIGIN lookup
    can find those aliases.
    """
    global _OPENEYE_COMPAT_EXTENSION_DIR, _OPENEYE_COMPAT_PRELOAD_PATHS

    _OPENEYE_COMPAT_PRELOAD_PATHS = []
    _OPENEYE_COMPAT_EXTENSION_DIR = None

    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = set(_runtime_shared_library_names(
        getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    ))
    expected_libs.update(_extension_runtime_library_names(os.path.dirname(__file__)))
    expected_libs = sorted(expected_libs)
    if not expected_libs:
        return False

    oe_lib_dir = _find_openeye_runtime_lib_dir(expected_libs)
    if oe_lib_dir is None:
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    runtime_version = _runtime_openeye_version()
    cache_dir = (
        _user_cache_root()
        / "openeye-libs"
        / _cache_key(oe_lib_dir, expected_libs, build_version, runtime_version)
    )

    preload_paths = []
    needs_cached_origin = False
    for expected_name in expected_libs:
        actual_path, needs_alias = _compatible_library_path(oe_lib_dir, expected_name)
        if needs_alias:
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError as exc:
                raise ImportError(
                    f"Could not create OpenEye compatibility cache directory "
                    f"{cache_dir}: {exc}"
                ) from exc
            alias_path = _ensure_cache_alias(cache_dir, expected_name, actual_path)
            preload_paths.append(str(alias_path))
            needs_cached_origin = True
        else:
            preload_paths.append(actual_path)

    _OPENEYE_COMPAT_PRELOAD_PATHS = preload_paths
    if needs_cached_origin:
        _OPENEYE_COMPAT_EXTENSION_DIR = cache_dir

    return needs_cached_origin


def _extension_suffixes():
    """Return extension-module suffixes for the active Python interpreter."""
    return tuple(importlib.machinery.EXTENSION_SUFFIXES)


def _find_extension_module_path(pkg_dir):
    """Find the installed _oeio extension file."""
    for suffix in _extension_suffixes():
        candidate = Path(pkg_dir) / f"_oeio{suffix}"
        if candidate.is_file():
            return candidate
    for candidate in Path(pkg_dir).glob("_oeio*"):
        if candidate.is_file() and str(candidate).endswith(_extension_suffixes()):
            return candidate
    return None


def _copy_if_stale(source_path, target_path):
    """Copy a file into the cache when size or mtime changed."""
    if (
        target_path.exists()
        and target_path.stat().st_size == source_path.stat().st_size
        and target_path.stat().st_mtime_ns == source_path.stat().st_mtime_ns
    ):
        return
    shutil.copy2(source_path, target_path)


def _copy_package_shared_sidecars(pkg_dir, cache_dir, extension_path):
    """Copy package-local shared library sidecars needed by cached extension."""
    for candidate in Path(pkg_dir).iterdir():
        name = candidate.name
        if not candidate.is_file() or candidate == extension_path:
            continue
        if (
            ".so" not in name
            and not name.endswith(".dylib")
            and not name.endswith(".dll")
            and not name.endswith(".pyd")
        ):
            continue
        _copy_if_stale(candidate, cache_dir / name)


def _load_cached_extension_if_needed():
    """Load _oeio from the cache when OpenEye aliases live there."""
    cache_dir = _OPENEYE_COMPAT_EXTENSION_DIR
    if cache_dir is None:
        return

    module_name = f"{__name__}._oeio"
    if module_name in sys.modules:
        return

    pkg_dir = os.path.dirname(__file__)
    extension_path = _find_extension_module_path(pkg_dir)
    if extension_path is None:
        return

    cached_extension_path = cache_dir / extension_path.name
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        _copy_if_stale(extension_path, cached_extension_path)
        _copy_package_shared_sidecars(pkg_dir, cache_dir, extension_path)
    except OSError as exc:
        raise ImportError(
            f"Could not prepare cached oeio extension in {cache_dir}: {exc}"
        ) from exc

    spec = importlib.util.spec_from_file_location(module_name, cached_extension_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import spec for {cached_extension_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise


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

    expected_libs = _runtime_shared_library_names(
        getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    )
    if not expected_libs:
        return

    oe_lib_dir = _find_openeye_runtime_lib_dir(expected_libs)
    if oe_lib_dir is None:
        return

    if not os.path.isdir(oe_lib_dir):
        return

    paths = _OPENEYE_COMPAT_PRELOAD_PATHS
    if not paths:
        paths = [
            os.path.join(oe_lib_dir, lib_name)
            for lib_name in expected_libs
            if os.path.exists(os.path.join(oe_lib_dir, lib_name))
        ]

    for path in paths:
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
        from importlib import metadata
        runtime_version = metadata.version("openeye-toolkits")
    except metadata.PackageNotFoundError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "Install with: pip install openeye-toolkits",
            ImportWarning
        )
        return

    build_parts = build_version.split('.')[:2]
    runtime_parts = runtime_version.split('.')[:2]
    if build_parts != runtime_parts:
        warnings.warn(
            f"OpenEye version mismatch: oeio was built with "
            f"OpenEye Toolkits {build_version} but runtime has {runtime_version}. "
            f"This may cause compatibility issues.",
            RuntimeWarning
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
_load_cached_extension_if_needed()
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
    _mol_to_string,
    _mol_from_string,
    _mol_to_bytes,
    _mol_from_bytes,
    _mol_to_bytes_bridged,
    _mol_to_string_bridged,
    _scalar_generic_data,
    _resolve_format,
    _NEEDS_DATA_REATTACH,
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
    :param reader_factory: Callable(path) -> iterable of oechem molecules
        (any type; not required to be OEGraphMol), or None.
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
    :param reader_factory: Callable(path) -> iterable of oechem molecules
        (any type; the iterable is returned to the caller as-is and is not
        required to expose as_type/read_into), or None.
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

    .. note::
        The default ``OEMol`` iteration, :meth:`Reader.as_type`, and
        :meth:`Reader.read_into` apply to the built-in SWIG-backed
        :class:`Reader`. When a registered plugin handler matches ``path``,
        this returns the plugin's ``reader_factory`` result instead, which is
        governed by that plugin's own contract and may not expose ``as_type``/
        ``read_into`` or default to ``OEMol``.

    :param path: Path to a molecular file.
    :param config: Optional handler-specific configuration.
    :returns: A :class:`Reader` (built-in formats) or a plugin-specific iterable.
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


# ============================================================================
# In-memory molecule serialization
# ============================================================================


def _as_format_code(format):
    """Resolve a format (str token or OEFormat_* int) to an OEFormat int.

    :param format: A format token (``"oeb"``, ``".sdf"``) or an ``OEFormat_*``
        integer.
    :returns: The integer OEFormat code.
    :raises FormatError: If a string token is not a recognized OEChem format.
    """
    if isinstance(format, str):
        return _resolve_format(format)
    return int(format)


def _new_mol(mol_type):
    """Validate ``mol_type`` and return a fresh instance.

    :param mol_type: An ``oechem`` molecule class.
    :returns: A new ``mol_type`` instance.
    :raises TypeError: If ``mol_type`` is not an ``OEMolBase`` subclass.
    """
    from openeye import oechem

    if not (isinstance(mol_type, type) and issubclass(mol_type, oechem.OEMolBase)):
        raise TypeError(
            "mol_type must be an oechem OEMolBase subclass, got {!r}".format(mol_type))
    return mol_type()


def _reattach_after_from(mol):
    """Re-key C++-set scalar generic data under Python's registry (static builds).

    On a static build ``mol_from_*`` sets named scalar data under oeio's tag
    registry, which is invisible to Python by name. Re-key each item under
    Python's registry and delete the oeio-keyed copy (a move, not a copy, so no
    duplication). No-op on shared builds, where the two registries coincide.

    :param mol: The molecule just populated by a ``from_*`` call.
    """
    if not _NEEDS_DATA_REATTACH:
        return
    for name, value, tag in _scalar_generic_data(mol):
        if not mol.HasData(name):
            mol.SetData(name, value)
            mol.DeleteData(tag)


def _has_scalar_data(mol):
    """Return ``True`` if the molecule carries any molecule-level generic data.

    A cheap peek used to decide whether the static-build serialize bridge is
    needed; data-free molecules skip the bridge entirely.

    :param mol: The molecule to inspect.
    :returns: ``True`` if any generic-data item is present.
    """
    for _ in mol.GetDataIter():
        return True
    return False


def _scalar_pairs_for_bridge(mol):
    """Return parallel ``(names, values)`` lists of the molecule's scalar data.

    Names are enumerated through Python's ``oechem`` registry so they are
    correct at the caller side; the C++ bridge re-sets them under oeio's
    registry. Non-scalar values are out of contract and filtered by the bridge.

    :param mol: The molecule whose scalar generic data to collect.
    :returns: A ``(names, values)`` tuple of parallel lists.
    """
    from openeye import oechem

    names, values = [], []
    for dp in mol.GetDataIter():
        name = oechem.OEGetTag(dp.GetTag())
        if not name:
            continue
        if not mol.HasData(name):
            continue
        names.append(name)
        values.append(mol.GetData(name))
    return names, values


def to_bytes(mol, format="oeb", flavor=None, gzip=False):
    """Serialize a single molecule to binary ``bytes`` (default OEB).

    :param mol: The molecule to serialize (single- or multi-conformer).
    :param format: Format token or ``OEFormat_*`` int; defaults to ``"oeb"``.
    :param flavor: Output flavor; ``None`` uses the per-format default.
    :param gzip: Whether to gzip-compress the output.
    :returns: The serialized bytes.
    :raises FormatError: If the format is not writeable.
    """
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    if _NEEDS_DATA_REATTACH and _has_scalar_data(mol):
        names, values = _scalar_pairs_for_bridge(mol)
        return _mol_to_bytes_bridged(mol, names, values, fmt, fl, bool(gzip))
    return _mol_to_bytes(mol, fmt, fl, bool(gzip))


def from_bytes(data, format="oeb", flavor=None, gzip=False, mol_type=None):
    """Deserialize a single molecule from binary ``bytes`` (default OEB).

    :param data: The bytes to parse (first record only).
    :param format: Format token or ``OEFormat_*`` int; defaults to ``"oeb"``.
    :param flavor: Input flavor; ``None`` uses the per-format default.
    :param gzip: Whether the input is gzip-compressed.
    :param mol_type: Molecule class to return; defaults to ``oechem.OEMol``.
    :returns: A new molecule of ``mol_type``.
    :raises FormatError: If parsing fails or the format is not readable.
    :raises TypeError: If ``mol_type`` is not an ``OEMolBase`` subclass.
    """
    from openeye import oechem

    mol = _new_mol(oechem.OEMol if mol_type is None else mol_type)
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    if not _mol_from_bytes(mol, data, fmt, fl, bool(gzip)):
        raise FormatError("oeio: failed to parse molecule from bytes")
    _reattach_after_from(mol)
    return mol


def from_bytes_into(mol, data, format="oeb", flavor=None, gzip=False):
    """Deserialize into a caller-provided molecule (zero-allocation).

    :param mol: The molecule to populate.
    :returns: ``True`` if a molecule was read, ``False`` otherwise.
    """
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    ok = _mol_from_bytes(mol, data, fmt, fl, bool(gzip))
    if ok:
        _reattach_after_from(mol)
    return ok


def to_string(mol, format, flavor=None):
    """Serialize a single molecule to a text ``str`` (uncompressed).

    :param mol: The molecule to serialize.
    :param format: Format token or ``OEFormat_*`` int (required).
    :param flavor: Output flavor; ``None`` uses the per-format default.
    :returns: The serialized text.
    :raises FormatError: If the format is not writeable.
    """
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    if _NEEDS_DATA_REATTACH and _has_scalar_data(mol):
        names, values = _scalar_pairs_for_bridge(mol)
        return _mol_to_string_bridged(mol, names, values, fmt, fl)
    return _mol_to_string(mol, fmt, fl)


def from_string(data, format, flavor=None, mol_type=None):
    """Deserialize a single molecule from a text ``str``.

    :param data: The text to parse (first record only).
    :param format: Format token or ``OEFormat_*`` int (required).
    :param flavor: Input flavor; ``None`` uses the per-format default.
    :param mol_type: Molecule class to return; defaults to ``oechem.OEMol``.
    :returns: A new molecule of ``mol_type``.
    :raises FormatError: If parsing fails or the format is not readable.
    """
    from openeye import oechem

    mol = _new_mol(oechem.OEMol if mol_type is None else mol_type)
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    if not _mol_from_string(mol, data, fmt, fl):
        raise FormatError("oeio: failed to parse molecule from string")
    _reattach_after_from(mol)
    return mol


def from_string_into(mol, data, format, flavor=None):
    """Deserialize text into a caller-provided molecule.

    :returns: ``True`` if a molecule was read, ``False`` otherwise.
    """
    fmt = _as_format_code(format)
    fl = 0 if flavor is None else int(flavor)
    ok = _mol_from_string(mol, data, fmt, fl)
    if ok:
        _reattach_after_from(mol)
    return ok


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
    "to_bytes",
    "from_bytes",
    "from_bytes_into",
    "to_string",
    "from_string",
    "from_string_into",
    "Reader",
    "Error",
    "FormatError",
    "FileError",
    "FormatInfo",
    "ReaderConfig",
    "WriterConfig",
]
