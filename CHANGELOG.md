# Changelog

All notable changes to `oeio` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.1] - 2026-07-09

### Changed

- **Faster CUBE reads.** The volumetric block is now parsed with a
  line-buffered `strtod` loop instead of `std::istream::operator>>`, measured
  ~6x faster end-to-end on large grids. Output is bit-identical and peak memory
  stays bounded (parsing proceeds line-by-line rather than buffering the whole
  volumetric tail).
- **Lower-memory FCHK reads.** The reader now streams the file with a
  single-line cursor instead of loading every line into memory. For large FCHKs
  dominated by unconsumed arrays (MO coefficients, densities) this cuts peak
  memory from O(file) to O(1) — ~2.4x lower resident memory on a large fixture —
  and runs modestly faster. Parsing behavior and error reporting are unchanged.

## [0.4.0] - 2026-07-06

### Added

- **Gaussian CUBE file support** (`.cube`, `.cub`) — a
  C++-first built-in handler that reads/writes a molecule plus N volumetric
  `OEScalarGrid` datasets (N>1 for multi-orbital "MO cubes").
  - `reader.with_grids()` iterates `(molecule, grids)` records.
  - `reader.read_into(mol, *grids)` fills caller-owned grids (returns N, or
    None at EOF); `reader.read_into(mol)` reads the molecule only.
  - `writer.append(mol, *grids)` writes a CUBE (at least one grid required).
  - New grid overloads on the C++ `MolSource`/`MolSink`/`MolRange`/`Writer`
    interfaces (additive; non-grid formats unchanged).
- `Reader.__next__`, so `next(reader)` returns the next molecule.

### Limitations

- CUBE grids must be axis-aligned with uniform (cubic-voxel) spacing;
  rotated/skewed/anisotropic grids raise `FormatError`. Atom coordinates and
  grid geometry are converted to Ångström on read, Bohr on write.

## [0.3.0] - 2026-07-04

### Changed (Breaking)

- **The default molecule type when reading is now `OEMol` (multi-conformer)
  instead of `OEGraphMol` (single-conformer)**, in both C++ and Python. This
  aligns with OpenEye Orion's default molecule type. Bare iteration now yields
  `OEMol`:

  ```python
  for mol in oeio.read("molecules.oeb"):   # mol is now oechem.OEMol
      ...
  ```
  ```cpp
  for (auto& mol : oeio::read("molecules.oeb")) { /* OEMol */ }
  ```

  `OEMol` is **not** a subclass of `OEGraphMol`, so `isinstance`/`dynamic_cast`
  checks that assumed `OEGraphMol` will behave differently. Multi-conformer
  records (e.g. from `.oeb`) are now assembled into a single `OEMol` rather
  than flattened into one single-conformer molecule per conformer.

  **Migration — to restore the previous single-conformer behavior**, request
  the type explicitly:

  ```python
  with oeio.read("ligands.sdf") as reader:
      for mol in reader.as_type(oechem.OEGraphMol):   # OEGraphMol
          ...
  ```
  ```cpp
  for (auto& mol : oeio::read("ligands.sdf").as<OEChem::OEGraphMol>()) { /* OEGraphMol */ }
  ```

- **The C++ `oeio::transform` callback signature widened from
  `void(OEChem::OEGraphMol&)` to `void(OEChem::OEMolBase&)`.** Existing
  `transform` callbacks typed as `OEGraphMol&` must be updated to `OEMolBase&`
  (`SetTitle`, `GetTitle`, `NumAtoms`, etc. are all available on `OEMolBase`):

  ```cpp
  oeio::transform(range, [](OEChem::OEMolBase& mol) { /* ... */ });
  ```

### Added

- **Typed iteration.** Choose the molecule type for iteration without changing
  the default:
  - C++: `oeio::read(path).as<T>()` returns a `TypedMolRange<T>` yielding `T&`
    (any type convertible to `OEChem::OEMolBase&`, e.g. `OEMol`, `OEGraphMol`,
    `OEQMol`).
  - Python: `reader.as_type(cls)` yields fresh `cls` instances.
- **Zero-copy by-reference reading** via `read_into(mol)` (C++
  `MolRange::read_into(OEChem::OEMolBase&)` and Python `Reader.read_into(mol)`),
  which populates a caller-owned molecule without allocating a new one per
  record. Useful for large files.
- **Python `filter`/`transform` now preserve the input molecule type**
  (they yield `type(mol)(mol)` rather than always producing `OEGraphMol`), so an
  `OEMol` stream stays `OEMol` — and multi-conformer records survive filtering
  and transformation.

### Fixed

- **Multi-conformer conformers are now preserved on the write path.** Writing an
  `OEMol` (e.g. `oeio::read("mc.oeb") | oeio::write(...)`) previously flattened
  it to a single conformer on disk because the write dispatched on the static
  `OEMolBase&` overload. The handler now dispatches on the dynamic type so
  multi-conformer molecules write all their conformers, while single-conformer
  molecules retain molecule-level title/data.
- **Reading into an `OEMol`/`OEQMol` now dispatches to the correct
  `OEReadMolecule` overload** (previously the single-conformer `OEMolBase&`
  overload was always selected, flattening multi-conformer records).
- `Reader.__iter__`, `Reader.as_type`, and `Reader.read_into` raise `ValueError`
  (rather than an internal `AttributeError`) if the reader is closed mid-iteration.

## [0.2.5]

- Prior release. See git history for details.
