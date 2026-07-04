# Changelog

All notable changes to `oeio` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
