# oeio

A modular molecule I/O library for the [OpenEye Toolkits](https://www.eyesopen.com/).
oeio provides range-based reading, RAII writing, and pipeline operators for
molecular file formats, with a plugin architecture that lets you register new
format handlers within your own libraries.

Both C++ and Python interfaces are provided (via SWIG)

## Supported Formats

The built-in OEChem handler supports all OEChem-native formats:

| Format | Extensions |
|--------|------------|
| SDF | `.sdf`, `.sdf.gz` |
| MOL | `.mol` |
| MOL2 | `.mol2`, `.mol2.gz` |
| PDB | `.pdb`, `.pdb.gz`, `.ent`, `.ent.gz` |
| OEB | `.oeb`, `.oeb.gz`, `.oez` |
| SMILES | `.smi`, `.smi.gz`, `.ism`, `.can` |
| CSV | `.csv`, `.csv.gz` |
| XYZ | `.xyz`, `.xyz.gz` |
| Other | `.fasta`, `.cif`, `.mmcif`, `.mopac`, `.cdx` |

Additional formats can be added through the plugin system (see
[Registering Format Plugins](#registering-format-plugins)).

## Installation

### From a Wheel

```bash
pip install oeio
```

Requires `openeye-toolkits` to be installed separately.

### From Source

Prerequisites:

- OpenEye C++ SDK
- CMake >= 3.16
- SWIG >= 4.0
- Python >= 3.10

Configure and build:

```bash
cmake --preset debug
cmake --build build-debug
```

Install the Python package for development:

```bash
pip install --config-settings editable_mode=compat -e python/
```

## Usage

### Python

**Reading molecules:**

```python
import oeio

for mol in oeio.read("molecules.sdf"):
    print(mol.GetTitle(), mol.NumAtoms())
```

**Writing molecules:**

```python
import oeio

with oeio.write("output.sdf") as writer:
    for mol in oeio.read("input.sdf"):
        writer.append(mol)
```

**Filtering:**

```python
import oeio

large_mols = oeio.filter(
    oeio.read("input.sdf"),
    lambda mol: mol.NumAtoms() > 10
)
for mol in large_mols:
    print(mol.GetTitle())
```

**Transforming:**

```python
import oeio
from openeye import oechem

prepared = oeio.transform(
    oeio.read("input.sdf"),
    lambda mol: oechem.OEAddExplicitHydrogens(mol)
)
with oeio.write("prepared.sdf") as writer:
    for mol in prepared:
        writer.append(mol)
```

**Listing registered formats:**

```python
for fmt in oeio.formats():
    print(f"{fmt.name}: {', '.join(fmt.extensions)}")
```

### C++

Include the umbrella header:

```cpp
#include <oeio/oeio.h>
```

**Reading molecules:**

```cpp
for (auto& mol : oeio::read("molecules.sdf")) {
    std::cout << mol.GetTitle() << "\n";
}
```

**Pipeline: filter, transform, write:**

```cpp
oeio::read("input.sdf")
    | oeio::filter([](const OEChem::OEMolBase& mol) {
        return mol.NumAtoms() > 10;
    })
    | oeio::transform([](OEChem::OEGraphMol& mol) {
        OEChem::OEAddExplicitHydrogens(mol);
    })
    | oeio::write("output.sdf");
```

The pipe operator (`|`) connects a `MolRange` to a `Writer`, streaming molecules
without intermediate storage.

**Reading with configuration:**

```cpp
oeio::oechem::ReaderConfig cfg;
cfg.iflavor_format = OEChem::OEFormat::SDF;
cfg.iflavor = OEChem::OEIFlavor::SDF::Default;

for (auto& mol : oeio::read("input.sdf", std::any(cfg))) {
    // ...
}
```

## Configuration

Both `ReaderConfig` and `WriterConfig` allow fine-grained control over format
behavior:

| Field | Type | Description |
|-------|------|-------------|
| `format` | `unsigned int` | OEFormat constant (0 = auto-detect from extension) |
| `iflavor_format` / `oflavor_format` | `unsigned int` | Format for flavor setting |
| `iflavor` / `oflavor` | `unsigned int` | OEIFlavor/OEOFlavor value |
| `num_threads` | `unsigned int` | Thread count (0 = single-threaded) |

In most cases the defaults work -- format is detected from the file extension.

## Running Tests

C++ tests:

```bash
cmake --build build-debug
cd build-debug && ctest --output-on-failure
```

Python tests:

```bash
pytest tests/python/ -v
```

## Registering Format Plugins

oeio uses a plugin architecture that allows you to add support for new molecular
file formats. A plugin consists of three parts: a `MolSource` for reading, a
`MolSink` for writing, and a `FormatHandler` that ties them together.

### Step 1: Implement MolSource and MolSink

```cpp
#include <oeio/format_handler.h>

namespace myformat {

class MyReader : public oeio::MolSource {
public:
    MyReader(const std::string& path) {
        // Open file, initialize parser
    }

    bool next(OEChem::OEGraphMol& mol) override {
        // Read next molecule into mol.
        // Return true on success, false at end-of-file.
    }

private:
    // File handle, parser state, etc.
};

class MyWriter : public oeio::MolSink {
public:
    MyWriter(const std::string& path) {
        // Open file for writing
    }

    bool write(const OEChem::OEMolBase& mol) override {
        // Write molecule. Return true on success.
    }

    void close() override {
        // Flush and close.
    }
};

} // namespace myformat
```

### Step 2: Implement FormatHandler

The handler provides metadata and factory methods for creating readers and writers:

```cpp
#include <oeio/format_handler.h>

namespace myformat {

class MyHandler : public oeio::FormatHandler {
public:
    oeio::FormatInfo info() const override {
        return {
            "MyFormat",                      // name
            {".myf", ".myf.gz"},             // extensions
            "My custom molecular format",    // description
            true,                            // supports_read
            true,                            // supports_write
            false,                           // supports_threaded_read
            false                            // supports_threaded_write
        };
    }

    std::unique_ptr<oeio::MolSource> make_reader(
            const std::string& path,
            const std::any& config) const override {
        return std::make_unique<MyReader>(path);
    }

    std::unique_ptr<oeio::MolSink> make_writer(
            const std::string& path,
            const std::any& config) const override {
        return std::make_unique<MyWriter>(path);
    }
};

} // namespace myformat
```

### Step 3: Register the Handler

Use the `OEIO_REGISTER_FORMAT` macro at file scope to register the handler when
the library loads:

```cpp
#include <oeio/format_registry.h>

OEIO_REGISTER_FORMAT(myformat::MyHandler)
```

Once registered, all oeio functions automatically recognize your format:

```cpp
// These just work -- oeio looks up ".myf" in the registry
for (auto& mol : oeio::read("data.myf")) { ... }
oeio::read("data.myf") | oeio::write("output.myf");
```

The `config` parameter in `make_reader` and `make_writer` carries format-specific
options as a type-erased `std::any`. You can define your own config struct and
extract it with `std::any_cast`.

## Building Wheels

```bash
python scripts/build_python.py --openeye-root /path/to/openeye/sdk
```

Options: `--clean`, `--upload`, `--test-upload`, `--verbose`. If `--openeye-root`
is not provided, the script checks `OPENEYE_ROOT`, then `CMakePresets.json`.

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `OEIO_BUILD_TESTS` | ON | Build C++ tests |
| `OEIO_BUILD_PYTHON` | ON | Build Python SWIG bindings |
| `OEIO_BUILD_BENCHMARKS` | OFF | Build performance benchmarks |
| `OEIO_UNIVERSAL2` | OFF | Build macOS universal2 binary |
| `OEIO_USE_STABLE_ABI` | ON | Use Python stable ABI (abi3) |

## License

MIT
