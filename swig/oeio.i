// swig/oeio.i
// SWIG interface file for oeio Python bindings
%module _oeio

%{
#include "oeio/oeio.h"
#include <oechem.h>
#include <oegrid.h>
#include <oesystem.h>

using namespace oeio;
%}

// ============================================================================
// Forward declarations for cross-module SWIG type resolution
// ============================================================================
// These enable typemaps for OpenEye types whose definitions live in the
// OpenEye SWIG runtime (v4). Only types you actually use in your wrapped API
// need full #include — forward declarations suffice for the typemaps.

namespace OEChem {
    class OEMolBase;
    class OEMCMolBase;
    class OEMol;
    class OEGraphMol;
    class OEAtomBase;
    class OEBondBase;
    class OEConfBase;
    class OEMatchBase;
    class OEMolDatabase;
    class oemolistream;
    class oemolostream;
    class OEQMol;
    class OEResidue;
    class OEUniMolecularRxn;
}

namespace OEBio {
    class OEDesignUnit;
    class OEHierView;
    class OEHierResidue;
    class OEHierFragment;
    class OEHierChain;
    class OEInteractionHint;
    class OEInteractionHintContainer;
}

namespace OEDocking {
    class OEReceptor;
}

namespace OEPlatform {
    class oeifstream;
    class oeofstream;
    class oeisstream;
    class oeosstream;
}

namespace OESystem {
    class OEScalarGrid;
    class OERecord;
    class OEMolRecord;
}

// ============================================================================
// Cross-runtime SWIG compatibility layer
// ============================================================================
// OpenEye's Python bindings use SWIG runtime v4; our module uses v5.
// Since the runtimes are separate, SWIG_TypeQuery cannot access OpenEye types.
// We use Python isinstance for type safety and directly extract the void*
// pointer from the SwigPyObject struct layout (stable across SWIG versions).
//
// This approach enables passing OpenEye objects between Python and C++ without
// serialization. The macros below generate the boilerplate for each type.

%{
// Minimal SwigPyObject layout compatible across SWIG runtime versions.
// The actual struct may have more fields, but ptr is always first after
// PyObject_HEAD.
struct _SwigPyObjectCompat {
    PyObject_HEAD
    void *ptr;
};

static void* _oeio_extract_swig_ptr(PyObject* obj) {
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        PyErr_Clear();
        return NULL;
    }
    void* ptr = ((_SwigPyObjectCompat*)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return ptr;
}

// ---- Type checker generator macro ----
// Generates a cached isinstance checker for an OpenEye Python type.
// TAG:    identifier suffix (e.g., oemolbase)
// MODULE: Python module string (e.g., "openeye.oechem")
// CLASS:  Python class name string (e.g., "OEMolBase")
#define DEFINE_OE_TYPE_CHECKER(TAG, MODULE, CLASS) \
    static PyObject* _oeio_oe_##TAG##_type = NULL; \
    static bool _oeio_is_##TAG(PyObject* obj) { \
        if (!_oeio_oe_##TAG##_type) { \
            PyObject* mod = PyImport_ImportModule(MODULE); \
            if (mod) { \
                _oeio_oe_##TAG##_type = PyObject_GetAttrString(mod, CLASS); \
                Py_DECREF(mod); \
            } \
            if (!_oeio_oe_##TAG##_type) return false; \
        } \
        return PyObject_IsInstance(obj, _oeio_oe_##TAG##_type) == 1; \
    }

// ---- Molecule types (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolbase,    "openeye.oechem", "OEMolBase")
DEFINE_OE_TYPE_CHECKER(oemcmolbase,  "openeye.oechem", "OEMCMolBase")
DEFINE_OE_TYPE_CHECKER(oemol,        "openeye.oechem", "OEMol")
DEFINE_OE_TYPE_CHECKER(oegraphmol,   "openeye.oechem", "OEGraphMol")
DEFINE_OE_TYPE_CHECKER(oeqmol,       "openeye.oechem", "OEQMol")

// ---- Atom / bond / conformer / residue (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeatombase,   "openeye.oechem", "OEAtomBase")
DEFINE_OE_TYPE_CHECKER(oebondbase,   "openeye.oechem", "OEBondBase")
DEFINE_OE_TYPE_CHECKER(oeconfbase,   "openeye.oechem", "OEConfBase")
DEFINE_OE_TYPE_CHECKER(oeresidue,    "openeye.oechem", "OEResidue")
DEFINE_OE_TYPE_CHECKER(oematchbase,  "openeye.oechem", "OEMatchBase")

// ---- Molecule I/O (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolistream, "openeye.oechem", "oemolistream")
DEFINE_OE_TYPE_CHECKER(oemolostream, "openeye.oechem", "oemolostream")
DEFINE_OE_TYPE_CHECKER(oemoldatabase,"openeye.oechem", "OEMolDatabase")

// ---- Reactions (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeunimolecularrxn, "openeye.oechem", "OEUniMolecularRxn")

// ---- Platform streams (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeifstream,   "openeye.oechem", "oeifstream")
DEFINE_OE_TYPE_CHECKER(oeofstream,   "openeye.oechem", "oeofstream")
DEFINE_OE_TYPE_CHECKER(oeisstream,   "openeye.oechem", "oeisstream")
DEFINE_OE_TYPE_CHECKER(oeosstream,   "openeye.oechem", "oeosstream")

// ---- Records (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oerecord,     "openeye.oechem", "OERecord")
DEFINE_OE_TYPE_CHECKER(oemolrecord,  "openeye.oechem", "OEMolRecord")

// ---- Bio / hierarchy (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oedesignunit, "openeye.oechem", "OEDesignUnit")
DEFINE_OE_TYPE_CHECKER(oehierview,   "openeye.oechem", "OEHierView")
DEFINE_OE_TYPE_CHECKER(oehierresidue,"openeye.oechem", "OEHierResidue")
DEFINE_OE_TYPE_CHECKER(oehierfragment,"openeye.oechem","OEHierFragment")
DEFINE_OE_TYPE_CHECKER(oehierchain,  "openeye.oechem", "OEHierChain")
DEFINE_OE_TYPE_CHECKER(oeinteractionhint,          "openeye.oechem", "OEInteractionHint")
DEFINE_OE_TYPE_CHECKER(oeinteractionhintcontainer, "openeye.oechem", "OEInteractionHintContainer")

// ---- Grid (openeye.oegrid) ----
DEFINE_OE_TYPE_CHECKER(oescalargrid, "openeye.oegrid", "OEScalarGrid")

// ---- Docking (openeye.oedocking) ----
DEFINE_OE_TYPE_CHECKER(oereceptor,   "openeye.oedocking", "OEReceptor")

#undef DEFINE_OE_TYPE_CHECKER

// ---- OEScalarGrid return-type helper (zero-copy pointer swap) ----
static PyObject* _oeio_wrap_as_oe_grid(OESystem::OEScalarGrid* grid) {
    if (!grid) {
        Py_RETURN_NONE;
    }
    PyObject* oegrid_mod = PyImport_ImportModule("openeye.oegrid");
    if (!oegrid_mod) {
        delete grid;
        return NULL;
    }
    PyObject* grid_cls = PyObject_GetAttrString(oegrid_mod, "OEScalarGrid");
    Py_DECREF(oegrid_mod);
    if (!grid_cls) {
        delete grid;
        return NULL;
    }
    PyObject* oe_grid = PyObject_CallNoArgs(grid_cls);
    Py_DECREF(grid_cls);
    if (!oe_grid) {
        delete grid;
        return NULL;
    }
    PyObject* thisAttr = PyObject_GetAttrString(oe_grid, "this");
    if (!thisAttr) {
        PyErr_Clear();
        Py_DECREF(oe_grid);
        delete grid;
        return NULL;
    }
    _SwigPyObjectCompat* swig_this = (_SwigPyObjectCompat*)thisAttr;
    delete reinterpret_cast<OESystem::OEScalarGrid*>(swig_this->ptr);
    swig_this->ptr = grid;
    Py_DECREF(thisAttr);
    return oe_grid;
}
%}

// ============================================================================
// Typemap generator macros
// ============================================================================

// Generate const-ref and non-const-ref typemaps for a cross-runtime OpenEye type.
// CPP_TYPE: fully qualified C++ type (e.g., OEChem::OEMolBase)
// CHECKER:  isinstance checker function name
// ERR_MSG:  error message on type mismatch
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oeio_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oeio_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            if (CHECKER($input)) {
                argp = _oeio_extract_swig_ptr($input);
                if (argp) res = SWIG_OK;
            }
        }
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const CPP_TYPE* {
    if ($input == Py_None) {
        $1 = 1;
    } else {
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, 0);
        $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
    }
}

%enddef

// ============================================================================
// Typemap declarations for all OpenEye types
// ============================================================================
// Each type gets const-ref and non-const-ref typemaps. Types that commonly
// appear as optional parameters also get nullable-pointer typemaps.
// These are inert until a wrapped function signature uses the type.

// ---- Molecule hierarchy (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolBase,    _oeio_is_oemolbase,    "Expected OEMolBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMCMolBase,  _oeio_is_oemcmolbase,  "Expected OEMCMolBase-derived object (OEMCMolBase or OEMol).")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMol,        _oeio_is_oemol,        "Expected OEMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEGraphMol,   _oeio_is_oegraphmol,   "Expected OEGraphMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEQMol,       _oeio_is_oeqmol,       "Expected OEQMol object.")

// ---- Atom / bond / conformer / residue / match (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEAtomBase,   _oeio_is_oeatombase,   "Expected OEAtomBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEBondBase,   _oeio_is_oebondbase,   "Expected OEBondBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEConfBase,   _oeio_is_oeconfbase,   "Expected OEConfBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEResidue,    _oeio_is_oeresidue,    "Expected OEResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMatchBase,  _oeio_is_oematchbase,  "Expected OEMatchBase-derived object.")

// ---- Molecule I/O (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolistream,  _oeio_is_oemolistream, "Expected oemolistream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolostream,  _oeio_is_oemolostream, "Expected oemolostream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolDatabase, _oeio_is_oemoldatabase,"Expected OEMolDatabase object.")

// ---- Reactions (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEUniMolecularRxn, _oeio_is_oeunimolecularrxn, "Expected OEUniMolecularRxn object.")

// ---- Platform streams (OEPlatform) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeifstream, _oeio_is_oeifstream, "Expected oeifstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeofstream, _oeio_is_oeofstream, "Expected oeofstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeisstream, _oeio_is_oeisstream, "Expected oeisstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeosstream, _oeio_is_oeosstream, "Expected oeosstream object.")

// ---- Records (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OERecord,    _oeio_is_oerecord,    "Expected OERecord object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEMolRecord, _oeio_is_oemolrecord, "Expected OEMolRecord object.")

// ---- Bio / hierarchy (OEBio) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEDesignUnit,   _oeio_is_oedesignunit, "Expected OEDesignUnit object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierView,     _oeio_is_oehierview,   "Expected OEHierView object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierResidue,  _oeio_is_oehierresidue,"Expected OEHierResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierFragment,  _oeio_is_oehierfragment,"Expected OEHierFragment object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierChain,    _oeio_is_oehierchain,  "Expected OEHierChain object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHint,          _oeio_is_oeinteractionhint,          "Expected OEInteractionHint object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHintContainer, _oeio_is_oeinteractionhintcontainer, "Expected OEInteractionHintContainer object.")

// ---- Grid (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEScalarGrid, _oeio_is_oescalargrid, "Expected OEScalarGrid-derived object.")
OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(OESystem::OEScalarGrid, _oeio_is_oescalargrid, "Expected OEScalarGrid or None.")

// OEScalarGrid return-type typemap (wraps C++ grid as native openeye.oegrid object)
%typemap(out) OESystem::OEScalarGrid* {
    $result = _oeio_wrap_as_oe_grid($1);
    if (!$result) SWIG_fail;
}

// ---- Docking (OEDocking) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _oeio_is_oereceptor, "Expected OEReceptor object.")

// ============================================================================
// Include STL typemaps
// ============================================================================
%include "std_string.i"
%include "stdint.i"
%include "exception.i"

// ============================================================================
// Exception handling
// ============================================================================
// Python exception objects populated at module-load time via %init below.
// The %exception block raises the most-specific Python subclass that matches
// the C++ type, falling back to RuntimeError for foreign std::exceptions.
%{
static PyObject* _oeio_py_Error        = nullptr;
static PyObject* _oeio_py_FormatError  = nullptr;
static PyObject* _oeio_py_FileError    = nullptr;

/// Install the Python exception classes so %exception can raise them.
/// Called once from %pythoncode at module-import time.
static void _oeio_set_exception_types(PyObject* err,
                                      PyObject* fmt_err,
                                      PyObject* file_err) {
    Py_XINCREF(err);      _oeio_py_Error       = err;
    Py_XINCREF(fmt_err);  _oeio_py_FormatError = fmt_err;
    Py_XINCREF(file_err); _oeio_py_FileError   = file_err;
}
%}

%inline %{
void _install_exception_types(PyObject* err,
                              PyObject* fmt_err,
                              PyObject* file_err) {
    _oeio_set_exception_types(err, fmt_err, file_err);
}
%}

%exception {
    try {
        $action
    } catch (const oeio::FormatError& e) {
        PyErr_SetString(_oeio_py_FormatError, e.what());
        SWIG_fail;
    } catch (const oeio::FileError& e) {
        PyErr_SetString(_oeio_py_FileError, e.what());
        SWIG_fail;
    } catch (const oeio::Error& e) {
        PyErr_SetString(_oeio_py_Error, e.what());
        SWIG_fail;
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception");
    }
}

// ============================================================================
// Additional STL typemaps
// ============================================================================
%include "std_vector.i"

// ============================================================================
// Version macros
// ============================================================================
#define OEIO_VERSION_MAJOR 0
#define OEIO_VERSION_MINOR 2
#define OEIO_VERSION_PATCH 2

// ============================================================================
// Header includes for SWIG compilation
// ============================================================================
%{
#include "oeio/oeio.h"
#include "oeio/format_handler.h"
#include "oeio/format_registry.h"
#include "oeio/oechem_config.h"
using namespace oeio;

// Force the linker to include oechem_handler.o from liboeio.a
// so that the OEIO_REGISTER_FORMAT static initializer runs.
namespace oeio { void oeio_force_link_oechem_handler(); }
static struct _OeioForceLink {
    _OeioForceLink() { oeio::oeio_force_link_oechem_handler(); }
} _oeio_force_link;
%}

// ============================================================================
// FormatInfo struct
// ============================================================================
namespace oeio {

struct FormatInfo {
    std::string name;
    std::vector<std::string> extensions;
    std::string description;
    bool supports_read;
    bool supports_write;
    bool supports_threaded_read;
    bool supports_threaded_write;
};

}  // namespace oeio

// Template instantiations for STL containers
%template(FormatInfoVector) std::vector<oeio::FormatInfo>;
%template(StringVector) std::vector<std::string>;

// ============================================================================
// FormatRegistry (limited API for Python)
// ============================================================================
namespace oeio {

class FormatRegistry {
public:
    static FormatRegistry& instance();
    std::vector<FormatInfo> formats() const;

private:
    FormatRegistry();
};

}  // namespace oeio

// ============================================================================
// OEChem config structs
// ============================================================================
namespace oeio {
namespace oechem {

struct ReaderConfig {
    unsigned int format;
    unsigned int iflavor_format;
    unsigned int iflavor;
    unsigned int num_threads;
};

struct WriterConfig {
    unsigned int format;
    unsigned int oflavor_format;
    unsigned int oflavor;
    unsigned int num_threads;
};

}  // namespace oechem
}  // namespace oeio

// ============================================================================
// Reader and Writer handles (Python-specific wrappers)
// ============================================================================
%newobject oeio::_open_reader;
%newobject oeio::_open_writer;

// Declare the classes for SWIG without constructors (Python never calls them directly)
// Note: next() takes OEMolBase& (not OEGraphMol&) because the cross-SWIG-runtime
// pointer extracted from Python openeye.oechem.OEGraphMol objects points to an
// internal OEMolBase implementation, NOT to an OEGraphMol handle.
namespace oeio {
class _ReaderHandle {
public:
    bool next(OEChem::OEMolBase& mol);
private:
    _ReaderHandle();
};

class _WriterHandle {
public:
    bool append(const OEChem::OEMolBase& mol);
    void close();
    ~_WriterHandle();
private:
    _WriterHandle();
};
_ReaderHandle* _open_reader(const std::string& path);
_WriterHandle* _open_writer(const std::string& path);

}  // namespace oeio

// Full C++ definitions and factory functions (compiled but not parsed by SWIG)
%{
namespace oeio {

class _ReaderHandle {
public:
    _ReaderHandle(std::unique_ptr<MolSource> source)
        : source_(std::move(source)) {}

    /// Read the next molecule into an OEMolBase reference.
    ///
    /// Dispatches to MolSource::next(OEMolBase&) which is overridden by
    /// OEChemMolSource for zero-copy reading directly into the Python
    /// molecule's OEMolBase implementation. The source is responsible
    /// for clearing the molecule before reading.
    bool next(OEChem::OEMolBase& mol) {
        if (!source_) return false;
        return source_->next(mol);
    }

private:
    std::unique_ptr<MolSource> source_;
};

class _WriterHandle {
public:
    _WriterHandle(std::unique_ptr<MolSink> sink)
        : sink_(std::move(sink)) {}

    bool append(const OEChem::OEMolBase& mol) {
        return sink_ ? sink_->write(mol) : false;
    }

    void close() {
        if (sink_) {
            sink_->close();
            sink_.reset();
        }
    }

    ~_WriterHandle() {
        close();
    }

private:
    std::unique_ptr<MolSink> sink_;
};

/// Create a reader handle for a file path.
_ReaderHandle* _open_reader(const std::string& path) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto source = handler->make_reader(path, std::any{});
    if (!source) {
        throw FileError(
            "oeio: failed to create reader for '" + path + "'");
    }
    return new _ReaderHandle(std::move(source));
}

/// Create a writer handle for a file path.
_WriterHandle* _open_writer(const std::string& path) {
    auto* handler = FormatRegistry::instance().lookup(path);
    if (!handler) {
        throw FormatError(
            "oeio: unrecognized file extension for '" + path + "'");
    }
    auto sink = handler->make_writer(path, std::any{});
    if (!sink) {
        throw FileError(
            "oeio: failed to create writer for '" + path + "'");
    }
    return new _WriterHandle(std::move(sink));
}

}  // namespace oeio
%}

// ============================================================================
// Context manager protocol for _WriterHandle
// ============================================================================
%extend oeio::_WriterHandle {
%pythoncode %{
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
%}
}

// ============================================================================
// Module-level Python convenience API
// ============================================================================
%pythoncode %{
class Error(RuntimeError):
    """Base class for all oeio exceptions."""


class FormatError(Error):
    """Raised when a file extension or format hint is not registered."""


class FileError(Error):
    """Raised when a file cannot be opened or a reader/writer cannot be created."""


class Reader:
    """Iterable, closeable molecule reader.

    Returned by :func:`oeio.read`. Can be iterated directly or used as a
    context manager for deterministic cleanup::

        with oeio.read("mols.sdf") as reader:
            for mol in reader:
                ...

    :ivar _handle: Underlying ``_ReaderHandle`` or ``None`` after close.
    :ivar _closed: ``True`` once :meth:`close` has been called.
    """

    def __init__(self, handle):
        self._handle = handle
        self._closed = False

    def __iter__(self):
        from openeye import oechem

        if self._closed:
            raise ValueError("I/O operation on closed reader")
        mol = oechem.OEGraphMol()
        while self._handle.next(mol):
            yield mol
            mol = oechem.OEGraphMol()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        """Release the underlying reader. Idempotent."""
        if not self._closed:
            self._handle = None
            self._closed = True


def read(path, config=None):
    """Open a molecule reader for ``path``.

    :param path: Path to a molecular file.
    :param config: Optional handler-specific configuration.
    :returns: A :class:`Reader` that is both iterable and a context manager.

    Example::

        with oeio.read("input.sdf") as reader:
            for mol in reader:
                print(mol.GetTitle())
    """
    handle = _open_reader(str(path))
    return Reader(handle)


def write(path, config=None):
    """Create a molecule writer for a file path.

    Returns a context manager that wraps a _WriterHandle.

    :param path: Path to write to.
    :param config: Optional handler-specific configuration.
    :returns: _WriterHandle context manager.

    Example::

        with oeio.write("output.sdf") as writer:
            for mol in oeio.read("input.sdf"):
                writer.append(mol)
    """
    return _open_writer(str(path))


def filter(iterable, predicate):
    """Filter molecules from an iterable by a predicate.

    :param iterable: An iterable of molecules (e.g., from read()).
    :param predicate: Function that takes an OEMolBase and returns bool.
    :returns: Generator yielding molecules that satisfy the predicate.

    Example::

        heavy = oeio.filter(oeio.read("in.sdf"),
            lambda mol: mol.NumAtoms() > 10)
    """
    from openeye import oechem

    for mol in iterable:
        if predicate(mol):
            yield oechem.OEGraphMol(mol)


def transform(iterable, func):
    """Transform molecules from an iterable in-place.

    :param iterable: An iterable of molecules (e.g., from read()).
    :param func: Function that takes an OEGraphMol and modifies it in-place.
    :returns: Generator yielding transformed molecules.

    Example::

        prepared = oeio.transform(oeio.read("in.sdf"),
            lambda mol: oechem.OEAddExplicitHydrogens(mol))
    """
    from openeye import oechem

    for mol in iterable:
        func(mol)
        yield oechem.OEGraphMol(mol)


def formats():
    """List all registered molecular file formats.

    :returns: List of FormatInfo objects.
    """
    return list(FormatRegistry.instance().formats())


__version__ = "0.1.0"

_install_exception_types(Error, FormatError, FileError)
%}
