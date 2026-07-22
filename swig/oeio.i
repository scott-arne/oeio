// swig/oeio.i
// SWIG interface file for oeio Python bindings
%module _oeio

%{
#include <limits>

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

// ---- Grid-list input typemaps (cross-runtime) ----
// Convert a Python sequence of native openeye.oegrid.OEScalarGrid objects into a
// std::vector of (const) OEScalarGrid* by extracting each element's underlying
// pointer with the same proven primitives the single-grid typemaps use. We do
// NOT use SWIG %template vectors here: SWIG has no native knowledge of the
// cross-runtime OEScalarGrid type, and the OEScalarGrid* out-typemap transfers
// ownership (wrong for a borrowed vector element). Building the vector by hand
// keeps ownership with Python and lets the reader fill the caller's grids
// in-place (fill-through).
%typemap(in) const std::vector<OESystem::OEScalarGrid*>&
        (std::vector<OESystem::OEScalarGrid*> tmp, PyObject* keepalive = NULL) {
    // Materialize the sequence into a real tuple (keepalive) so every element
    // stays alive for the duration of the C++ call. Extracting a raw C++ pointer
    // and DECREF-ing the element first would be unsafe for a sequence whose
    // __getitem__ returns a fresh owning wrapper each call (the only reference
    // would drop and the pointer would dangle). PySequence_Tuple also rejects
    // non-sequences and reports a size error for us.
    keepalive = PySequence_Tuple($input);
    if (!keepalive) SWIG_exception_fail(SWIG_TypeError, "Expected a sequence of OEScalarGrid.");
    Py_ssize_t n = PyTuple_GET_SIZE(keepalive);
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PyTuple_GET_ITEM(keepalive, i);  // borrowed (keepalive owns it)
        if (!_oeio_is_oescalargrid(item)) {
            Py_DECREF(keepalive); keepalive = NULL;
            SWIG_exception_fail(SWIG_TypeError, "Expected OEScalarGrid in sequence.");
        }
        void* p = _oeio_extract_swig_ptr(item);
        if (!p) {
            Py_DECREF(keepalive); keepalive = NULL;
            SWIG_exception_fail(SWIG_RuntimeError, "failed to extract grid pointer");
        }
        tmp.push_back(reinterpret_cast<OESystem::OEScalarGrid*>(p));
    }
    $1 = &tmp;
}
// Free the keepalive tuple after the wrapped call has consumed the pointers.
// Null-out on the in-typemap error paths above makes this a safe no-op there.
%typemap(freearg) const std::vector<OESystem::OEScalarGrid*>& {
    Py_XDECREF(keepalive$argnum);
}

%typemap(in) const std::vector<const OESystem::OEScalarGrid*>&
        (std::vector<const OESystem::OEScalarGrid*> tmp, PyObject* keepalive = NULL) {
    // See the non-const overload above: hold a tuple reference so no element is
    // freed before the C++ call reads its extracted pointer.
    keepalive = PySequence_Tuple($input);
    if (!keepalive) SWIG_exception_fail(SWIG_TypeError, "Expected a sequence of OEScalarGrid.");
    Py_ssize_t n = PyTuple_GET_SIZE(keepalive);
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PyTuple_GET_ITEM(keepalive, i);  // borrowed (keepalive owns it)
        if (!_oeio_is_oescalargrid(item)) {
            Py_DECREF(keepalive); keepalive = NULL;
            SWIG_exception_fail(SWIG_TypeError, "Expected OEScalarGrid in sequence.");
        }
        void* p = _oeio_extract_swig_ptr(item);
        if (!p) {
            Py_DECREF(keepalive); keepalive = NULL;
            SWIG_exception_fail(SWIG_RuntimeError, "failed to extract grid pointer");
        }
        tmp.push_back(reinterpret_cast<const OESystem::OEScalarGrid*>(p));
    }
    $1 = &tmp;
}
%typemap(freearg) const std::vector<const OESystem::OEScalarGrid*>& {
    Py_XDECREF(keepalive$argnum);
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
#define OEIO_VERSION_MINOR 5
#define OEIO_VERSION_PATCH 0

// ============================================================================
// Header includes for SWIG compilation
// ============================================================================
%{
#include "oeio/oeio.h"
#include "oeio/format_handler.h"
#include "oeio/format_registry.h"
#include "oeio/oechem_config.h"
using namespace oeio;

// Force the linker to include oechem_handler.o and cube_handler.o from
// liboeio.a so that their OEIO_REGISTER_FORMAT static initializers run.
namespace oeio { void oeio_force_link_oechem_handler(); void oeio_force_link_cube_handler(); void oeio_force_link_fchk_handler(); }
static struct _OeioForceLink {
    _OeioForceLink() {
        oeio::oeio_force_link_oechem_handler();
        oeio::oeio_force_link_cube_handler();
        oeio::oeio_force_link_fchk_handler();
    }
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
// In-memory serialization (Capability A)
// ============================================================================
%{
#include "oeio/serialize.h"
%}

// Binary-safe payload type: a std::string carrying raw bytes. Distinct typemaps
// convert it to/from Python `bytes`, never `str`. These attach only to the
// specially-named oeio::Bytes typedef, so the default std::string<->str mapping
// used elsewhere in this module is untouched.
%typemap(out) oeio::Bytes {
    $result = PyBytes_FromStringAndSize($1.data(), (Py_ssize_t)$1.size());
}
%typemap(in) const oeio::Bytes& (oeio::Bytes tmp) {
    char* buf = nullptr; Py_ssize_t len = 0;
    if (PyBytes_AsStringAndSize($input, &buf, &len) != 0) SWIG_fail;
    // The copy runs in the input typemap, before the %exception try/catch wraps
    // the wrapped call, so a std::bad_alloc/length_error from a huge payload
    // would otherwise escape uncaught through the Python C API and abort. Convert
    // it to a Python exception here instead.
    try {
        tmp.assign(buf, (size_t)len);
    } catch (const std::bad_alloc&) {
        PyErr_NoMemory();
        SWIG_fail;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    }
    $1 = &tmp;
}
// The `in` typemap above uses a stack-local buffer, so nothing needs freeing.
// Override the default std_string freearg (reached via the Bytes->std::string
// typedef) which would reference an undeclared res$argnum for this type.
%typemap(freearg) const oeio::Bytes& "";

// Compile-only helpers (not SWIG-wrapped, like the static helpers above) for the
// static-build serialize bridge. Defined ahead of the %inline bridge entry
// points that call them.
%{
namespace oeio {
// Clear a temp's copied (stale-registry) generic data, then set the caller-
// supplied scalar pairs under THIS module's OpenEye tag registry (by name), so
// the subsequent serialize keys names through the same registry that will read
// them back. Operates on the abstract OEMolBase so it is independent of the
// temp's concrete type. `names`/`values` are parallel Python sequences.
//
// Fails closed: the pairs are validated and extracted into C++ holders BEFORE
// the temp is mutated, so a malformed call (non-sequence, length mismatch,
// non-string name, or a conversion/overflow failure) clears any pending Python
// error, releases every owned reference, and throws oeio::Error — which the
// %exception block maps to Python. The wrapped bridge functions therefore never
// serialize or emit partial output after a C-API error, and never return with a
// pending Python exception. PySequence_Size/GetItem are used (not the
// PySequence_Fast_GET_* macros, which are excluded under the Python limited/
// stable ABI this module targets).
static void _set_bridged_data(OEChem::OEMolBase& temp, PyObject* names,
                              PyObject* values) {
    Py_ssize_t nn = PySequence_Size(names);
    Py_ssize_t nv = PySequence_Size(values);
    if (nn < 0 || nv < 0 || nn != nv) {
        PyErr_Clear();
        throw oeio::Error("oeio: bridged data name/value length mismatch");
    }

    // Default-initialized so an unused alternative is never read.
    struct Pair {
        std::string name;
        int kind = -1;  // -1 skip, 0 bool, 1 int, 2 double, 3 string
        bool b = false;
        int i = 0;
        double d = 0.0;
        std::string s;
    };
    std::vector<Pair> pairs;
    pairs.reserve(static_cast<size_t>(nn));
    for (Py_ssize_t idx = 0; idx < nn; ++idx) {
        PyObject* nm = PySequence_GetItem(names, idx);   // new ref (or NULL)
        PyObject* vv = PySequence_GetItem(values, idx);  // new ref (or NULL)
        // Release owned refs, clear the pending PyErr, and throw so the wrapped
        // call never proceeds to serialize.
        auto fail = [&](const char* msg) {
            Py_XDECREF(nm); Py_XDECREF(vv);
            PyErr_Clear();
            throw oeio::Error(msg);
        };
        if (!nm || !vv) fail("oeio: bridged data element access failed");
        if (!PyUnicode_Check(nm)) fail("oeio: bridged data name is not a string");
        // SWIG_PyUnicode_AsUTF8AndSize (not PyUnicode_AsUTF8) is used so this
        // compiles under the limited ABI; the returned char* is valid only while
        // its companion bytes object lives, so copy before releasing it.
        PyObject* nmbytes = nullptr;
        const char* cname = SWIG_PyUnicode_AsUTF8AndSize(nm, nullptr, &nmbytes);
        if (!cname) { Py_XDECREF(nmbytes); fail("oeio: bridged data name is not valid UTF-8"); }
        Pair p;
        p.name.assign(cname);
        Py_XDECREF(nmbytes);
        // PyBool_Check precedes PyLong_Check because bool is a subclass of int.
        if (PyBool_Check(vv)) {
            p.kind = 0; p.b = (vv == Py_True);
        } else if (PyLong_Check(vv)) {
            long lv = PyLong_AsLong(vv);
            if (lv == -1 && PyErr_Occurred()) fail("oeio: bridged data int value out of range");
            // The oeio generic-data int model is int32 (SetData/GetIntData);
            // fail closed rather than silently truncate a value that fits C long
            // but not int (reachable only via an out-of-contract direct call).
            if (lv < std::numeric_limits<int>::min() ||
                lv > std::numeric_limits<int>::max())
                fail("oeio: bridged int value out of 32-bit range");
            p.kind = 1; p.i = static_cast<int>(lv);
        } else if (PyFloat_Check(vv)) {
            double dv = PyFloat_AsDouble(vv);
            if (dv == -1.0 && PyErr_Occurred()) fail("oeio: bridged data float value is invalid");
            p.kind = 2; p.d = dv;
        } else if (PyUnicode_Check(vv)) {
            PyObject* vbytes = nullptr;
            const char* sval = SWIG_PyUnicode_AsUTF8AndSize(vv, nullptr, &vbytes);
            if (!sval) { Py_XDECREF(vbytes); fail("oeio: bridged data string value is not valid UTF-8"); }
            p.kind = 3; p.s.assign(sval); Py_XDECREF(vbytes);
        }
        // Non-scalar values leave p.kind == -1 (out of contract) and are skipped.
        Py_XDECREF(nm); Py_XDECREF(vv);
        if (p.kind >= 0) pairs.push_back(std::move(p));
    }

    // All Python work validated; mutate the temp now (no Python errors possible
    // past this point). Clear stale data (collect tags first; cannot delete while
    // iterating), then set the validated pairs.
    std::vector<unsigned int> stale;
    for (OESystem::OEIter<OESystem::OEBaseData> it = temp.GetDataIter(); it; ++it)
        stale.push_back(it->GetTag());
    for (unsigned int t : stale) temp.DeleteData(t);
    for (const Pair& p : pairs) {
        switch (p.kind) {
            case 0: temp.SetData(p.name.c_str(), p.b); break;
            case 1: temp.SetData(p.name.c_str(), p.i); break;
            case 2: temp.SetData(p.name.c_str(), p.d); break;
            case 3: temp.SetData(p.name.c_str(), p.s); break;
            default: break;
        }
    }
}

// Build a TYPE-PRESERVING temp carrying the caller-supplied scalar data under
// this module's registry, then serialize it via `serialize`. The temp type
// mirrors write_dispatch's own threshold (src/serialize.cpp): a genuine
// multi-conformer source (OEMCMolBase with NumConfs()>1) uses an OEMol temp
// copied via the OEMCMolBase overload so all conformers survive; every other
// molecule uses an OEGraphMol temp so the molecule TITLE survives an OEB
// single-conformer write (an OEMol temp drops the title — the same
// conformer/title quirk write_dispatch guards, and SetTitle on the OEMol temp
// does NOT fix it). No transient duplication on any live molecule: only the
// throwaway temp is mutated, and its stale data is cleared before the correct
// pairs are set.
template <class Serialize>
static auto _with_bridged_temp(OEChem::OEMolBase& mol, PyObject* names,
                               PyObject* values, Serialize&& serialize) {
    auto* mc = dynamic_cast<OEChem::OEMCMolBase*>(&mol);
    if (mc && mc->NumConfs() > 1) {
        OEChem::OEMol temp;
        OEChem::OECopyMol(temp, *mc);             // structure + all conformers
        _set_bridged_data(temp, names, values);
        return serialize(temp);
    }
    OEChem::OEGraphMol temp;
    OEChem::OECopyMol(temp, mol);                 // structure + title (single conf)
    _set_bridged_data(temp, names, values);
    return serialize(temp);
}
}  // namespace oeio
%}

%inline %{
namespace oeio {
// Alias so the typemaps above bind only to the bytes entry points.
using Bytes = std::string;

// --- string (text) entry points: default std::string<->str typemaps apply ---
std::string _mol_to_string(const OEChem::OEMolBase& mol, unsigned fmt,
                           unsigned flavor) {
    return mol_to_string(mol, fmt, flavor);
}
bool _mol_from_string(OEChem::OEMolBase& mol, const std::string& data,
                      unsigned fmt, unsigned flavor) {
    return mol_from_string(mol, data, fmt, flavor);
}

// --- bytes entry points: Bytes typemaps apply (binary-safe) ---
Bytes _mol_to_bytes(const OEChem::OEMolBase& mol, unsigned fmt,
                    unsigned flavor, bool gzip) {
    return mol_to_bytes(mol, fmt, flavor, gzip);
}
bool _mol_from_bytes(OEChem::OEMolBase& mol, const Bytes& data, unsigned fmt,
                     unsigned flavor, bool gzip) {
    return mol_from_bytes(mol, data, fmt, flavor, gzip);
}

unsigned _resolve_format(const std::string& fmt) { return resolve_format(fmt); }

// Enumerate molecule-level scalar generic data as (name, value, tag) triples,
// names resolved through THIS module's OpenEye tag registry. Returns a new
// reference to a Python list. Mirrors the historical _ReaderHandle body but is
// callable without a reader handle; generic_data_pairs now projects its
// (name, value) 2-tuples from this function. Only scalar data
// (int/double/float/bool/string) is surfaced; other data is skipped.
PyObject* _scalar_generic_data(OEChem::OEMolBase& mol) {
    PyObject* lst = PyList_New(0);
    if (!lst) return NULL;
    for (OESystem::OEIter<OESystem::OEBaseData> it = mol.GetDataIter(); it; ++it) {
        const unsigned int tag = it->GetTag();
        const char* name = OESystem::OEGetTag(tag);
        if (!name || !name[0]) continue;
        const void* dtype = it->GetDataType();
        PyObject* val = NULL;
        if (dtype == OESystem::OEGetDataType<int>()) {
            val = PyLong_FromLong(mol.GetIntData(tag));
        } else if (dtype == OESystem::OEGetDataType<double>()) {
            val = PyFloat_FromDouble(mol.GetDoubleData(tag));
        } else if (dtype == OESystem::OEGetDataType<float>()) {
            val = PyFloat_FromDouble(static_cast<double>(mol.GetFloatData(tag)));
        } else if (dtype == OESystem::OEGetDataType<bool>()) {
            val = PyBool_FromLong(mol.GetBoolData(tag) ? 1 : 0);
        } else if (dtype == OESystem::OEGetDataType<std::string>()) {
            val = PyUnicode_FromString(mol.GetStringData(tag).c_str());
        } else {
            continue;  // non-scalar generic data: leave as-is
        }
        if (!val) { Py_DECREF(lst); return NULL; }
        PyObject* pyname = PyUnicode_FromString(name);
        PyObject* pytag = PyLong_FromUnsignedLong(tag);
        if (!pyname || !pytag) { Py_XDECREF(pyname); Py_XDECREF(pytag);
                                 Py_DECREF(val); Py_DECREF(lst); return NULL; }
        PyObject* tup = PyTuple_Pack(3, pyname, val, pytag);
        Py_DECREF(pyname); Py_DECREF(val); Py_DECREF(pytag);
        if (!tup) { Py_DECREF(lst); return NULL; }
        if (PyList_Append(lst, tup) != 0) { Py_DECREF(tup); Py_DECREF(lst); return NULL; }
        Py_DECREF(tup);
    }
    return lst;
}

// Build a type-preserving temp carrying the caller-supplied scalar data under
// this module's registry, then serialize to bytes. Used on static builds.
Bytes _mol_to_bytes_bridged(OEChem::OEMolBase& mol, PyObject* names,
                            PyObject* values, unsigned fmt, unsigned flavor,
                            bool gzip) {
    return _with_bridged_temp(mol, names, values,
        [&](OEChem::OEMolBase& temp) {
            return mol_to_bytes(temp, fmt, flavor, gzip);
        });
}

// Text counterpart of _mol_to_bytes_bridged (no gzip).
std::string _mol_to_string_bridged(OEChem::OEMolBase& mol, PyObject* names,
                                   PyObject* values, unsigned fmt,
                                   unsigned flavor) {
    return _with_bridged_temp(mol, names, values,
        [&](OEChem::OEMolBase& temp) {
            return mol_to_string(temp, fmt, flavor);
        });
}
}  // namespace oeio
%}

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
    // Read mol + up to grids.size() grids (filled in place via the grid-list
    // input typemap). Returns N (>=0), the record's grid count, or -1 at EOF.
    int next_grids(OEChem::OEMolBase& mol,
                   const std::vector<OESystem::OEScalarGrid*>& grids);
    // Read mol + all N grids; returns a Python tuple of native OEScalarGrid
    // objects, or None at end-of-stream. Used by with_grids().
    PyObject* next_grid_tuple(OEChem::OEMolBase& mol);
    // Enumerate the molecule's C++-set scalar generic data as (name, value)
    // tuples, with names resolved through oeio's OpenEye tag registry. Lets a
    // static build hand values back to Python for re-attach under Python's
    // registry. Returns a new reference to a list.
    PyObject* generic_data_pairs(OEChem::OEMolBase& mol);
private:
    _ReaderHandle();
};

class _WriterHandle {
public:
    // Molecule-only write (renamed from append so the Python `append` added via
    // %extend can dispatch without hiding this path). CUBE raises here.
    bool append_mol(const OEChem::OEMolBase& mol);
    // Write mol + grids (converted from a Python sequence via the input typemap).
    bool append_grids(const OEChem::OEMolBase& mol,
                      const std::vector<const OESystem::OEScalarGrid*>& grids);
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

    /// Read mol + fill up to grids.size() caller-owned grids in place.
    /// Returns N (the record's grid count, >= 0) or -1 at end-of-stream.
    int next_grids(OEChem::OEMolBase& mol,
                   const std::vector<OESystem::OEScalarGrid*>& grids) {
        if (!source_) return -1;
        int n = 0;
        if (!source_->next(mol, grids, &n)) return -1;
        return n;
    }

    /// Read mol + all N grids and wrap them as a Python tuple of native
    /// openeye.oegrid.OEScalarGrid objects. Returns a new reference to the
    /// tuple, or a new reference to Py_None at end-of-stream. Ownership of each
    /// heap grid transfers to the Python object via _oeio_wrap_as_oe_grid.
    PyObject* next_grid_tuple(OEChem::OEMolBase& mol) {
        if (!source_) Py_RETURN_NONE;
        std::vector<OESystem::OEScalarGrid> owned;
        if (!source_->next(mol, owned)) Py_RETURN_NONE;
        PyObject* tup = PyTuple_New(static_cast<Py_ssize_t>(owned.size()));
        if (!tup) return NULL;
        for (std::size_t i = 0; i < owned.size(); ++i) {
            // Copy each grid onto the heap; _oeio_wrap_as_oe_grid adopts it. The
            // copy allocates, so guard against a throw leaving the partially
            // built tuple (and its already-inserted wrappers) leaked.
            OESystem::OEScalarGrid* heap = nullptr;
            try {
                heap = new OESystem::OEScalarGrid(owned[i]);
            } catch (...) {
                Py_DECREF(tup);
                throw;
            }
            PyObject* py = _oeio_wrap_as_oe_grid(heap);
            if (!py) { Py_DECREF(tup); return NULL; }
            PyTuple_SET_ITEM(tup, static_cast<Py_ssize_t>(i), py);  // steals ref
        }
        return tup;
    }

    /// Enumerate the molecule's scalar generic data as (name, value) tuples.
    ///
    /// Names are resolved through THIS module's OpenEye tag registry — the one
    /// that set the data. When oeio is linked statically, that registry differs
    /// from the Python oechem module's, so a string tag set from C++ is not
    /// resolvable by name from Python; the value is physically on the molecule
    /// but addressed by the other registry's integer. Handing (name, value) back
    /// lets Python re-attach it under its own registry. Only scalar data
    /// (int/double/float/bool/string) is surfaced; other data is left untouched.
    /// Returns a new reference to a Python list of (str, value) tuples.
    PyObject* generic_data_pairs(OEChem::OEMolBase& mol) {
        // Project the (name, value, tag) triples from _scalar_generic_data down
        // to the historical (name, value) 2-tuple contract, so
        // Reader._reattach_cpp_data stays unchanged.
        PyObject* triples = oeio::_scalar_generic_data(mol);
        if (!triples) return NULL;
        Py_ssize_t n = PyList_Size(triples);
        PyObject* out = PyList_New(0);
        if (!out) { Py_DECREF(triples); return NULL; }
        for (Py_ssize_t i = 0; i < n; ++i) {
            PyObject* t = PyList_GetItem(triples, i);           // borrowed
            PyObject* pair = PyTuple_Pack(2, PyTuple_GetItem(t, 0),
                                             PyTuple_GetItem(t, 1));
            if (!pair) { Py_DECREF(triples); Py_DECREF(out); return NULL; }
            if (PyList_Append(out, pair) != 0) { Py_DECREF(pair);
                Py_DECREF(triples); Py_DECREF(out); return NULL; }
            Py_DECREF(pair);
        }
        Py_DECREF(triples);
        return out;
    }

private:
    std::unique_ptr<MolSource> source_;
};

class _WriterHandle {
public:
    _WriterHandle(std::unique_ptr<MolSink> sink)
        : sink_(std::move(sink)) {}

    bool append_mol(const OEChem::OEMolBase& mol) {
        return sink_ ? sink_->write(mol) : false;
    }

    bool append_grids(const OEChem::OEMolBase& mol,
                      const std::vector<const OESystem::OEScalarGrid*>& grids) {
        return sink_ ? sink_->write(mol, grids) : false;
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

    def append(self, mol, *grids):
        """Write a molecule, optionally with one or more scalar grids.

        Without grids, this writes the molecule via the molecule-only path
        (formats with no grid support write as before; the CUBE format, which
        requires at least one grid, raises ``FormatError``). With grids, it
        writes the molecule plus grids in one record.

        :param mol: An ``oechem`` molecule to write.
        :param grids: Zero or more ``oegrid.OEScalarGrid`` to write alongside.
        :returns: ``True`` on success.
        :raises FormatError: If the format cannot represent the request (e.g.
            a grid-less CUBE write).
        """
        if not grids:
            return self.append_mol(mol)
        return self.append_grids(mol, list(grids))
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


def _needs_data_reattach():
    """True when C++-set generic data must be re-attached under Python's registry.

    On a shared OpenEye link, oeio and Python's ``oechem`` share one
    ``liboesystem`` tag registry, so data set from C++ is already resolvable by
    name from Python and re-attaching is unnecessary — skipping it keeps the read
    hot path free of per-molecule enumeration (measured ~16% faster on a
    SD-tagged SDF read). On a static link the registries differ and re-attach is
    required. Unknown builds default to re-attaching (correctness over speed).

    This is evaluated once at import, not per molecule.
    """
    try:
        from . import _build_info
    except Exception:
        return True
    return getattr(_build_info, "OPENEYE_LIBRARY_TYPE", "STATIC") != "SHARED"


_NEEDS_DATA_REATTACH = _needs_data_reattach()


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

    def _reattach_cpp_data(self, mol):
        """Make C++-set scalar generic data visible from Python by string tag.

        oeio's C++ readers set typed scalars (e.g. an FCHK ``Total Energy``)
        through their own OpenEye tag registry. On a shared OpenEye link that is
        the same registry Python's ``oechem`` uses, so the data already resolves
        by name and this is a no-op (``HasData`` is already True for each). On a
        static link the two registries differ, so the value sits on the molecule
        under the other registry's integer tag and ``GetData("Total Energy")``
        misses; here each value is surfaced by name (resolved through oeio's
        registry) and re-attached under Python's registry when not already
        visible. Idempotent, so it is safe to call after every read.

        Gated on a shared build (see :func:`_needs_data_reattach`) so the read
        hot path pays nothing where re-attach would be a no-op.
        """
        if not _NEEDS_DATA_REATTACH:
            return
        for name, value in self._handle.generic_data_pairs(mol):
            if not mol.HasData(name):
                mol.SetData(name, value)

    def __iter__(self):
        from openeye import oechem

        if self._closed:
            raise ValueError("I/O operation on closed reader")
        mol = oechem.OEMol()
        while True:
            if self._closed:
                raise ValueError("I/O operation on closed reader")
            if not self._handle.next(mol):
                break
            self._reattach_cpp_data(mol)
            yield mol
            mol = oechem.OEMol()

    def __next__(self):
        """Return the next molecule, raising ``StopIteration`` at end-of-stream.

        Lets a reader be advanced one molecule at a time with the built-in
        ``next()``. Each call allocates and returns a fresh ``OEMol``.

        :returns: The next ``OEMol``.
        :raises StopIteration: At end-of-stream.
        :raises ValueError: If the reader is closed.
        """
        from openeye import oechem

        if self._closed:
            raise ValueError("I/O operation on closed reader")
        mol = oechem.OEMol()
        if not self._handle.next(mol):
            raise StopIteration
        self._reattach_cpp_data(mol)
        return mol

    def read_into(self, mol, *grids):
        """Read the next molecule (and optionally grids) into caller-owned objects.

        Without grids, this reuses ``mol`` as the buffer and returns ``True``
        while records remain, ``False`` at end-of-stream. With one or more
        grids, it additionally fills ``min(len(grids), N)`` of them and returns
        ``N`` (the record's grid count) on success, or ``None`` at
        end-of-stream. ``N == 0`` is a valid non-EOF read (a molecule with no
        grids), distinct from ``None``.

        :param mol: An ``oechem`` molecule to populate.
        :param grids: Zero or more ``oegrid.OEScalarGrid`` to fill in place.
        :returns: ``bool`` when no grids are passed; ``int`` or ``None`` when
            grids are passed.
        :raises ValueError: If the reader is closed.
        """
        if self._closed:
            raise ValueError("I/O operation on closed reader")
        if not grids:
            ok = self._handle.next(mol)
            if ok:
                self._reattach_cpp_data(mol)
            return ok
        n = self._handle.next_grids(mol, list(grids))
        if n >= 0:
            self._reattach_cpp_data(mol)
        return None if n < 0 else n

    def with_grids(self):
        """Iterate ``(molecule, grids)`` records, ``grids`` being all N grids.

        Each record yields a fresh ``OEMol`` and a tuple of native
        ``oegrid.OEScalarGrid`` objects (empty for a molecule with no grids).
        Only valid directly on a ``read(...)`` range, not after
        ``filter``/``transform`` (which yield molecules only).

        :returns: Generator of ``(OEMol, tuple-of-OEScalarGrid)``.
        :raises ValueError: If the reader is closed.
        """
        from openeye import oechem

        if self._closed:
            raise ValueError("I/O operation on closed reader")
        while True:
            if self._closed:
                raise ValueError("I/O operation on closed reader")
            mol = oechem.OEMol()
            grids = self._handle.next_grid_tuple(mol)
            if grids is None:
                break
            self._reattach_cpp_data(mol)
            yield mol, grids

    def as_type(self, cls):
        """Iterate molecules as instances of ``cls``.

        :param cls: An ``oechem`` molecule class (``OEMol``, ``OEGraphMol``,
            or ``OEQMol``).
        :returns: A generator yielding fresh ``cls`` instances.
        :raises ValueError: If the reader is closed.
        :raises TypeError: If ``cls`` is not an ``OEMolBase`` subclass.
        """
        from openeye import oechem

        if self._closed:
            raise ValueError("I/O operation on closed reader")
        if not (isinstance(cls, type) and issubclass(cls, oechem.OEMolBase)):
            raise TypeError(
                "as_type() requires an oechem OEMolBase subclass, "
                "got {!r}".format(cls))
        mol = cls()
        while True:
            if self._closed:
                raise ValueError("I/O operation on closed reader")
            if not self._handle.next(mol):
                break
            self._reattach_cpp_data(mol)
            yield mol
            mol = cls()

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
    for mol in iterable:
        if predicate(mol):
            yield type(mol)(mol)


def transform(iterable, func):
    """Transform molecules from an iterable in-place.

    :param iterable: An iterable of molecules (e.g., from read()).
    :param func: Function that takes a molecule (OEMolBase) and modifies it in-place.
    :returns: Generator yielding transformed molecules.

    Example::

        prepared = oeio.transform(oeio.read("in.sdf"),
            lambda mol: oechem.OEAddExplicitHydrogens(mol))
    """
    for mol in iterable:
        func(mol)
        yield type(mol)(mol)


def formats():
    """List all registered molecular file formats.

    :returns: List of FormatInfo objects.
    """
    return list(FormatRegistry.instance().formats())


__version__ = "0.4.0"

_install_exception_types(Error, FormatError, FileError)
%}
