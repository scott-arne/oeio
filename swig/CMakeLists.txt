include(OpenEyeSWIG)

openeye_add_swig_module(
    NAME oeio
    SWIG_FILE oeio.i
    LINK_LIBS oeio
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oeio
    STABLE_ABI ${OEIO_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen"
    EXPECTED_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
        OEMATH_LIBRARY OEZSTD_LIBRARY
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oeio/__init__.py
)
