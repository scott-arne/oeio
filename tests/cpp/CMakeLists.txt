include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(oeio_tests
    test_registry.cpp
    test_extension_matching.cpp
    test_oechem_handler.cpp
    test_cube_handler.cpp
    test_fchk_handler.cpp
    test_mol_range.cpp
    test_pipeline.cpp
    test_exceptions.cpp
    test_serialize.cpp
    test_read_status.cpp
)

target_link_libraries(oeio_tests
    PRIVATE
        oeio
        GTest::gtest_main
)

# Workaround: tests use OEAddExplicitHydrogens which is in OEQuacPac
if(OPENEYE_USE_SHARED AND EXISTS "/Users/johnss51/Applications/micromamba/envs/main/lib/python3.13/site-packages/openeye/libs/python3-osx-universal-clang++/liboequacpac-2.3.0.0.dylib")
    target_link_libraries(oeio_tests PRIVATE "/Users/johnss51/Applications/micromamba/envs/main/lib/python3.13/site-packages/openeye/libs/python3-osx-universal-clang++/liboequacpac-2.3.0.0.dylib")
endif()

target_compile_definitions(oeio_tests PRIVATE
    OEIO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")

include(GoogleTest)
gtest_discover_tests(oeio_tests)
