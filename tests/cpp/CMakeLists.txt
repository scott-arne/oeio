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
)

target_link_libraries(oeio_tests
    PRIVATE
        oeio
        GTest::gtest_main
)

target_compile_definitions(oeio_tests PRIVATE
    OEIO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")

include(GoogleTest)
gtest_discover_tests(oeio_tests)
