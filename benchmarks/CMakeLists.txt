include(FetchContent)

FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.1
    GIT_SHALLOW TRUE
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark)

# Data generation utility
add_executable(oeio_generate_bench_data cpp/generate_bench_data.cpp)
target_link_libraries(oeio_generate_bench_data PRIVATE oeio)

# C++ benchmarks
add_executable(oeio_bench_read cpp/bench_read.cpp)
target_link_libraries(oeio_bench_read PRIVATE oeio benchmark::benchmark_main)

add_executable(oeio_bench_pipeline cpp/bench_pipeline.cpp)
target_link_libraries(oeio_bench_pipeline PRIVATE oeio benchmark::benchmark_main)

add_executable(oeio_bench_registry cpp/bench_registry.cpp)
target_link_libraries(oeio_bench_registry PRIVATE oeio benchmark::benchmark_main)
