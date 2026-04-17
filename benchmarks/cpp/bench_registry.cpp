/// \file bench_registry.cpp
/// \brief Benchmark FormatRegistry lookup latency.

#include <benchmark/benchmark.h>

#include "oeio/oeio.h"

namespace oeio { void oeio_force_link_oechem_handler(); }
static struct _BenchForceLink {
    _BenchForceLink() { oeio::oeio_force_link_oechem_handler(); }
} _bench_force_link;

static void BM_LookupSimpleExt(benchmark::State& state) {
    const auto& registry = oeio::FormatRegistry::instance();
    for (auto _ : state) {
        auto* handler = registry.lookup("molecules.sdf");
        benchmark::DoNotOptimize(handler);
    }
}
BENCHMARK(BM_LookupSimpleExt);

static void BM_LookupCompoundExt(benchmark::State& state) {
    const auto& registry = oeio::FormatRegistry::instance();
    for (auto _ : state) {
        auto* handler = registry.lookup("molecules.sdf.gz");
        benchmark::DoNotOptimize(handler);
    }
}
BENCHMARK(BM_LookupCompoundExt);

static void BM_LookupMiss(benchmark::State& state) {
    const auto& registry = oeio::FormatRegistry::instance();
    for (auto _ : state) {
        auto* handler = registry.lookup("molecules.unknown");
        benchmark::DoNotOptimize(handler);
    }
}
BENCHMARK(BM_LookupMiss);

static void BM_LookupExtDirect(benchmark::State& state) {
    const auto& registry = oeio::FormatRegistry::instance();
    for (auto _ : state) {
        auto* handler = registry.lookup_ext(".sdf");
        benchmark::DoNotOptimize(handler);
    }
}
BENCHMARK(BM_LookupExtDirect);
