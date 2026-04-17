/// \file bench_read.cpp
/// \brief Benchmark oeio::read() vs raw OEChem molecule reading.

#include <benchmark/benchmark.h>

#include "oeio/oeio.h"

#include <oechem.h>

#include <cstdlib>
#include <filesystem>
#include <string>

// Force-link the OEChem handler for static builds
namespace oeio { void oeio_force_link_oechem_handler(); }
static struct _BenchForceLink {
    _BenchForceLink() { oeio::oeio_force_link_oechem_handler(); }
} _bench_force_link;

static std::string get_bench_file() {
    const char* path = std::getenv("OEIO_BENCH_FILE");
    if (!path) {
        static bool warned = false;
        if (!warned) {
            fprintf(stderr, "Set OEIO_BENCH_FILE to a molecule file for benchmarking.\n"
                            "Example: OEIO_BENCH_FILE=/tmp/bench_10k.sdf ./oeio_bench_read\n");
            warned = true;
        }
        return "";
    }
    return path;
}

static void BM_OeioRead(benchmark::State& state) {
    std::string path = get_bench_file();
    if (path.empty()) {
        state.SkipWithError("OEIO_BENCH_FILE not set");
        return;
    }
    for (auto _ : state) {
        int count = 0;
        for (auto& mol : oeio::read(path)) {
            benchmark::DoNotOptimize(mol);
            ++count;
        }
        state.counters["molecules"] = count;
    }
}
BENCHMARK(BM_OeioRead)->Unit(benchmark::kMillisecond);

static void BM_RawOEChemRead(benchmark::State& state) {
    std::string path = get_bench_file();
    if (path.empty()) {
        state.SkipWithError("OEIO_BENCH_FILE not set");
        return;
    }
    for (auto _ : state) {
        OEChem::oemolistream ifs;
        ifs.open(path);
        OEChem::OEGraphMol mol;
        int count = 0;
        while (OEChem::OEReadMolecule(ifs, mol)) {
            benchmark::DoNotOptimize(mol);
            ++count;
        }
        state.counters["molecules"] = count;
    }
}
BENCHMARK(BM_RawOEChemRead)->Unit(benchmark::kMillisecond);
