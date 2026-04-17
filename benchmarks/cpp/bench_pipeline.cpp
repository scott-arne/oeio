/// \file bench_pipeline.cpp
/// \brief Benchmark oeio pipeline vs hand-written OEChem loop.

#include <benchmark/benchmark.h>

#include "oeio/oeio.h"

#include <oechem.h>

#include <cstdlib>
#include <filesystem>
#include <string>

namespace oeio { void oeio_force_link_oechem_handler(); }
static struct _BenchForceLink {
    _BenchForceLink() { oeio::oeio_force_link_oechem_handler(); }
} _bench_force_link;

static std::string get_bench_file() {
    const char* path = std::getenv("OEIO_BENCH_FILE");
    return path ? path : "";
}

static void BM_OeioPipeline(benchmark::State& state) {
    std::string in_path = get_bench_file();
    if (in_path.empty()) {
        state.SkipWithError("OEIO_BENCH_FILE not set");
        return;
    }
    auto tmp = std::filesystem::temp_directory_path() / "oeio_bench_pipeline.sdf";
    std::string out_path = tmp.string();

    for (auto _ : state) {
        oeio::transform(
            oeio::filter(oeio::read(in_path), [](const OEChem::OEMolBase& mol) {
                return mol.NumAtoms() > 3;
            }),
            [](OEChem::OEGraphMol& mol) {
                mol.SetTitle("processed");
            }
        ) | oeio::write(out_path);
    }
    std::filesystem::remove(tmp);
}
BENCHMARK(BM_OeioPipeline)->Unit(benchmark::kMillisecond);

static void BM_RawOEChemPipeline(benchmark::State& state) {
    std::string in_path = get_bench_file();
    if (in_path.empty()) {
        state.SkipWithError("OEIO_BENCH_FILE not set");
        return;
    }
    auto tmp = std::filesystem::temp_directory_path() / "oeio_bench_raw_pipeline.sdf";
    std::string out_path = tmp.string();

    for (auto _ : state) {
        OEChem::oemolistream ifs;
        ifs.open(in_path);
        OEChem::oemolostream ofs;
        ofs.open(out_path);
        OEChem::OEGraphMol mol;
        while (OEChem::OEReadMolecule(ifs, mol)) {
            if (mol.NumAtoms() > 3) {
                mol.SetTitle("processed");
                OEChem::OEWriteMolecule(ofs, mol);
            }
        }
        ofs.close();
    }
    std::filesystem::remove(tmp);
}
BENCHMARK(BM_RawOEChemPipeline)->Unit(benchmark::kMillisecond);
