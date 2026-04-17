"""Benchmark oeio.read() vs raw oechem.oemolistream.

Usage:
    python bench_read.py <molecule_file> [--iterations N]
"""

import argparse
import time

from openeye import oechem


def bench_raw_oechem(path, iterations=5):
    """Benchmark raw OEChem molecule reading."""
    times = []
    for _ in range(iterations):
        ifs = oechem.oemolistream()
        ifs.open(path)
        mol = oechem.OEGraphMol()
        start = time.perf_counter()
        count = 0
        while oechem.OEReadMolecule(ifs, mol):
            count += 1
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return times, count


def bench_oeio_read(path, iterations=5):
    """Benchmark oeio.read() molecule reading."""
    import oeio

    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        count = 0
        for mol in oeio.read(path):
            count += 1
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return times, count


def main():
    parser = argparse.ArgumentParser(description="Benchmark oeio.read() vs raw OEChem")
    parser.add_argument("path", help="Path to molecule file")
    parser.add_argument("--iterations", type=int, default=5, help="Number of iterations")
    args = parser.parse_args()

    print(f"Benchmarking: {args.path}")
    print(f"Iterations: {args.iterations}")
    print()

    raw_times, count = bench_raw_oechem(args.path, args.iterations)
    print(f"Raw OEChem:  {count} molecules")
    print(f"  mean: {sum(raw_times)/len(raw_times)*1000:.1f} ms")
    print(f"  min:  {min(raw_times)*1000:.1f} ms")
    print()

    oeio_times, count = bench_oeio_read(args.path, args.iterations)
    print(f"oeio.read(): {count} molecules")
    print(f"  mean: {sum(oeio_times)/len(oeio_times)*1000:.1f} ms")
    print(f"  min:  {min(oeio_times)*1000:.1f} ms")
    print()

    ratio = min(oeio_times) / min(raw_times) if min(raw_times) > 0 else float('inf')
    print(f"Ratio (oeio / raw): {ratio:.2f}x")


if __name__ == "__main__":
    main()
