"""Benchmark oeio.batch_read() vs oeio.read() vs raw OEChem.

Usage:
    python bench_batch_read.py <molecule_file> [--iterations N]
"""

import argparse
import time

from openeye import oechem


def bench_raw_oechem(path, iterations=5):
    """Benchmark raw OEChem reading."""
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
    """Benchmark oeio.read()."""
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


def bench_oeio_batch_read(path, batch_size=1000, iterations=5):
    """Benchmark oeio.batch_read()."""
    import oeio

    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        count = 0
        for mol in oeio.batch_read(path, batch_size=batch_size):
            count += 1
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return times, count


def main():
    parser = argparse.ArgumentParser(description="Benchmark batch_read vs read vs raw OEChem")
    parser.add_argument("path", help="Path to molecule file")
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=1000)
    args = parser.parse_args()

    print(f"Benchmarking: {args.path}")
    print(f"Iterations: {args.iterations}, Batch size: {args.batch_size}")
    print()

    raw_times, count = bench_raw_oechem(args.path, args.iterations)
    raw_min = min(raw_times)
    print(f"Raw OEChem:       {count} molecules, min: {raw_min*1000:.1f} ms")

    read_times, _ = bench_oeio_read(args.path, args.iterations)
    read_min = min(read_times)
    print(f"oeio.read():      {count} molecules, min: {read_min*1000:.1f} ms  ({read_min/raw_min:.2f}x)")

    batch_times, _ = bench_oeio_batch_read(args.path, args.batch_size, args.iterations)
    batch_min = min(batch_times)
    print(f"oeio.batch_read():{count} molecules, min: {batch_min*1000:.1f} ms  ({batch_min/raw_min:.2f}x)")


if __name__ == "__main__":
    main()
