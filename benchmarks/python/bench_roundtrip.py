"""Benchmark oeio read-write roundtrip vs raw OEChem.

Usage:
    python bench_roundtrip.py <molecule_file> [--iterations N]
"""

import argparse
import os
import tempfile
import time

from openeye import oechem


def bench_raw_roundtrip(path, iterations=5):
    """Benchmark raw OEChem read-write roundtrip."""
    times = []
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        out_path = f.name
    try:
        for _ in range(iterations):
            ifs = oechem.oemolistream()
            ifs.open(path)
            ofs = oechem.oemolostream()
            ofs.open(out_path)
            mol = oechem.OEGraphMol()
            start = time.perf_counter()
            count = 0
            while oechem.OEReadMolecule(ifs, mol):
                oechem.OEWriteMolecule(ofs, mol)
                count += 1
            ofs.close()
            elapsed = time.perf_counter() - start
            times.append(elapsed)
    finally:
        os.unlink(out_path)
    return times, count


def bench_oeio_roundtrip(path, iterations=5):
    """Benchmark oeio read-write roundtrip."""
    import oeio

    times = []
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        out_path = f.name
    try:
        for _ in range(iterations):
            start = time.perf_counter()
            count = 0
            with oeio.write(out_path) as writer:
                for mol in oeio.read(path):
                    writer.append(mol)
                    count += 1
            elapsed = time.perf_counter() - start
            times.append(elapsed)
    finally:
        os.unlink(out_path)
    return times, count


def main():
    parser = argparse.ArgumentParser(description="Benchmark oeio roundtrip vs raw OEChem")
    parser.add_argument("path", help="Path to molecule file")
    parser.add_argument("--iterations", type=int, default=5, help="Number of iterations")
    args = parser.parse_args()

    print(f"Benchmarking roundtrip: {args.path}")
    print(f"Iterations: {args.iterations}")
    print()

    raw_times, count = bench_raw_roundtrip(args.path, args.iterations)
    print(f"Raw OEChem:    {count} molecules")
    print(f"  mean: {sum(raw_times)/len(raw_times)*1000:.1f} ms")
    print(f"  min:  {min(raw_times)*1000:.1f} ms")
    print()

    oeio_times, count = bench_oeio_roundtrip(args.path, args.iterations)
    print(f"oeio roundtrip: {count} molecules")
    print(f"  mean: {sum(oeio_times)/len(oeio_times)*1000:.1f} ms")
    print(f"  min:  {min(oeio_times)*1000:.1f} ms")
    print()

    ratio = min(oeio_times) / min(raw_times) if min(raw_times) > 0 else float('inf')
    print(f"Ratio (oeio / raw): {ratio:.2f}x")


if __name__ == "__main__":
    main()
