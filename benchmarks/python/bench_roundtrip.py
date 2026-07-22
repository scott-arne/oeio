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


def bench_inmemory_bytes(path, iterations=20000):
    """Benchmark oeio.to_bytes/from_bytes (OEB) vs oemolistream openstring.

    Uses the first molecule from ``path`` and times a single-molecule OEB
    round-trip (serialize + deserialize) for each mechanism.
    """
    import oeio

    ifs = oechem.oemolistream()
    ifs.open(path)
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()

    OEB = oechem.OEFormat_OEB

    def fast():
        b = oeio.to_bytes(mol)  # OEB bytes
        oeio.from_bytes(b)

    def openstring():
        oms = oechem.oemolostream()
        oms.SetFormat(OEB)
        oms.openstring()
        oechem.OEWriteMolecule(oms, mol)
        s = oms.GetString()
        ims = oechem.oemolistream()
        ims.SetFormat(OEB)
        ims.openstring(s)
        d = oechem.OEGraphMol()
        oechem.OEReadMolecule(ims, d)

    results = {}
    for label, fn in (("oeio.to/from_bytes (OEB)", fast), ("openstring (OEB)", openstring)):
        fn()  # warmup
        start = time.perf_counter()
        for _ in range(iterations):
            fn()
        results[label] = (time.perf_counter() - start) * 1e6 / iterations
    return results


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
    print()

    print("In-memory single-molecule bytes round-trip (first molecule):")
    mem = bench_inmemory_bytes(args.path)
    for label, us in mem.items():
        print(f"  {label:28s}: {us:7.1f} us/round-trip")
    fast_us = mem.get("oeio.to/from_bytes (OEB)")
    open_us = mem.get("openstring (OEB)")
    if fast_us and open_us and fast_us > 0:
        print(f"  speedup (openstring / oeio-bytes): {open_us / fast_us:.2f}x")


if __name__ == "__main__":
    main()
