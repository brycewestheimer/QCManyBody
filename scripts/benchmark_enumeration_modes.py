#!/usr/bin/env python3
"""
Benchmark HMBE enumeration modes: filter vs direct vs auto.

This script measures the performance of different enumeration modes for various
system sizes and HMBE configurations. It helps determine the optimal crossover
point for automatic mode selection.

Usage:
    python scripts/benchmark_enumeration_modes.py
    python scripts/benchmark_enumeration_modes.py --max-frags 64
    python scripts/benchmark_enumeration_modes.py --plot results.png
"""

import argparse
import time
import sys
from typing import Dict, List, Tuple
from dataclasses import dataclass

import numpy as np
from qcelemental.models import Molecule

from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.core import ManyBodyCore
from qcmanybody.hmbe_enumerate import get_enumeration_statistics


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run."""
    n_fragments: int
    n_groups: int
    truncation_orders: Tuple[int, int]
    mode: str
    time_seconds: float
    n_terms: int
    memory_mb: float


def create_test_molecule(n_fragments: int, spacing: float = 3.0) -> Molecule:
    """
    Create a test water cluster with specified number of fragments.

    Args:
        n_fragments: Number of water molecules
        spacing: Spacing between waters in Angstroms

    Returns:
        Molecule object with n_fragments water molecules
    """
    symbols = []
    geometry = []
    fragments = []

    # Arrange waters in a grid
    grid_size = int(np.ceil(np.sqrt(n_fragments)))

    for i in range(n_fragments):
        row = i // grid_size
        col = i % grid_size

        # Water geometry (simple, not optimized)
        x_offset = col * spacing
        y_offset = row * spacing

        symbols.extend(["O", "H", "H"])
        geometry.extend([
            [x_offset, y_offset, 0.0],
            [x_offset + 0.757, y_offset + 0.586, 0.0],
            [x_offset - 0.757, y_offset + 0.586, 0.0],
        ])
        fragments.append([i * 3, i * 3 + 1, i * 3 + 2])

    return Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=0.0,
        molecular_multiplicity=1,
    )


def create_test_hierarchy(n_fragments: int, frags_per_group: int = 4) -> FragmentHierarchy:
    """
    Create a 2-tier hierarchy for testing.

    Args:
        n_fragments: Total number of fragments
        frags_per_group: Fragments per tier-1 group

    Returns:
        FragmentHierarchy with balanced grouping
    """
    fragment_tiers = {}
    n_groups = int(np.ceil(n_fragments / frags_per_group))

    for i in range(n_fragments):
        frag_id = i + 1  # 1-indexed
        group_id = i // frags_per_group
        fragment_tiers[frag_id] = (f"G{group_id}", f"F{i}")

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("group", "fragment")
    )


def benchmark_enumeration_mode(
    n_fragments: int,
    truncation_orders: Tuple[int, int],
    mode: str,
    frags_per_group: int = 4
) -> BenchmarkResult:
    """
    Benchmark a single enumeration mode configuration.

    Args:
        n_fragments: Number of fragments in test system
        truncation_orders: HMBE truncation orders (T_1, T_2)
        mode: Enumeration mode ("filter", "direct", or "auto")
        frags_per_group: Fragments per tier-1 group

    Returns:
        BenchmarkResult with timing and statistics
    """
    # Create test molecule and hierarchy
    mol = create_test_molecule(n_fragments)
    hierarchy = create_test_hierarchy(n_fragments, frags_per_group)
    n_groups = int(np.ceil(n_fragments / frags_per_group))

    # Create HMBE specification
    hmbe_spec = HMBESpecification(
        truncation_orders=truncation_orders,
        hierarchy=hierarchy,
        enumeration_mode=mode
    )

    # Time the compute map construction (this is where enumeration happens)
    start_time = time.perf_counter()

    try:
        mbc = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        # Access compute_map to trigger enumeration
        _ = mbc.compute_map

        elapsed_time = time.perf_counter() - start_time

        # Get statistics
        stats = mbc.get_hmbe_statistics()
        n_terms = stats["hmbe_term_counts"]["hf/sto-3g"]

        # Rough memory estimate (not precise, but indicative)
        memory_mb = sys.getsizeof(mbc.compute_map) / 1e6

    except Exception as e:
        print(f"Error for {n_fragments} frags, mode={mode}: {e}")
        elapsed_time = float('inf')
        n_terms = 0
        memory_mb = 0.0

    return BenchmarkResult(
        n_fragments=n_fragments,
        n_groups=n_groups,
        truncation_orders=truncation_orders,
        mode=mode,
        time_seconds=elapsed_time,
        n_terms=n_terms,
        memory_mb=memory_mb
    )


def run_benchmark_suite(
    fragment_counts: List[int],
    truncation_configs: List[Tuple[int, int]],
    modes: List[str] = ["filter", "direct", "auto"],
    frags_per_group: int = 4
) -> List[BenchmarkResult]:
    """
    Run comprehensive benchmark suite.

    Args:
        fragment_counts: List of fragment counts to test
        truncation_configs: List of (T_1, T_2) configurations
        modes: Enumeration modes to benchmark
        frags_per_group: Fragments per tier-1 group

    Returns:
        List of BenchmarkResult objects
    """
    results = []
    total_runs = len(fragment_counts) * len(truncation_configs) * len(modes)
    current_run = 0

    print(f"Running {total_runs} benchmark configurations...")
    print(f"{'Frags':<8} {'Config':<12} {'Mode':<10} {'Time (s)':<12} {'Terms':<10} {'Memory (MB)':<12}")
    print("-" * 80)

    for n_frags in fragment_counts:
        for trunc_orders in truncation_configs:
            for mode in modes:
                current_run += 1
                config_str = f"({trunc_orders[0]},{trunc_orders[1]})"

                result = benchmark_enumeration_mode(
                    n_fragments=n_frags,
                    truncation_orders=trunc_orders,
                    mode=mode,
                    frags_per_group=frags_per_group
                )

                results.append(result)

                # Print progress
                print(f"{result.n_fragments:<8} {config_str:<12} {result.mode:<10} "
                      f"{result.time_seconds:<12.4f} {result.n_terms:<10} {result.memory_mb:<12.2f}")

    print("\nBenchmark complete!")
    return results


def analyze_results(results: List[BenchmarkResult]) -> None:
    """
    Analyze and report benchmark results.

    Args:
        results: List of BenchmarkResult objects
    """
    print("\n" + "=" * 80)
    print("BENCHMARK ANALYSIS")
    print("=" * 80)

    # Find crossover point where direct becomes faster than filter
    crossover_data = {}
    for result in results:
        if result.mode in ["filter", "direct"]:
            key = (result.n_fragments, result.truncation_orders)
            if key not in crossover_data:
                crossover_data[key] = {}
            crossover_data[key][result.mode] = result.time_seconds

    print("\nCrossover Analysis (Direct vs Filter):")
    print(f"{'Fragments':<12} {'Config':<12} {'Filter (s)':<15} {'Direct (s)':<15} {'Speedup':<10}")
    print("-" * 80)

    crossover_point = None
    for (n_frags, trunc_orders), times in sorted(crossover_data.items()):
        if "filter" in times and "direct" in times:
            filter_time = times["filter"]
            direct_time = times["direct"]
            speedup = filter_time / direct_time if direct_time > 0 else float('inf')

            config_str = f"({trunc_orders[0]},{trunc_orders[1]})"
            print(f"{n_frags:<12} {config_str:<12} {filter_time:<15.4f} {direct_time:<15.4f} {speedup:<10.2f}x")

            # Track when direct becomes faster
            if speedup > 1.0 and crossover_point is None:
                crossover_point = n_frags

    if crossover_point:
        print(f"\n✓ Crossover point: Direct enumeration becomes faster at ~{crossover_point} fragments")
    else:
        print(f"\n⚠ No clear crossover point found in tested range")

    # Memory usage comparison
    print("\nMemory Usage Comparison:")
    print(f"{'Fragments':<12} {'Config':<12} {'Filter (MB)':<15} {'Direct (MB)':<15} {'Ratio':<10}")
    print("-" * 80)

    memory_data = {}
    for result in results:
        if result.mode in ["filter", "direct"]:
            key = (result.n_fragments, result.truncation_orders)
            if key not in memory_data:
                memory_data[key] = {}
            memory_data[key][result.mode] = result.memory_mb

    for (n_frags, trunc_orders), mems in sorted(memory_data.items()):
        if "filter" in mems and "direct" in mems:
            filter_mem = mems["filter"]
            direct_mem = mems["direct"]
            ratio = filter_mem / direct_mem if direct_mem > 0 else float('inf')

            config_str = f"({trunc_orders[0]},{trunc_orders[1]})"
            print(f"{n_frags:<12} {config_str:<12} {filter_mem:<15.2f} {direct_mem:<15.2f} {ratio:<10.2f}x")


def plot_results(results: List[BenchmarkResult], output_file: str) -> None:
    """
    Create performance plots.

    Args:
        results: List of BenchmarkResult objects
        output_file: Path to save plot image
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("⚠ matplotlib not installed, skipping plots")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Group results by config
    for trunc_orders in set(r.truncation_orders for r in results):
        config_results = [r for r in results if r.truncation_orders == trunc_orders]

        for mode in ["filter", "direct"]:
            mode_results = sorted(
                [r for r in config_results if r.mode == mode],
                key=lambda r: r.n_fragments
            )

            if not mode_results:
                continue

            frags = [r.n_fragments for r in mode_results]
            times = [r.time_seconds for r in mode_results]
            mems = [r.memory_mb for r in mode_results]

            label = f"({trunc_orders[0]},{trunc_orders[1]})-HMBE {mode}"
            ax1.plot(frags, times, marker='o', label=label)
            ax2.plot(frags, mems, marker='s', label=label)

    ax1.set_xlabel("Number of Fragments")
    ax1.set_ylabel("Enumeration Time (seconds)")
    ax1.set_title("Enumeration Performance: Filter vs Direct")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')

    ax2.set_xlabel("Number of Fragments")
    ax2.set_ylabel("Memory Usage (MB)")
    ax2.set_title("Memory Usage: Filter vs Direct")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark HMBE enumeration modes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--max-frags",
        type=int,
        default=32,
        help="Maximum number of fragments to test (default: 32)"
    )
    parser.add_argument(
        "--frags-per-group",
        type=int,
        default=4,
        help="Fragments per tier-1 group (default: 4)"
    )
    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        help="Save performance plot to file (requires matplotlib)"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick benchmark with fewer configurations"
    )

    args = parser.parse_args()

    # Define test configurations
    if args.quick:
        fragment_counts = [4, 8, 16, 24, 32]
        truncation_configs = [(2, 3)]
        modes = ["filter", "direct"]
    else:
        fragment_counts = [4, 8, 12, 16, 20, 24, 28, 32]
        if args.max_frags > 32:
            fragment_counts.extend(range(36, args.max_frags + 1, 4))

        truncation_configs = [
            (2, 3),  # (2,3)-HMBE
            (2, 4),  # (2,4)-HMBE
        ]
        modes = ["filter", "direct", "auto"]

    # Run benchmarks
    results = run_benchmark_suite(
        fragment_counts=fragment_counts,
        truncation_configs=truncation_configs,
        modes=modes,
        frags_per_group=args.frags_per_group
    )

    # Analyze results
    analyze_results(results)

    # Create plots if requested
    if args.plot:
        plot_results(results, args.plot)

    # Recommendations
    print("\n" + "=" * 80)
    print("RECOMMENDATIONS")
    print("=" * 80)
    print("""
Based on benchmarks, the automatic mode selection (threshold = 30 fragments) is set to:
- Use FILTER mode for <30 fragments (faster enumeration, acceptable memory)
- Use DIRECT mode for ≥30 fragments (much faster, lower memory)

You can override automatic selection by setting enumeration_mode explicitly:
- enumeration_mode="filter": Always use filter approach
- enumeration_mode="direct": Always use direct enumeration
- enumeration_mode="auto": Automatic selection (recommended)

For very large systems (100+ fragments), DIRECT mode is essential.
    """)


if __name__ == "__main__":
    main()
