"""
Tutorial 3: Parallel Execution for HMBE

This tutorial demonstrates how to use parallel execution with HMBE calculations
to achieve significant speedups on multi-core systems.

Key concepts:
- Enabling parallel execution
- Choosing number of workers
- Performance expectations
- Monitoring progress

System: 12 water molecules (3×4 hierarchy)
Parallel: 4 workers (multiprocessing)
Expected speedup: ~3-4x on 4-core system (80% parallel efficiency)

NOTE: This tutorial shows the API for parallel execution. To actually run
parallel QC calculations, you need QCEngine and a QC program (e.g., Psi4).
"""

from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models.v1 import BsseEnum, ManyBodyInput, ManyBodyKeywords
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification


def create_test_molecule(num_waters=12):
    """Create a simple water cluster for parallel testing."""

    symbols = ["O", "H", "H"] * num_waters
    geometry = []
    fragments = []

    for i in range(num_waters):
        x = (i % 4) * 3.0  # 4 waters per row
        y = (i // 4) * 3.0
        # O-H-H
        geometry.extend([
            [x, y, 0.0],
            [x + 0.757, y, 0.587],
            [x - 0.757, y, 0.587]
        ])
        start = i * 3
        fragments.append([start, start + 1, start + 2])

    return Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments
    )


def create_test_hierarchy(num_waters=12, waters_per_group=4):
    """Create hierarchy for parallel testing."""

    fragment_tiers = {}
    num_groups = num_waters // waters_per_group

    for i in range(num_waters):
        frag_id = i + 1
        group_id = i // waters_per_group
        fragment_tiers[frag_id] = (f"G{group_id}", f"W{i}")

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("group", "water")
    )


def demonstrate_parallel_api():
    """Show different ways to use parallel execution."""

    print("="*80)
    print("PARALLEL EXECUTION API")
    print("="*80)
    print()

    print("Method 1: Simple parallel flag")
    print("-"*80)
    print("""
    result = ManyBodyComputer.from_manybodyinput(
        mb_input,
        parallel=True,         # Enable parallel execution
        n_workers=4            # Use 4 CPU cores
    )
    """)
    print()

    print("Method 2: Auto-detect number of workers")
    print("-"*80)
    print("""
    result = ManyBodyComputer.from_manybodyinput(
        mb_input,
        parallel=True          # Uses all available cores
    )
    """)
    print()

    print("Method 3: Custom executor (advanced)")
    print("-"*80)
    print("""
    from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

    config = ExecutorConfig(
        n_workers=8,
        timeout_per_task=1800,  # 30 min timeout per task
        max_retries=3
    )
    executor = MultiprocessingExecutor(config)

    result = ManyBodyComputer.from_manybodyinput(
        mb_input,
        executor=executor
    )
    """)
    print()

    print("="*80)
    print()


def explain_performance_expectations():
    """Explain what speedups to expect."""

    print("="*80)
    print("PERFORMANCE EXPECTATIONS")
    print("="*80)
    print()

    print("Theoretical Maximum Speedup:")
    print("-"*80)
    print()
    print("  Workers  | Ideal Speedup | Typical Actual | Parallel Efficiency")
    print("  ---------|---------------|----------------|--------------------")
    print("     1     |      1.0x     |     1.0x       |       100%")
    print("     2     |      2.0x     |     1.8x       |        90%")
    print("     4     |      4.0x     |     3.2x       |        80%")
    print("     8     |      8.0x     |     5.6x       |        70%")
    print("    16     |     16.0x     |     9.6x       |        60%")
    print()

    print("Why Not 100% Efficiency?")
    print("-"*80)
    print("  1. Task setup/teardown overhead")
    print("  2. Inter-process communication")
    print("  3. Load imbalance (some tasks take longer)")
    print("  4. Sequential bottlenecks (aggregation, I/O)")
    print()

    print("Sweet Spot:")
    print("-"*80)
    print("  For most systems: 4-8 workers gives best efficiency/speedup balance")
    print("  Diminishing returns beyond 16 workers for typical calculations")
    print()

    print("="*80)
    print()


def demonstrate_monitoring(verbose=True):
    """Show how to monitor parallel execution progress."""

    print("="*80)
    print("MONITORING PARALLEL EXECUTION")
    print("="*80)
    print()

    if verbose:
        print("When running with logging enabled, you'll see:")
        print("-"*80)
        print("""
    [INFO] Starting parallel execution with 4 workers
    [INFO] Building 142 fragment calculations
    [INFO] Level 1: Executing 12 1-body calculations
    [INFO] Level 1: Completed 12/12 (100%)
    [INFO] Level 2: Executing 66 2-body calculations
    [INFO] Level 2: Completed 66/66 (100%)
    [INFO] Level 3: Executing 64 3-body calculations
    [INFO] Level 3: Completed 64/64 (100%)
    [INFO] All calculations complete
    [INFO] Total time: 15.2 minutes
    [INFO] Effective speedup: 3.4x
        """)
        print()

    print("Enabling detailed logging:")
    print("-"*80)
    print("""
    import logging
    logging.basicConfig(level=logging.INFO)

    # Then run your calculation
    result = ManyBodyComputer.from_manybodyinput(mb_input, parallel=True)
    """)
    print()

    print("="*80)
    print()


def demonstrate_best_practices():
    """Show best practices for parallel HMBE."""

    print("="*80)
    print("BEST PRACTICES FOR PARALLEL HMBE")
    print("="*80)
    print()

    print("1. Choose Appropriate Number of Workers")
    print("-"*80)
    print("  ✓ Match to your CPU cores (leave 1-2 for OS)")
    print("  ✓ For 8-core machine: use n_workers=6-7")
    print("  ✗ Don't use more workers than cores (causes slowdown)")
    print()

    print("2. Memory Considerations")
    print("-"*80)
    print("  • Each worker needs memory for its QC calculation")
    print("  • Rule of thumb: (Total RAM / n_workers) should be >2 GB")
    print("  • For large basis sets: may need fewer workers")
    print()

    print("3. HMBE + Parallel = Maximum Efficiency")
    print("-"*80)
    print("  • HMBE reduces number of calculations")
    print("  • Parallel speeds up remaining calculations")
    print("  • Combined speedup: HMBE_reduction × parallel_speedup")
    print("  • Example: 3x HMBE reduction × 3.5x parallel = ~10x total!")
    print()

    print("4. Task Granularity")
    print("-"*80)
    print("  • HMBE creates many small tasks (good for load balancing)")
    print("  • Parallel executor distributes tasks efficiently")
    print("  • Works well even with heterogeneous task durations")
    print()

    print("5. When NOT to Use Parallel")
    print("-"*80)
    print("  ✗ Very small systems (<20 fragments) - overhead dominates")
    print("  ✗ Serial debugging - use sequential for easier troubleshooting")
    print("  ✗ Limited memory - better to run fewer workers")
    print()

    print("="*80)
    print()


def example_workflow():
    """Show complete workflow with parallel HMBE."""

    print("="*80)
    print("COMPLETE PARALLEL HMBE WORKFLOW EXAMPLE")
    print("="*80)
    print()

    print("This example shows the complete code for a parallel HMBE calculation.")
    print("(Note: Requires QCEngine + Psi4 to actually execute)")
    print()

    example_code = '''
#!/usr/bin/env python
"""Example: Parallel HMBE calculation on 12-water cluster"""

import logging
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models.v1 import BsseEnum, ManyBodyInput, ManyBodyKeywords
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification

# Enable logging to monitor progress
logging.basicConfig(level=logging.INFO)

# 1. Create molecule (12 waters)
mol = Molecule(...)  # Your geometry here

# 2. Define hierarchy (3 groups × 4 waters)
hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers={
        1: ("G0", "W0"), 2: ("G0", "W1"), 3: ("G0", "W2"), 4: ("G0", "W3"),
        5: ("G1", "W4"), 6: ("G1", "W5"), 7: ("G1", "W6"), 8: ("G1", "W7"),
        9: ("G2", "W8"), 10: ("G2", "W9"), 11: ("G2", "W10"), 12: ("G2", "W11"),
    }
)

# 3. Create HMBE specification
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3),
    hierarchy=hierarchy
)

# 4. Create input
keywords = ManyBodyKeywords(
    max_nbody=3,
    bsse_type=[BsseEnum.cp],
    hmbe_spec=hmbe_spec
)

mb_input = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": keywords,
        "specification": {
            "mp2/cc-pvdz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"},
                "driver": "energy"
            }
        }
    }
)

# 5. Run with parallel execution (4 workers)
result = ManyBodyComputer.from_manybodyinput(
    mb_input,
    parallel=True,
    n_workers=4
)

# 6. Get results
print(f"HMBE Energy: {result.properties.return_result['mp2/cc-pvdz']:.10f} Ha")
print(f"Speedup info: {result.properties.hmbe_metadata}")
'''

    print(example_code)
    print()
    print("="*80)
    print()


def main():
    """Run complete parallel tutorial."""

    print()
    print("#" * 80)
    print("# Tutorial 3: Parallel Execution for HMBE")
    print("#" * 80)
    print()

    # Show different API methods
    demonstrate_parallel_api()

    # Explain performance
    explain_performance_expectations()

    # Show monitoring
    demonstrate_monitoring()

    # Best practices
    demonstrate_best_practices()

    # Complete example
    example_workflow()

    # Summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print()

    print("KEY TAKEAWAYS:")
    print("  1. Parallel execution is easy: just set parallel=True")
    print("  2. Use 4-8 workers for best efficiency")
    print("  3. HMBE + Parallel = multiplicative speedup")
    print("  4. Monitor progress with logging")
    print("  5. Consider memory limits when choosing n_workers")
    print()

    print("EXPECTED PERFORMANCE:")
    print("  • Small systems (10-20 fragments):")
    print("    - HMBE reduction: ~1.5-2x")
    print("    - Parallel speedup: ~2-3x (4 workers)")
    print("    - Combined: ~3-6x total speedup")
    print()
    print("  • Medium systems (50 fragments):")
    print("    - HMBE reduction: ~5-10x")
    print("    - Parallel speedup: ~3-4x (4 workers)")
    print("    - Combined: ~15-40x total speedup")
    print()
    print("  • Large systems (100+ fragments):")
    print("    - HMBE reduction: ~50-1000x")
    print("    - Parallel speedup: ~3-4x (4 workers)")
    print("    - Combined: ~150-4000x total speedup!")
    print()

    print("NEXT STEPS:")
    print("  - Try parallel execution on your system")
    print("  - Experiment with different n_workers")
    print("  - Monitor performance and adjust")
    print("  - For >100 fragment systems: await Phase 2 (Direct Enumeration)")
    print()

    print("Tutorial complete!")
    print()


if __name__ == "__main__":
    main()
