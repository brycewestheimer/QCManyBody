#!/usr/bin/env python3
"""
Test direct enumeration mode vs filter mode for HMBE.
Ensures both methods produce identical term sets.
"""

import sys
from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.hmbe_enumerate import enumerate_hmbe_terms
from qcmanybody.hmbe_filter import filter_compute_list
from qcmanybody.builder import build_nbody_compute_list


def test_direct_vs_filter_enumeration():
    """Test that direct and filter enumeration produce identical results."""
    print("Testing direct vs filter enumeration...")

    # Create 6-fragment molecule (2-tier hierarchy)
    mol = Molecule(
        symbols=["H"] * 6,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [7.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3], [4], [5]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    # Test with (2,3)-HMBE
    hmbe_spec_filter = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        enumeration_mode="filter"
    )

    hmbe_spec_direct = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        enumeration_mode="direct"
    )

    # Create cores with both modes
    mbc_filter = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec_filter
    )

    mbc_direct = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec_direct
    )

    # Extract fragment tuples from both
    mc = "hf/sto-3g"

    filter_frags = set()
    for bsse_dict in mbc_filter.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                filter_frags.add(frag)

    direct_frags = set()
    for bsse_dict in mbc_direct.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                direct_frags.add(frag)

    print(f"  Filter mode terms: {len(filter_frags)}")
    print(f"  Direct mode terms: {len(direct_frags)}")

    # Check if sets are identical
    if filter_frags == direct_frags:
        print("  ✓ Filter and direct modes produce identical term sets")
    else:
        print("  ✗ MISMATCH between filter and direct modes!")
        only_filter = filter_frags - direct_frags
        only_direct = direct_frags - filter_frags

        if only_filter:
            print(f"    Only in filter ({len(only_filter)}): {sorted(only_filter)[:5]}...")
        if only_direct:
            print(f"    Only in direct ({len(only_direct)}): {sorted(only_direct)[:5]}...")

        raise AssertionError("Direct and filter enumeration produced different term sets")

    # Verify completeness for both
    from itertools import combinations
    for name, frag_set in [("filter", filter_frags), ("direct", direct_frags)]:
        missing = []
        for frag_tuple in frag_set:
            n = len(frag_tuple)
            if n == 1:
                continue
            for sub_size in range(1, n):
                for sub_cluster in combinations(frag_tuple, sub_size):
                    if sub_cluster not in frag_set:
                        missing.append((frag_tuple, sub_cluster))

        if missing:
            print(f"  ✗ Completeness violation in {name} mode: {len(missing)} missing sub-clusters")
            raise AssertionError(f"Completeness violation in {name} mode")

    print("  ✓ Both modes satisfy completeness requirement")
    print("✓ Direct vs filter enumeration test passed!\n")


def test_direct_enumeration_function():
    """Test the standalone enumerate_hmbe_terms function."""
    print("Testing standalone enumerate_hmbe_terms...")

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
    )

    # Call direct enumeration
    hmbe_terms = enumerate_hmbe_terms(hmbe_spec)

    print(f"  Enumerated {len(hmbe_terms)} HMBE terms")

    # Verify these are valid HMBE terms
    from qcmanybody.hmbe_filter import passes_hmbe_filter

    for frag_tuple in hmbe_terms:
        if not passes_hmbe_filter(frag_tuple, hmbe_spec):
            print(f"  ✗ Term {frag_tuple} should not pass HMBE filter!")
            raise AssertionError(f"Invalid HMBE term: {frag_tuple}")

    print("  ✓ All enumerated terms pass HMBE filter")

    # Check completeness
    from itertools import combinations
    missing = []
    for frag_tuple in hmbe_terms:
        n = len(frag_tuple)
        if n == 1:
            continue
        for sub_size in range(1, n):
            for sub_cluster in combinations(frag_tuple, sub_size):
                if sub_cluster not in hmbe_terms:
                    missing.append((frag_tuple, sub_cluster))

    if missing:
        print(f"  ✗ Completeness violation: {len(missing)} missing sub-clusters")
        for cluster, sub in missing[:5]:
            print(f"    {cluster} missing {sub}")
        raise AssertionError("Completeness violation in direct enumeration")

    print("  ✓ Direct enumeration satisfies completeness")
    print("✓ Standalone enumeration test passed!\n")


def test_larger_system_direct():
    """Test direct enumeration on a larger system (12 fragments)."""
    print("Testing direct enumeration on larger system (12 fragments)...")

    # Create 12-fragment molecule (3 groups of 4 fragments each)
    mol = Molecule(
        symbols=["H"] * 12,
        geometry=[[float(i), 0.0, 0.0] for i in range(12)],
        fragments=[[i] for i in range(12)],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"), 3: ("G0", "F2"), 4: ("G0", "F3"),
            5: ("G1", "F4"), 6: ("G1", "F5"), 7: ("G1", "F6"), 8: ("G1", "F7"),
            9: ("G2", "F8"), 10: ("G2", "F9"), 11: ("G2", "F10"), 12: ("G2", "F11"),
        },
        tier_names=("group", "fragment")
    )

    # (2,4)-HMBE: max 2 groups, max 4-body
    hmbe_spec_direct = HMBESpecification(
        truncation_orders=(2, 4),
        hierarchy=hierarchy,
        enumeration_mode="direct"
    )

    # Create core with direct mode (this would timeout with filter for larger systems)
    mbc = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec_direct
    )

    # Get statistics
    stats = mbc.get_hmbe_statistics()
    print(f"  MBE terms: {stats['mbe_term_counts']}")
    print(f"  HMBE terms: {stats['hmbe_term_counts']}")
    print(f"  Reduction: {stats['reduction_factors']}")
    print(f"  Enumeration mode: {stats['actual_enumeration_mode']}")

    # Verify the reduction is significant
    mc = "hf/sto-3g"
    reduction = stats['reduction_factors'][mc]
    assert reduction > 1.0, f"Expected reduction > 1.0x, got {reduction}"
    print(f"  ✓ Achieved {reduction:.1f}x reduction")

    print("✓ Larger system test passed!\n")


def test_direct_enum_with_schengen():
    """Test that direct enumeration works with Schengen terms."""
    print("Testing direct enumeration with Schengen terms...")

    mol = Molecule(
        symbols=["H"] * 6,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [7.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3], [4], [5]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    from qcmanybody.models.hierarchy import SchengenSpecification

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.3,
            distance_metric="R2"
        ),
        enumeration_mode="direct"
    )

    mbc = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    # Extract fragment tuples
    mc = "hf/sto-3g"
    frag_set = set()
    for bsse_dict in mbc.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                frag_set.add(frag)

    print(f"  Total terms with Schengen: {len(frag_set)}")

    # Validate completeness
    try:
        mbc._validate_hmbe_completeness(frag_set, mc, "nocp")
        print("  ✓ Completeness validated with direct + Schengen")
    except RuntimeError as e:
        print(f"  ✗ Completeness validation failed: {e}")
        raise

    print("✓ Direct enumeration with Schengen test passed!\n")


if __name__ == "__main__":
    try:
        test_direct_enumeration_function()
        test_direct_vs_filter_enumeration()
        test_direct_enum_with_schengen()
        test_larger_system_direct()

        print("=" * 60)
        print("ALL DIRECT ENUMERATION TESTS PASSED!")
        print("=" * 60)
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
