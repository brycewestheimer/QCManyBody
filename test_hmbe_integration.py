#!/usr/bin/env python3
"""
Integration test for HMBE term enumeration in ManyBodyCore.
Tests that terms are properly counted and completeness is validated.
"""

import sys
from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification


def test_hmbe_term_counts():
    """Test that HMBE properly reduces term counts vs full MBE."""
    print("Testing HMBE term count reduction...")

    # Create a 6-fragment molecule (3 groups of 2 fragments each)
    mol = Molecule(
        symbols=["H"] * 6,
        geometry=[
            [0.0, 0.0, 0.0],  # Frag 1, Group 0
            [1.0, 0.0, 0.0],  # Frag 2, Group 0
            [3.0, 0.0, 0.0],  # Frag 3, Group 1
            [4.0, 0.0, 0.0],  # Frag 4, Group 1
            [6.0, 0.0, 0.0],  # Frag 5, Group 2
            [7.0, 0.0, 0.0],  # Frag 6, Group 2
        ],
        fragments=[[0], [1], [2], [3], [4], [5]],
    )

    # Create 2-tier hierarchy: 3 groups, 6 fragments
    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    # (2,3)-HMBE: max 2 groups, max 3-body
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        enumeration_mode="filter"  # Use filter mode for testing
    )

    # Create ManyBodyCore WITHOUT HMBE (full MBE)
    mbc_mbe = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=None  # No HMBE
    )

    # Create ManyBodyCore WITH HMBE
    mbc_hmbe = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    # Get compute maps
    compute_map_mbe = mbc_mbe.compute_map
    compute_map_hmbe = mbc_hmbe.compute_map

    # Count terms
    mc = "hf/sto-3g"
    mbe_terms = set()
    for bsse_dict in compute_map_mbe[mc].values():
        for terms in bsse_dict.values():
            mbe_terms.update(terms)

    hmbe_terms = set()
    for bsse_dict in compute_map_hmbe[mc].values():
        for terms in bsse_dict.values():
            hmbe_terms.update(terms)

    print(f"  Full MBE terms: {len(mbe_terms)}")
    print(f"  HMBE terms: {len(hmbe_terms)}")
    print(f"  Reduction factor: {len(mbe_terms) / len(hmbe_terms):.2f}x")

    # For (2,3)-HMBE with 6 fragments:
    # - 1-body: 6 terms (all pass)
    # - 2-body: should include only pairs from same group or across 2 groups
    # - 3-body: should exclude terms spanning 3 groups

    # Full MBE should have: 6 + 15 + 20 = 41 terms (1b + 2b + 3b)
    # HMBE should significantly reduce 3-body terms
    assert len(hmbe_terms) < len(mbe_terms), "HMBE should reduce term count"
    print("  ✓ HMBE reduces term count")

    # Verify HMBE statistics
    stats = mbc_hmbe.get_hmbe_statistics()
    assert stats is not None, "HMBE stats should be available"
    print(f"  ✓ HMBE stats: {stats['hmbe_term_counts']}")

    print("✓ HMBE term count test passed!\n")


def test_hmbe_completeness_validation():
    """Test that completeness validation catches missing sub-clusters."""
    print("Testing HMBE completeness validation...")

    # Create a simple 4-fragment molecule
    mol = Molecule(
        symbols=["H"] * 4,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
        },
        tier_names=("group", "fragment")
    )

    # (2,2)-HMBE: max 2 groups, max 2-body
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        enumeration_mode="filter"
    )

    # Create ManyBodyCore with HMBE
    mbc = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    # Extract fragment tuples from compute map
    mc = "hf/sto-3g"
    fragment_tuples = set()
    for bsse_dict in mbc.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                fragment_tuples.add(frag)

    print(f"  Fragment tuples in HMBE: {sorted(fragment_tuples)}")

    # Manually validate completeness
    from itertools import combinations
    missing = []
    for frag_tuple in fragment_tuples:
        n = len(frag_tuple)
        if n == 1:
            continue
        for sub_size in range(1, n):
            for sub_cluster in combinations(frag_tuple, sub_size):
                if sub_cluster not in fragment_tuples:
                    missing.append((frag_tuple, sub_cluster))

    if missing:
        print(f"  ✗ Found {len(missing)} missing sub-clusters:")
        for cluster, sub in missing[:5]:
            print(f"    {cluster} missing {sub}")
        raise AssertionError("Completeness validation failed!")
    else:
        print("  ✓ All clusters have their sub-clusters present")

    # Test the validation method
    try:
        mbc._validate_hmbe_completeness(fragment_tuples, mc, "nocp")
        print("  ✓ _validate_hmbe_completeness() passed")
    except RuntimeError as e:
        print(f"  ✗ Validation failed: {e}")
        raise

    print("✓ Completeness validation test passed!\n")


def test_hmbe_with_schengen():
    """Test HMBE with Schengen terms to ensure sub-clusters are added."""
    print("Testing HMBE with Schengen terms...")

    # Create a 6-fragment molecule
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

    # (2,3)-HMBE with Schengen enabled
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.3,  # Select 30% of excluded terms
            distance_metric="R2"
        ),
        enumeration_mode="filter"
    )

    # Create ManyBodyCore with HMBE + Schengen
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
    fragment_tuples = set()
    for bsse_dict in mbc.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                fragment_tuples.add(frag)

    print(f"  Total terms with Schengen: {len(fragment_tuples)}")

    # Validate completeness (should pass with Schengen sub-cluster addition)
    try:
        mbc._validate_hmbe_completeness(fragment_tuples, mc, "nocp")
        print("  ✓ Completeness validated with Schengen sub-clusters")
    except RuntimeError as e:
        print(f"  ✗ Completeness validation failed: {e}")
        raise

    print("✓ HMBE with Schengen test passed!\n")


def test_hmbe_statistics():
    """Test that HMBE statistics are correctly computed."""
    print("Testing HMBE statistics computation...")

    mol = Molecule(
        symbols=["H"] * 4,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
        },
        tier_names=("group", "fragment")
    )

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        enumeration_mode="auto"
    )

    mbc = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.cp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    stats = mbc.get_hmbe_statistics()
    assert stats is not None, "Stats should not be None"
    assert "mbe_term_counts" in stats, "Missing mbe_term_counts"
    assert "hmbe_term_counts" in stats, "Missing hmbe_term_counts"
    assert "reduction_factors" in stats, "Missing reduction_factors"
    assert "truncation_orders" in stats, "Missing truncation_orders"
    assert "enumeration_mode" in stats, "Missing enumeration_mode"

    print(f"  MBE counts: {stats['mbe_term_counts']}")
    print(f"  HMBE counts: {stats['hmbe_term_counts']}")
    print(f"  Reduction: {stats['reduction_factors']}")
    print(f"  Truncation orders: {stats['truncation_orders']}")
    print(f"  Enumeration mode: {stats['enumeration_mode']} -> {stats['actual_enumeration_mode']}")
    print("  ✓ All statistics present")

    print("✓ HMBE statistics test passed!\n")


if __name__ == "__main__":
    try:
        test_hmbe_term_counts()
        test_hmbe_completeness_validation()
        test_hmbe_with_schengen()
        test_hmbe_statistics()

        print("=" * 60)
        print("ALL INTEGRATION TESTS PASSED!")
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
