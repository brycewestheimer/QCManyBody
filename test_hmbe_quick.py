#!/usr/bin/env python3
"""
Quick lightweight test for HMBE term enumeration fixes.
Tests key functionality without heavy pytest infrastructure.
"""

import sys
from qcmanybody.hmbe_filter import generate_all_subclusters, select_schengen_terms
from qcelemental.models import Molecule
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification


def test_generate_subclusters():
    """Test sub-cluster generation."""
    print("Testing generate_all_subclusters...")

    # Test 1: Single 3-body cluster
    result = generate_all_subclusters({(1, 2, 3)})
    expected = {(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)}
    assert result == expected, f"Failed: got {result}, expected {expected}"
    print("  ✓ Single 3-body cluster")

    # Test 2: Multiple clusters
    result = generate_all_subclusters({(1, 5, 9), (2, 3)})
    expected = {
        (1,), (2,), (3,), (5,), (9,),
        (1, 5), (1, 9), (5, 9), (2, 3),
        (1, 5, 9)
    }
    assert result == expected, f"Failed: got {result}, expected {expected}"
    print("  ✓ Multiple clusters")

    # Test 3: Empty set
    result = generate_all_subclusters(set())
    assert result == set(), f"Failed: got {result}, expected empty set"
    print("  ✓ Empty set")

    # Test 4: Monomers only (no sub-clusters)
    result = generate_all_subclusters({(1,), (2,), (3,)})
    expected = {(1,), (2,), (3,)}
    assert result == expected, f"Failed: got {result}, expected {expected}"
    print("  ✓ Monomers only")

    print("✓ All sub-cluster tests passed!\n")


def test_schengen_returns_tuple():
    """Test that select_schengen_terms returns (selected, required_subs) tuple."""
    print("Testing select_schengen_terms return signature...")

    # Create simple 6-fragment molecule
    mol = Molecule(
        symbols=["H"] * 6,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3], [4], [5]],
    )

    # Create hierarchy (2-tier, 3 groups of 2 fragments each)
    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    # Create (2,3)-HMBE spec with Schengen enabled
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.5,
            distance_metric="R2"
        )
    )

    # Test with candidate that spans 3 groups (would be excluded by (2,3)-HMBE)
    candidates = {(1, 3, 5)}

    result = select_schengen_terms(candidates, mol, hmbe_spec)

    # Should return tuple
    assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
    assert len(result) == 2, f"Expected tuple of length 2, got {len(result)}"
    print("  ✓ Returns tuple of length 2")

    selected, required_subs = result

    # Selected should contain the candidate
    assert (1, 3, 5) in selected, f"Expected (1,3,5) in selected, got {selected}"
    print("  ✓ Selected contains expected term")

    # Required subs should contain all proper sub-clusters
    expected_subs = {(1,), (3,), (5,), (1, 3), (1, 5), (3, 5)}
    assert required_subs == expected_subs, f"Expected {expected_subs}, got {required_subs}"
    print("  ✓ Required sub-clusters computed correctly")

    print("✓ Schengen return signature test passed!\n")


def test_schengen_disabled():
    """Test that disabled Schengen returns empty sets."""
    print("Testing Schengen disabled...")

    mol = Molecule(
        symbols=["H"] * 3,
        geometry=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        fragments=[[0], [1], [2]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={1: ("G0", "F0"), 2: ("G0", "F1"), 3: ("G1", "F2")},
        tier_names=("group", "fragment")
    )

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(enabled=False, selection_fraction=0.1)
    )

    candidates = {(1, 2, 3)}
    result = select_schengen_terms(candidates, mol, hmbe_spec)

    assert result == (set(), set()), f"Expected (set(), set()), got {result}"
    print("  ✓ Returns empty sets when disabled")
    print("✓ Schengen disabled test passed!\n")


def test_schengen_empty_candidates():
    """Test that empty candidates returns empty sets."""
    print("Testing empty candidates...")

    mol = Molecule(
        symbols=["H"] * 3,
        geometry=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        fragments=[[0], [1], [2]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={1: ("G0", "F0"), 2: ("G0", "F1"), 3: ("G1", "F2")},
        tier_names=("group", "fragment")
    )

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(enabled=True, selection_fraction=0.1)
    )

    candidates = set()
    result = select_schengen_terms(candidates, mol, hmbe_spec)

    assert result == (set(), set()), f"Expected (set(), set()), got {result}"
    print("  ✓ Returns empty sets for empty candidates")
    print("✓ Empty candidates test passed!\n")


if __name__ == "__main__":
    try:
        test_generate_subclusters()
        test_schengen_returns_tuple()
        test_schengen_disabled()
        test_schengen_empty_candidates()

        print("=" * 60)
        print("ALL TESTS PASSED!")
        print("=" * 60)
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
