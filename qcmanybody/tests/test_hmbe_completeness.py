"""
Unit tests for HMBE completeness validation and sub-cluster generation.

These tests verify that the completeness requirement for Möbius inversion
is properly enforced when using HMBE with Schengen terms.
"""

import pytest
from qcmanybody.hmbe_filter import generate_all_subclusters, select_schengen_terms
from qcelemental.models import Molecule


def test_generate_all_subclusters_single():
    """Test sub-cluster generation for single cluster."""
    result = generate_all_subclusters({(1, 2, 3)})
    expected = {(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)}
    assert result == expected


def test_generate_all_subclusters_multiple():
    """Test sub-cluster generation for multiple clusters."""
    result = generate_all_subclusters({(1, 5, 9), (2, 3)})
    expected = {
        (1,), (2,), (3,), (5,), (9,),
        (1, 5), (1, 9), (5, 9), (2, 3),
        (1, 5, 9)
    }
    assert result == expected


def test_generate_all_subclusters_empty():
    """Test sub-cluster generation for empty set."""
    result = generate_all_subclusters(set())
    assert result == set()


def test_generate_all_subclusters_monomers_only():
    """Test sub-cluster generation for monomers (no sub-clusters)."""
    result = generate_all_subclusters({(1,), (2,), (3,)})
    expected = {(1,), (2,), (3,)}
    assert result == expected


def test_generate_all_subclusters_mixed():
    """Test sub-cluster generation for mixed cluster sizes."""
    result = generate_all_subclusters({(1,), (2, 3), (4, 5, 6)})
    expected = {
        (1,),
        (2,), (3,), (2, 3),
        (4,), (5,), (6,), (4, 5), (4, 6), (5, 6), (4, 5, 6)
    }
    assert result == expected


def test_select_schengen_terms_returns_tuple():
    """Test that select_schengen_terms returns a tuple of (selected, required_subs)."""
    from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification

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

    # Create hierarchy
    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
            5: ("G2", "F4"), 6: ("G2", "F5"),
        },
        tier_names=("group", "fragment")
    )

    # Create HMBE spec with Schengen
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.5,
            distance_metric="R2"
        )
    )

    # Test with candidate that spans 3 groups
    candidates = {(1, 3, 5)}

    result = select_schengen_terms(candidates, mol, hmbe_spec)

    # Should return tuple
    assert isinstance(result, tuple)
    assert len(result) == 2

    selected, required_subs = result

    # Selected should contain the candidate
    assert (1, 3, 5) in selected

    # Required subs should contain all proper sub-clusters
    expected_subs = {(1,), (3,), (5,), (1, 3), (1, 5), (3, 5)}
    assert required_subs == expected_subs


def test_select_schengen_terms_disabled():
    """Test that disabled Schengen returns empty sets."""
    from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification

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

    assert result == (set(), set())


def test_select_schengen_terms_empty_candidates():
    """Test that empty candidates returns empty sets."""
    from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification

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

    assert result == (set(), set())


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
