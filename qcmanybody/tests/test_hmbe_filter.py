"""
Unit tests for HMBE filtering logic.

Tests the core HMBE filtering algorithm, data structures, and Schengen term selection.
"""

import pytest
import numpy as np
from qcelemental.models import Molecule

from qcmanybody.models.hierarchy import (
    FragmentHierarchy,
    HMBESpecification,
    SchengenSpecification,
)
from qcmanybody.hmbe_filter import (
    passes_hmbe_filter,
    filter_compute_list,
    get_schengen_candidates,
    select_schengen_terms,
    compute_distance_metric,
)


class TestFragmentHierarchy:
    """Tests for FragmentHierarchy data structure."""

    def test_2tier_hierarchy_creation(self):
        """Test creating a valid 2-tier hierarchy."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("GroupA", "SubgroupA1"),
                2: ("GroupA", "SubgroupA2"),
                3: ("GroupB", "SubgroupB1"),
                4: ("GroupB", "SubgroupB2"),
            },
            tier_names=("group", "subgroup"),
        )

        assert hierarchy.num_tiers == 2
        assert hierarchy.fragment_tiers[1] == ("GroupA", "SubgroupA1")
        assert hierarchy.tier_names == ("group", "subgroup")

    def test_3tier_hierarchy_creation(self):
        """Test creating a valid 3-tier hierarchy."""
        hierarchy = FragmentHierarchy(
            num_tiers=3,
            fragment_tiers={
                1: ("Domain1", "Helix1", "Residue1"),
                2: ("Domain1", "Helix1", "Residue2"),
            },
        )

        assert hierarchy.num_tiers == 3
        assert len(hierarchy.fragment_tiers[1]) == 3

    def test_tier_depth_validation(self):
        """Test that mismatched tier depths raise error."""
        with pytest.raises(ValueError, match="expected 2"):
            FragmentHierarchy(
                num_tiers=2,
                fragment_tiers={
                    1: ("GroupA", "SubgroupA1"),
                    2: ("GroupB",),  # Missing tier-2!
                },
            )

    def test_tier_names_validation(self):
        """Test that tier_names length must match num_tiers."""
        with pytest.raises(ValueError, match="expected 2"):
            FragmentHierarchy(
                num_tiers=2,
                fragment_tiers={1: ("A", "A1"), 2: ("B", "B1")},
                tier_names=("tier1",),  # Only 1 name, need 2
            )

    def test_get_tier_groups(self):
        """Test get_tier_groups method."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("GroupA", "SubA1"),
                2: ("GroupA", "SubA2"),
                3: ("GroupB", "SubB1"),
                4: ("GroupB", "SubB2"),
            },
        )

        # Get tier-1 (coarse) groups
        tier1_groups = hierarchy.get_tier_groups(0)
        assert tier1_groups["GroupA"] == frozenset({1, 2})
        assert tier1_groups["GroupB"] == frozenset({3, 4})

        # Get tier-2 (fine) groups
        tier2_groups = hierarchy.get_tier_groups(1)
        assert tier2_groups["SubA1"] == frozenset({1})
        assert tier2_groups["SubA2"] == frozenset({2})

    def test_count_distinct_groups(self):
        """Test count_distinct_groups_at_tier method."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("GroupA", "SubA1"),
                2: ("GroupA", "SubA2"),
                3: ("GroupB", "SubB1"),
            },
        )

        # Fragments 1 and 2 both in GroupA
        assert hierarchy.count_distinct_groups_at_tier((1, 2), 0) == 1

        # Fragments 1 and 3 in different tier-1 groups
        assert hierarchy.count_distinct_groups_at_tier((1, 3), 0) == 2

        # All 3 fragments at tier-2 level
        assert hierarchy.count_distinct_groups_at_tier((1, 2, 3), 1) == 3


class TestHMBESpecification:
    """Tests for HMBESpecification validation."""

    def test_valid_23_hmbe(self):
        """Test valid (2,3)-HMBE specification."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("A", "A2")},
        )

        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy)

        assert spec.truncation_orders == (2, 3)
        assert spec.max_nbody == 3
        assert spec.num_tiers == 2
        assert not spec.is_standard_mbe()

    def test_valid_33_mbe(self):
        """Test that (3,3)-HMBE is recognized as standard MBE."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("A", "A2")},
        )

        spec = HMBESpecification(truncation_orders=(3, 3), hierarchy=hierarchy)

        assert spec.is_standard_mbe()  # All orders equal

    def test_monotonicity_validation(self):
        """Test that non-monotonic truncation orders are rejected."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("A", "A2")},
        )

        # T_1 > T_2 should fail
        with pytest.raises(ValueError, match="non-decreasing"):
            HMBESpecification(truncation_orders=(3, 2), hierarchy=hierarchy)

    def test_truncation_order_count_validation(self):
        """Test that truncation orders must match num_tiers."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("A", "A2")},
        )

        # 3 truncation orders for 2-tier hierarchy
        with pytest.raises(ValueError, match="must match"):
            HMBESpecification(truncation_orders=(2, 3, 4), hierarchy=hierarchy)


class TestPassesHMBEFilter:
    """Tests for the core HMBE filtering function."""

    def setup_method(self):
        """Create test hierarchies for filtering tests."""
        # 4x4 = 16 fragments, 4 tier-1 groups
        self.fragment_tiers_4x4 = {}
        for i in range(4):
            for j in range(4):
                frag_id = i * 4 + j + 1  # 1-indexed
                self.fragment_tiers_4x4[frag_id] = (f"G{i}", f"G{i}_S{j}")

        self.hierarchy_4x4 = FragmentHierarchy(
            num_tiers=2, fragment_tiers=self.fragment_tiers_4x4
        )

    def test_23_hmbe_2body_same_tier1_group(self):
        """Test (2,3)-HMBE: 2-body from same tier-1 group should PASS."""
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=self.hierarchy_4x4)

        # Fragments 1 and 2 both in G0
        assert passes_hmbe_filter((1, 2), spec) is True

    def test_23_hmbe_2body_different_tier1_groups(self):
        """Test (2,3)-HMBE: 2-body from 2 different tier-1 groups should PASS."""
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=self.hierarchy_4x4)

        # Fragments 1 (G0) and 5 (G1): 2 tier-1 groups <= T_1=2
        assert passes_hmbe_filter((1, 5), spec) is True

    def test_23_hmbe_3body_from_3_tier1_groups(self):
        """Test (2,3)-HMBE: 3-body from 3 tier-1 groups should FAIL."""
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=self.hierarchy_4x4)

        # Fragments 1 (G0), 5 (G1), 9 (G2): 3 tier-1 groups > T_1=2
        assert passes_hmbe_filter((1, 5, 9), spec) is False

    def test_23_hmbe_3body_from_2_tier1_groups(self):
        """Test (2,3)-HMBE: 3-body from 2 tier-1 groups should PASS if <= 3 tier-2."""
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=self.hierarchy_4x4)

        # Fragments 1 (G0_S0), 2 (G0_S1), 5 (G1_S0): 2 tier-1, 3 tier-2 <= T_2=3
        assert passes_hmbe_filter((1, 2, 5), spec) is True

    def test_23_hmbe_4body_exceeds_max(self):
        """Test (2,3)-HMBE: 4-body exceeds T_K=3 should FAIL."""
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=self.hierarchy_4x4)

        # Any 4-body term fails n > T_K check
        assert passes_hmbe_filter((1, 2, 3, 4), spec) is False

    def test_33_hmbe_same_as_mbe3(self):
        """Test (3,3)-HMBE: should behave like MBE-3."""
        spec = HMBESpecification(truncation_orders=(3, 3), hierarchy=self.hierarchy_4x4)

        # All 3-body terms should pass
        assert passes_hmbe_filter((1, 5, 9), spec) is True  # 3 tier-1 groups, ok for T_1=3
        assert passes_hmbe_filter((1, 2, 3), spec) is True

        # 4-body should still fail (exceeds max_nbody)
        assert passes_hmbe_filter((1, 2, 3, 4), spec) is False

    def test_24_hmbe_4body_from_2_tier1_groups(self):
        """Test (2,4)-HMBE: 4-body from 2 tier-1 groups can PASS."""
        spec = HMBESpecification(truncation_orders=(2, 4), hierarchy=self.hierarchy_4x4)

        # Fragments 1, 2, 3, 4 all from G0: 1 tier-1 group <= T_1=2
        # 4 tier-2 groups <= T_2=4
        assert passes_hmbe_filter((1, 2, 3, 4), spec) is True


class TestFilterComputeList:
    """Tests for filter_compute_list function."""

    def test_filter_basic_compute_list(self):
        """Test filtering a basic compute list."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "G0_0"),
                2: ("G0", "G0_1"),
                3: ("G1", "G1_0"),
                4: ("G1", "G1_1"),
            },
        )
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy)

        # Mock compute list (fragment_tuple, basis_tuple)
        compute_list = {
            2: {
                ((1, 2), (1, 2, 3, 4)),  # 2-body, 1 tier-1 group -> KEEP
                ((1, 3), (1, 2, 3, 4)),  # 2-body, 2 tier-1 groups -> KEEP
            },
            3: {
                ((1, 2, 3), (1, 2, 3, 4)),  # 3-body, 2 tier-1 groups -> KEEP
                ((1, 3, 4), (1, 2, 3, 4)),  # 3-body, 2 tier-1 groups -> KEEP
                # Would add a 3-tier-1-group term here if we had 6+ fragments
            },
        }

        filtered = filter_compute_list(compute_list, spec)

        # All terms should pass for this simple case
        assert 2 in filtered
        assert len(filtered[2]) == 2
        assert 3 in filtered
        assert len(filtered[3]) == 2

    def test_filter_removes_forbidden_terms(self):
        """Test that filter correctly removes forbidden terms."""
        # 6 fragments: G0(1,2), G1(3,4), G2(5,6)
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G0", "S1"),
                3: ("G1", "S2"),
                4: ("G1", "S3"),
                5: ("G2", "S4"),
                6: ("G2", "S5"),
            },
        )
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy)

        compute_list = {
            3: {
                ((1, 2, 3), (1, 2, 3, 4, 5, 6)),  # 2 tier-1 groups -> KEEP
                ((1, 3, 5), (1, 2, 3, 4, 5, 6)),  # 3 tier-1 groups -> REJECT
            }
        }

        filtered = filter_compute_list(compute_list, spec)

        assert len(filtered[3]) == 1
        assert ((1, 2, 3), (1, 2, 3, 4, 5, 6)) in filtered[3]
        assert ((1, 3, 5), (1, 2, 3, 4, 5, 6)) not in filtered[3]


class TestSchengenSelection:
    """Tests for Schengen term selection."""

    def test_schengen_disabled_returns_empty(self):
        """Test that disabled Schengen returns no terms."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("B", "B1")},
        )
        spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=SchengenSpecification(enabled=False),
        )

        candidates = {(1, 2)}
        mol = Molecule(
            symbols=["He", "He"],
            geometry=[[0, 0, 0], [0, 0, 3]],
            fragments=[[0], [1]],
        )

        selected, required_subs = select_schengen_terms(candidates, mol, spec)
        assert len(selected) == 0
        assert len(required_subs) == 0

    def test_schengen_selection_fraction(self):
        """Test that Schengen respects selection_fraction."""
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                i + 1: (f"G{i // 2}", f"S{i}") for i in range(10)
            },  # 10 fragments
        )
        spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=SchengenSpecification(
                enabled=True, selection_fraction=0.2, distance_metric="R2"
            ),
        )

        # Create molecule with 10 He atoms
        geometry = [[float(i), 0, 0] for i in range(10)]
        mol = Molecule(
            symbols=["He"] * 10, geometry=geometry, fragments=[[i] for i in range(10)]
        )

        # 10 candidates, 20% = 2 terms
        candidates = {tuple(range(1, 11))[i : i + 2] for i in range(8)}  # 8 candidates
        selected, required_subs = select_schengen_terms(candidates, mol, spec)

        expected_count = max(1, int(0.2 * len(candidates)))  # 20% of 8 = 1.6 -> 1
        assert len(selected) == expected_count
        # required_subs should contain sub-clusters of selected terms
        # For 2-body terms, each has 2 1-body sub-clusters
        assert len(required_subs) > 0

    def test_compute_distance_metric_r2(self):
        """Test R2 distance metric calculation."""
        mol = Molecule(
            symbols=["He", "He", "He"],
            geometry=[[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            fragments=[[0], [1], [2]],
        )

        # Fragments 1, 2: distance = 1
        # R2 = 1^2 = 1
        dist = compute_distance_metric((1, 2), mol, "R2")
        assert abs(dist - 1.0) < 1e-6

        # Fragments 1, 2, 3: 3 pairwise distances
        # (1,2): 1, (1,3): 1, (2,3): sqrt(2)
        # R2 = 1 + 1 + 2 = 4
        dist = compute_distance_metric((1, 2, 3), mol, "R2")
        assert abs(dist - 4.0) < 1e-6

    def test_compute_distance_metric_r(self):
        """Test R distance metric calculation."""
        mol = Molecule(
            symbols=["He", "He"],
            geometry=[[0, 0, 0], [3, 0, 0]],
            fragments=[[0], [1]],
        )

        dist = compute_distance_metric((1, 2), mol, "R")
        assert abs(dist - 3.0) < 1e-6

    def test_distance_metric_unknown(self):
        """Test that unknown metric raises error."""
        mol = Molecule(
            symbols=["He", "He"],
            geometry=[[0, 0, 0], [1, 0, 0]],
            fragments=[[0], [1]],
        )

        with pytest.raises(ValueError, match="Unknown distance metric"):
            compute_distance_metric((1, 2), mol, "unknown_metric")


class TestGetSchengenCandidates:
    """Tests for get_schengen_candidates function."""

    def test_get_candidates_basic(self):
        """Test identifying Schengen candidates."""
        # Full MBE list
        mbe_list = {
            2: {((1, 2), (1, 2)), ((1, 3), (1, 2, 3))},
            3: {((1, 2, 3), (1, 2, 3))},
        }

        # Base HMBE list (subset of MBE)
        hmbe_list = {
            2: {((1, 2), (1, 2))},  # Only one 2-body term kept
            # 3-body completely excluded
        }

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={1: ("A", "A1"), 2: ("B", "B1"), 3: ("C", "C1")},
        )
        spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy)

        candidates = get_schengen_candidates(mbe_list, hmbe_list, spec)

        # Should find (1,3) at 2-body and (1,2,3) at 3-body
        assert (1, 3) in candidates
        assert (1, 2, 3) in candidates
        assert (1, 2) not in candidates  # This one is in base HMBE


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
