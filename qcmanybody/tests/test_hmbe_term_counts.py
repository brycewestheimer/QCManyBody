"""
Comprehensive tests for HMBE term enumeration correctness.

These tests verify that HMBE enumeration (both direct and filter modes) produces
the correct number of unique fragment tuples without running actual QM calculations.

Tests cover:
- Common HMBE configurations: (2,3), (2,4), (3,4), (2,3,4)
- Both 2-tier and 3-tier hierarchies
- Comparison against full MBE term counts
- Verification that direct and filter modes produce identical results
- Different system sizes and hierarchy organizations
"""

import pytest
from typing import Set, Tuple, Dict
from qcelemental.models import Molecule

from qcmanybody import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.hmbe_enumerate import enumerate_hmbe_terms, estimate_mbe_term_count
from qcmanybody.builder import build_nbody_compute_list


def create_uniform_2tier_hierarchy(num_fragments: int, frags_per_group: int) -> FragmentHierarchy:
    """Create a uniform 2-tier hierarchy with equal-sized groups.

    Args:
        num_fragments: Total number of fragments
        frags_per_group: Number of fragments per tier-1 group

    Returns:
        FragmentHierarchy with uniform distribution
    """
    fragment_tiers = {}
    for i in range(num_fragments):
        frag_id = i + 1
        tier1_group = f"G{i // frags_per_group}"
        tier2_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_name)

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("group", "fragment")
    )


def create_uniform_3tier_hierarchy(num_fragments: int, tier1_size: int, tier2_size: int) -> FragmentHierarchy:
    """Create a uniform 3-tier hierarchy.

    Args:
        num_fragments: Total number of fragments
        tier1_size: Number of fragments per tier-1 group
        tier2_size: Number of fragments per tier-2 group (within tier-1)

    Returns:
        FragmentHierarchy with 3 tiers
    """
    fragment_tiers = {}
    for i in range(num_fragments):
        frag_id = i + 1
        tier1_group = f"T1_{i // tier1_size}"
        tier2_group = f"T2_{i // tier2_size}"
        tier3_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_group, tier3_name)

    return FragmentHierarchy(
        num_tiers=3,
        fragment_tiers=fragment_tiers,
        tier_names=("domain", "subdomain", "fragment")
    )


def extract_fragment_tuples_from_compute_map(compute_map: Dict) -> Set[Tuple[int, ...]]:
    """Extract all unique fragment tuples from a compute map.

    Args:
        compute_map: Compute map from ManyBodyCore

    Returns:
        Set of unique fragment tuples
    """
    all_frags = set()

    # Only iterate over the consolidated 'all' set to avoid redundant
    # scanning of identical data across bsse subdicts (e.g., 'all' vs 'nocp').
    for mc_dict in compute_map.values():
        bsse_all = mc_dict.get("all", {})
        for nbody_set in bsse_all.values():
            for frag_tuple, _ in nbody_set:
                all_frags.add(frag_tuple)

    return all_frags


class TestHMBETermCounts2Tier:
    """Test HMBE term counts for 2-tier hierarchies."""

    def test_23_hmbe_water6_3x2(self):
        """Test (2,3)-HMBE on 6-fragment 3×2 system."""
        # Setup: 6 fragments in 3 groups of 2
        mol = Molecule(
            symbols=["He"] * 6,
            geometry=[[float(i), 0.0, 0.0] for i in range(6)],
            fragments=[[i] for i in range(6)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=6, frags_per_group=2)
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        # Test filter mode
        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy, enumeration_mode="filter"),
        )

        # Test direct mode
        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy, enumeration_mode="direct"),
        )

        # Extract fragment tuples
        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        # Both modes should produce identical tuples
        assert filter_tuples == direct_tuples, "Filter and direct modes must produce identical fragment sets"

        # Calculate expected counts
        # MBE-3 for 6 fragments: C(6,1) + C(6,2) + C(6,3) = 6 + 15 + 20 = 41
        mbe_count = estimate_mbe_term_count(num_fragments=6, max_nbody=3)
        assert mbe_count == 41

        # HMBE should have fewer terms
        hmbe_count = len(filter_tuples)
        assert hmbe_count < mbe_count, "HMBE should reduce term count"
        assert hmbe_count > 0, "HMBE should have some terms"

        # Expected exact count for (2,3)-HMBE with 3 groups of 2:
        # 1-body: all 6 monomers → 6
        # 2-body: all pairs (constraint always satisfied) → C(6,2)=15
        # 3-body: choose any two groups (C(3,2)=3) and any 3 of the 4 fragments within them (C(4,3)=4) → 3*4=12
        # Total = 6 + 15 + 12 = 33
        expected_exact = 33
        assert hmbe_count == expected_exact, f"Expected {expected_exact} terms, got {hmbe_count}"

        # Verify statistics
        stats = mbcore_filter.get_hmbe_statistics()
        assert stats["hmbe_term_counts"]["hf/sto-3g"] == hmbe_count
        assert stats["mbe_term_counts"]["hf/sto-3g"] == mbe_count

    def test_24_hmbe_water8_4x2(self):
        """Test (2,4)-HMBE on 8-fragment 4×2 system."""
        mol = Molecule(
            symbols=["He"] * 8,
            geometry=[[float(i), 0.0, 0.0] for i in range(8)],
            fragments=[[i] for i in range(8)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=8, frags_per_group=2)

        # Test both modes
        hmbe_spec_filter = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=hierarchy,
            enumeration_mode="filter"
        )

        hmbe_spec_direct = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=hierarchy,
            enumeration_mode="direct"
        )

        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_filter,
        )

        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_direct,
        )

        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        # Verify modes agree
        assert filter_tuples == direct_tuples

        # MBE-4 for 8 fragments: C(8,1) + C(8,2) + C(8,3) + C(8,4) = 8 + 28 + 56 + 70 = 162
        mbe_count = estimate_mbe_term_count(num_fragments=8, max_nbody=4)
        assert mbe_count == 162

        hmbe_count = len(filter_tuples)
        assert hmbe_count < mbe_count, "HMBE should reduce term count"

        # Verify reduction factor is reasonable
        reduction = mbe_count / hmbe_count
        assert reduction > 1.5, f"Expected significant reduction, got {reduction:.2f}x"

    def test_34_hmbe_water12_4x3(self):
        """Test (3,4)-HMBE on 12-fragment 4×3 system."""
        mol = Molecule(
            symbols=["He"] * 12,
            geometry=[[float(i % 4), float(i // 4), 0.0] for i in range(12)],
            fragments=[[i] for i in range(12)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=12, frags_per_group=3)

        hmbe_spec_filter = HMBESpecification(
            truncation_orders=(3, 4),
            hierarchy=hierarchy,
            enumeration_mode="filter"
        )

        hmbe_spec_direct = HMBESpecification(
            truncation_orders=(3, 4),
            hierarchy=hierarchy,
            enumeration_mode="direct"
        )

        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_filter,
        )

        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_direct,
        )

        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        # Verify modes agree
        assert filter_tuples == direct_tuples

        # MBE-4 for 12 fragments: C(12,1) + C(12,2) + C(12,3) + C(12,4) = 12 + 66 + 220 + 495 = 793
        mbe_count = estimate_mbe_term_count(num_fragments=12, max_nbody=4)
        assert mbe_count == 793

        hmbe_count = len(filter_tuples)
        assert hmbe_count < mbe_count, "HMBE should reduce term count"

        # For (3,4)-HMBE, reduction should be modest but noticeable
        reduction = mbe_count / hmbe_count
        assert reduction > 1.2, f"Expected reduction, got {reduction:.2f}x"


class TestHMBETermCounts3Tier:
    """Test HMBE term counts for 3-tier hierarchies."""

    def test_234_hmbe_8frags_2x2x2(self):
        """Test (2,3,4)-HMBE on 8-fragment 2×2×2 system."""
        mol = Molecule(
            symbols=["He"] * 8,
            geometry=[[float(i % 2), float((i // 2) % 2), float(i // 4)] for i in range(8)],
            fragments=[[i] for i in range(8)]
        )

        hierarchy = create_uniform_3tier_hierarchy(
            num_fragments=8,
            tier1_size=4,  # 2 tier-1 groups
            tier2_size=2   # 4 tier-2 groups
        )

        hmbe_spec_filter = HMBESpecification(
            truncation_orders=(2, 3, 4),
            hierarchy=hierarchy,
            enumeration_mode="filter"
        )

        hmbe_spec_direct = HMBESpecification(
            truncation_orders=(2, 3, 4),
            hierarchy=hierarchy,
            enumeration_mode="direct"
        )

        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_filter,
        )

        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_direct,
        )

        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        # Verify modes agree
        assert filter_tuples == direct_tuples

        # MBE-4 for 8 fragments: 8 + 28 + 56 + 70 = 162
        mbe_count = estimate_mbe_term_count(num_fragments=8, max_nbody=4)
        assert mbe_count == 162

        hmbe_count = len(filter_tuples)
        assert hmbe_count < mbe_count, "3-tier HMBE should reduce term count"

        # 3-tier should give better reduction
        reduction = mbe_count / hmbe_count
        assert reduction > 1.5, f"Expected good reduction for 3-tier, got {reduction:.2f}x"


class TestHMBEVsMBEComparison:
    """Test HMBE vs MBE term count comparisons."""

    @pytest.mark.parametrize("truncation_orders,num_fragments,frags_per_group", [
        ((2, 3), 6, 2),   # (2,3)-HMBE, 6 fragments, 3 groups
        ((2, 3), 8, 2),   # (2,3)-HMBE, 8 fragments, 4 groups
        ((2, 4), 8, 2),   # (2,4)-HMBE, 8 fragments, 4 groups
        ((3, 4), 9, 3),   # (3,4)-HMBE, 9 fragments, 3 groups
        ((2, 3), 12, 3),  # (2,3)-HMBE, 12 fragments, 4 groups
        ((2, 4), 12, 3),  # (2,4)-HMBE, 12 fragments, 4 groups
    ])
    def test_hmbe_reduces_terms(self, truncation_orders, num_fragments, frags_per_group):
        """Test that HMBE reduces term count compared to standard MBE."""
        mol = Molecule(
            symbols=["He"] * num_fragments,
            geometry=[[float(i), 0.0, 0.0] for i in range(num_fragments)],
            fragments=[[i] for i in range(num_fragments)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments, frags_per_group)
        T_K = truncation_orders[-1]  # max_nbody

        hmbe_spec = HMBESpecification(
            truncation_orders=truncation_orders,
            hierarchy=hierarchy,
            enumeration_mode="direct"  # Use direct for efficiency
        )

        levels = {i: "hf/sto-3g" for i in range(1, T_K + 1)}

        # HMBE calculation
        mbcore_hmbe = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        # Standard MBE calculation (no HMBE)
        mbcore_mbe = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=None,
        )

        hmbe_tuples = extract_fragment_tuples_from_compute_map(mbcore_hmbe.compute_map)
        mbe_tuples = extract_fragment_tuples_from_compute_map(mbcore_mbe.compute_map)

        # HMBE should be a subset of MBE
        assert hmbe_tuples.issubset(mbe_tuples), "HMBE terms must be subset of MBE terms"

        # HMBE should have fewer terms (unless T_1 >= num_groups and T_K >= num_fragments)
        num_groups = (num_fragments + frags_per_group - 1) // frags_per_group
        T_1 = truncation_orders[0]

        if T_1 >= num_groups:
            # HMBE is equivalent to MBE in this case
            assert len(hmbe_tuples) == len(mbe_tuples), "When T_1 >= num_groups, HMBE = MBE"
        else:
            # HMBE should reduce terms
            assert len(hmbe_tuples) < len(mbe_tuples), "HMBE should reduce term count when T_1 < num_groups"

            reduction_factor = len(mbe_tuples) / len(hmbe_tuples)
            print(f"\n  {truncation_orders}-HMBE on {num_fragments} fragments: "
                  f"{len(hmbe_tuples)}/{len(mbe_tuples)} terms ({reduction_factor:.2f}x reduction)")


class TestDirectVsFilterModeEquivalence:
    """Test that direct and filter enumeration modes are exactly equivalent."""

    @pytest.mark.parametrize("num_fragments,frags_per_group,truncation_orders", [
        (6, 2, (2, 3)),
        (8, 2, (2, 4)),
        (9, 3, (3, 3)),
        (12, 3, (2, 4)),
        (12, 4, (3, 4)),
    ])
    def test_modes_produce_identical_tuples(self, num_fragments, frags_per_group, truncation_orders):
        """Test that direct and filter modes produce exactly the same fragment tuples."""
        mol = Molecule(
            symbols=["He"] * num_fragments,
            geometry=[[float(i), 0.0, 0.0] for i in range(num_fragments)],
            fragments=[[i] for i in range(num_fragments)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments, frags_per_group)
        T_K = truncation_orders[-1]
        levels = {i: "hf/sto-3g" for i in range(1, T_K + 1)}

        # Filter mode
        hmbe_spec_filter = HMBESpecification(
            truncation_orders=truncation_orders,
            hierarchy=hierarchy,
            enumeration_mode="filter"
        )

        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_filter,
        )

        # Direct mode
        hmbe_spec_direct = HMBESpecification(
            truncation_orders=truncation_orders,
            hierarchy=hierarchy,
            enumeration_mode="direct"
        )

        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_direct,
        )

        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        # Must be exactly identical
        assert filter_tuples == direct_tuples, (
            f"Direct and filter modes produced different tuples for {truncation_orders}-HMBE "
            f"on {num_fragments} fragments: filter={len(filter_tuples)}, direct={len(direct_tuples)}"
        )

        # Verify statistics also match
        stats_filter = mbcore_filter.get_hmbe_statistics()
        stats_direct = mbcore_direct.get_hmbe_statistics()

        assert stats_filter["hmbe_term_counts"] == stats_direct["hmbe_term_counts"]
        assert stats_filter["reduction_factors"] == stats_direct["reduction_factors"]


class TestAutoModeSelection:
    """Test that auto mode selects the appropriate enumeration method."""

    def test_auto_mode_selects_filter_for_small_systems(self):
        """Test that auto mode uses filter for <30 fragments."""
        mol = Molecule(
            symbols=["He"] * 12,
            geometry=[[float(i), 0.0, 0.0] for i in range(12)],
            fragments=[[i] for i in range(12)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=12, frags_per_group=3)

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            enumeration_mode="auto"
        )

        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        stats = mbcore.get_hmbe_statistics()
        assert stats["enumeration_mode"] == "auto"
        assert stats["actual_enumeration_mode"] == "filter", "Auto should select filter for 12 fragments"

    def test_auto_mode_selects_direct_for_large_systems(self):
        """Test that auto mode uses direct for >=30 fragments."""
        mol = Molecule(
            symbols=["He"] * 32,
            geometry=[[float(i % 8), float(i // 8), 0.0] for i in range(32)],
            fragments=[[i] for i in range(32)]
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=32, frags_per_group=4)

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            enumeration_mode="auto"
        )

        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        stats = mbcore.get_hmbe_statistics()
        assert stats["enumeration_mode"] == "auto"
        assert stats["actual_enumeration_mode"] == "direct", "Auto should select direct for 32 fragments"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
