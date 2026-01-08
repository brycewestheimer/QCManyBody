"""
Unit tests for direct HMBE enumeration (top-down approach).

These tests verify that the direct enumeration algorithm correctly generates
HMBE-allowed terms without first generating all MBE terms.
"""

import pytest
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.hmbe_enumerate import (
    enumerate_hmbe_terms,
    enumerate_hmbe_terms_2tier,
    enumerate_hmbe_terms_3tier,
    estimate_mbe_term_count,
    estimate_hmbe_reduction_factor,
    get_enumeration_statistics,
)
from qcmanybody.hmbe_filter import passes_hmbe_filter


def create_2tier_4x4_hierarchy():
    """Create a 2-tier 4×4 hierarchy (16 fragments, 4 tier-1 groups)."""
    fragment_tiers = {}
    for i in range(16):
        frag_id = i + 1
        tier1_group = f"G{i // 4}"
        tier2_name = f"W{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_name)

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("group", "water")
    )


def create_3tier_2x2x2_hierarchy():
    """Create a 3-tier 2×2×2 hierarchy (8 fragments)."""
    fragment_tiers = {}
    for i in range(8):
        frag_id = i + 1
        tier1_group = f"T1_{i // 4}"
        tier2_group = f"T2_{i // 2}"
        tier3_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_group, tier3_name)

    return FragmentHierarchy(
        num_tiers=3,
        fragment_tiers=fragment_tiers,
        tier_names=("domain", "subdomain", "fragment")
    )


class TestEnumerateHMBETerms2Tier:
    """Tests for 2-tier direct enumeration."""

    def test_2tier_23_hmbe_16frags(self):
        """Test (2,3)-HMBE on 16-fragment 4×4 system."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Verify terms are generated
        assert len(terms) > 0

        # Verify all terms pass HMBE filter (consistency check)
        for term in terms:
            assert passes_hmbe_filter(term, hmbe_spec), \
                f"Term {term} should pass HMBE filter but doesn't"

        # Verify term counts match expected
        # From water-16 validation: (2,3)-HMBE should have 440 terms
        assert len(terms) == 440, f"Expected 440 terms, got {len(terms)}"

    def test_2tier_24_hmbe_16frags(self):
        """Test (2,4)-HMBE on 16-fragment 4×4 system."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Verify all terms pass HMBE filter
        for term in terms:
            assert passes_hmbe_filter(term, hmbe_spec)

        # From water-16 validation: (2,4)-HMBE should have 852 terms
        assert len(terms) == 852, f"Expected 852 terms, got {len(terms)}"

    def test_2tier_34_hmbe_16frags(self):
        """Test (3,4)-HMBE on 16-fragment 4×4 system."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(3, 4),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Verify all terms pass HMBE filter
        for term in terms:
            assert passes_hmbe_filter(term, hmbe_spec)

        # From water-16 validation: (3,4)-HMBE should have 2260 terms
        assert len(terms) == 2260, f"Expected 2260 terms, got {len(terms)}"

    def test_2tier_single_tier1_group(self):
        """Test edge case: all fragments in single tier-1 group."""
        fragment_tiers = {
            1: ("G0", "F1"),
            2: ("G0", "F2"),
            3: ("G0", "F3"),
            4: ("G0", "F4"),
        }
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers=fragment_tiers
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # With single tier-1 group and T_1=2, only terms from 1 group are allowed
        # This should give all 1-body, 2-body, and 3-body terms
        # C(4,1) + C(4,2) + C(4,3) = 4 + 6 + 4 = 14
        assert len(terms) == 14

        # Verify all are from same tier-1 group
        for term in terms:
            groups = {hierarchy.fragment_tiers[f][0] for f in term}
            assert len(groups) == 1, "All terms should be from single tier-1 group"

    def test_2tier_all_1body_terms(self):
        """Test that all 1-body terms are always included."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Extract 1-body terms
        onebody_terms = {t for t in terms if len(t) == 1}

        # Should have all 16 monomers
        assert len(onebody_terms) == 16
        assert onebody_terms == {(i,) for i in range(1, 17)}

    def test_2tier_no_3group_3body_terms(self):
        """Test that 3-body terms spanning 3 tier-1 groups are excluded for (2,3)-HMBE."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Check 3-body terms
        threebody_terms = [t for t in terms if len(t) == 3]

        for term in threebody_terms:
            tier1_groups = {hierarchy.fragment_tiers[f][0] for f in term}
            # For (2,3)-HMBE, 3-body terms can span at most 2 tier-1 groups
            assert len(tier1_groups) <= 2, \
                f"Term {term} spans {len(tier1_groups)} tier-1 groups (should be ≤2)"


class TestEnumerateHMBETerms3Tier:
    """Tests for 3-tier direct enumeration."""

    def test_3tier_234_hmbe_8frags(self):
        """Test (2,3,4)-HMBE on 8-fragment 2×2×2 system."""
        hierarchy = create_3tier_2x2x2_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3, 4),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_3tier(hmbe_spec)

        # Verify terms are generated
        assert len(terms) > 0

        # Verify all terms pass HMBE filter
        for term in terms:
            assert passes_hmbe_filter(term, hmbe_spec), \
                f"Term {term} should pass HMBE filter but doesn't"

        # Verify all 1-body terms present
        onebody_terms = {t for t in terms if len(t) == 1}
        assert len(onebody_terms) == 8

    def test_3tier_all_tiers_constrained(self):
        """Test that all tier constraints are respected in 3-tier HMBE."""
        hierarchy = create_3tier_2x2x2_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 2, 3),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_3tier(hmbe_spec)

        for term in terms:
            # Count tier-1 groups
            tier1_groups = {hierarchy.fragment_tiers[f][0] for f in term}
            assert len(tier1_groups) <= 2, \
                f"Term {term} violates T_1=2 constraint"

            # Count tier-2 groups
            tier2_groups = {hierarchy.fragment_tiers[f][1] for f in term}
            assert len(tier2_groups) <= 2, \
                f"Term {term} violates T_2=2 constraint"

            # Check max n-body
            assert len(term) <= 3, \
                f"Term {term} violates T_3=3 constraint"


class TestEnumerateHMBETermsDispatcher:
    """Tests for the main dispatcher function."""

    def test_dispatcher_2tier(self):
        """Test dispatcher correctly routes to 2-tier implementation."""
        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        # Use dispatcher
        terms_dispatcher = enumerate_hmbe_terms(hmbe_spec)

        # Use direct 2-tier
        terms_direct = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Should be identical
        assert terms_dispatcher == terms_direct

    def test_dispatcher_3tier(self):
        """Test dispatcher correctly routes to 3-tier implementation."""
        hierarchy = create_3tier_2x2x2_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3, 4),
            hierarchy=hierarchy
        )

        # Use dispatcher
        terms_dispatcher = enumerate_hmbe_terms(hmbe_spec)

        # Use direct 3-tier
        terms_direct = enumerate_hmbe_terms_3tier(hmbe_spec)

        # Should be identical
        assert terms_dispatcher == terms_direct

    def test_dispatcher_invalid_tiers(self):
        """Test dispatcher raises error for unsupported tier counts."""
        # Try to create 4-tier hierarchy (not supported, max is 3)
        with pytest.raises(ValueError, match="ensure this value is less than or equal to 3"):
            hierarchy = FragmentHierarchy(
                num_tiers=4,
                fragment_tiers={
                    1: ("T1", "T2", "T3", "T4"),
                    2: ("T1", "T2", "T3", "T4"),
                }
            )


class TestDirectVsFilterEquivalence:
    """Tests comparing direct enumeration with filter approach."""

    def test_direct_equals_filter_2tier_23(self):
        """Test direct enumeration produces same results as filter for (2,3)-HMBE."""
        from qcmanybody.builder import build_nbody_compute_list
        from qcmanybody.hmbe_filter import filter_compute_list

        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        # Direct enumeration
        direct_terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Filter approach: generate all MBE-3 terms, then filter
        mbe_compute_list = build_nbody_compute_list(
            bsse_type=["nocp"],
            nfragments=16,
            nbodies=[1, 2, 3],
            return_total_data=False,
            supersystem_ie_only=False,
        )

        # Filter to HMBE
        filtered_dict = filter_compute_list(mbe_compute_list["all"], hmbe_spec)

        # Extract fragment tuples from filtered dict
        filter_terms = set()
        for nbody_dict in filtered_dict.values():
            for frag_tuple, _ in nbody_dict:
                filter_terms.add(frag_tuple)

        # Direct and filter should produce identical term sets
        assert direct_terms == filter_terms, \
            f"Direct enumeration ({len(direct_terms)} terms) differs from filter ({len(filter_terms)} terms)"

    def test_direct_equals_filter_2tier_24(self):
        """Test equivalence for (2,4)-HMBE."""
        from qcmanybody.builder import build_nbody_compute_list
        from qcmanybody.hmbe_filter import filter_compute_list

        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=hierarchy
        )

        # Direct enumeration
        direct_terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Filter approach
        mbe_compute_list = build_nbody_compute_list(
            bsse_type=["nocp"],
            nfragments=16,
            nbodies=[1, 2, 3, 4],
            return_total_data=False,
            supersystem_ie_only=False,
        )

        filtered_dict = filter_compute_list(mbe_compute_list["all"], hmbe_spec)

        filter_terms = set()
        for nbody_dict in filtered_dict.values():
            for frag_tuple, _ in nbody_dict:
                filter_terms.add(frag_tuple)

        assert direct_terms == filter_terms


class TestHelperFunctions:
    """Tests for estimation and statistics helper functions."""

    def test_estimate_mbe_term_count(self):
        """Test MBE term count estimation."""
        # For 16 fragments, MBE-3:
        # C(16,1) + C(16,2) + C(16,3) = 16 + 120 + 560 = 696
        count = estimate_mbe_term_count(num_fragments=16, max_nbody=3)
        assert count == 696

        # For 10 fragments, MBE-4:
        # C(10,1) + C(10,2) + C(10,3) + C(10,4) = 10 + 45 + 120 + 210 = 385
        count = estimate_mbe_term_count(num_fragments=10, max_nbody=4)
        assert count == 385

    def test_estimate_hmbe_reduction_factor(self):
        """Test HMBE reduction factor estimation."""
        # For 100 fragments, 10 groups, (2,4)-HMBE
        # Should give substantial reduction
        reduction = estimate_hmbe_reduction_factor(
            num_fragments=100,
            truncation_orders=(2, 4),
            frags_per_group=10
        )

        # Should be > 1 (speedup)
        assert reduction > 1.0

        # For very aggressive HMBE, should be large
        assert reduction > 10.0, \
            f"Expected large reduction for 100 fragments, got {reduction:.1f}x"

    def test_get_enumeration_statistics(self):
        """Test statistics generation."""
        from qcmanybody.builder import build_nbody_compute_list

        hierarchy = create_2tier_4x4_hierarchy()
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy
        )

        # Generate both MBE and HMBE terms
        mbe_compute_list = build_nbody_compute_list(
            bsse_type=["nocp"],
            nfragments=16,
            nbodies=[1, 2, 3],
            return_total_data=False,
            supersystem_ie_only=False,
        )

        # Extract MBE terms
        mbe_terms = set()
        for nbody_dict in mbe_compute_list["all"].values():
            mbe_terms.update(nbody_dict)

        # Get HMBE terms
        hmbe_terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Get statistics
        stats = get_enumeration_statistics(hmbe_spec, mbe_terms, hmbe_terms)

        # Verify statistics
        assert stats["num_fragments"] == 16
        assert stats["truncation_orders"] == (2, 3)
        assert stats["num_tiers"] == 2
        assert stats["mbe_term_count"] == len(mbe_terms)
        assert stats["hmbe_term_count"] == len(hmbe_terms)
        assert stats["reduction_factor"] == len(mbe_terms) / len(hmbe_terms)
        assert stats["terms_avoided"] == len(mbe_terms) - len(hmbe_terms)


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_minimal_2tier_system(self):
        """Test smallest possible 2-tier system (2 fragments, 2 groups)."""
        fragment_tiers = {
            1: ("G0", "F0"),
            2: ("G1", "F1"),
        }
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers=fragment_tiers
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 2),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Should have: 2 monomers + 1 dimer = 3 terms
        assert len(terms) == 3
        assert (1,) in terms
        assert (2,) in terms
        assert (1, 2) in terms

    def test_uniform_truncation_equals_mbe(self):
        """Test that when T_1 >= num_groups, HMBE equals standard MBE."""
        hierarchy = create_2tier_4x4_hierarchy()

        # (4,4)-HMBE with 4 tier-1 groups should include all MBE-4 terms
        # because we can use all 4 groups (T_1=4)
        hmbe_spec = HMBESpecification(
            truncation_orders=(4, 4),
            hierarchy=hierarchy
        )

        hmbe_terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Should equal full MBE-4 term count: C(16,1)+C(16,2)+C(16,3)+C(16,4) = 2516
        assert len(hmbe_terms) == 2516

    def test_very_restrictive_truncation(self):
        """Test very restrictive (1,2)-HMBE."""
        hierarchy = create_2tier_4x4_hierarchy()

        # (1,2)-HMBE: only 1 tier-1 group allowed, max 2-body
        hmbe_spec = HMBESpecification(
            truncation_orders=(1, 2),
            hierarchy=hierarchy
        )

        terms = enumerate_hmbe_terms_2tier(hmbe_spec)

        # Should have:
        # - All 16 monomers (always included)
        # - Dimers within each group: 4 groups × C(4,2) = 4 × 6 = 24
        # Total: 16 + 24 = 40
        assert len(terms) == 40

        # Verify no dimers span multiple tier-1 groups
        dimer_terms = [t for t in terms if len(t) == 2]
        for term in dimer_terms:
            tier1_groups = {hierarchy.fragment_tiers[f][0] for f in term}
            assert len(tier1_groups) == 1, \
                "With T_1=1, no dimers should span multiple tier-1 groups"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
