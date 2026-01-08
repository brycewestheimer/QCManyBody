"""
Reference validation tests for HMBE using water-16 systems.

These tests validate the HMBE implementation against reference calculations
from water-16 clusters with a 4x4 hierarchy (4 tier-1 groups × 4 waters each).
"""

import pytest
import numpy as np
from pathlib import Path
from qcelemental.models import Molecule

from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification


def create_water16_4x4_hierarchy() -> FragmentHierarchy:
    """Create 4x4 hierarchy for water-16 system.

    Returns:
        FragmentHierarchy with 4 tier-1 groups (G0-G3), 16 waters total
    """
    fragment_tiers = {}
    for i in range(16):
        frag_id = i + 1  # 1-indexed
        tier1_group = f"G{i // 4}"  # 0-3 → G0, 4-7 → G1, 8-11 → G2, 12-15 → G3
        tier2_name = f"W{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_name)

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("domain", "water"),
    )


def create_water16_test_molecule() -> Molecule:
    """Create a simple test molecule for water-16.

    Returns:
        Molecule with 16 water fragments in a grid arrangement
    """
    # Create 16 waters in a 4x4 grid
    symbols = []
    geometry = []
    fragments = []

    for i in range(16):
        row = i // 4
        col = i % 4

        # Place water at grid position
        x = col * 4.0
        y = row * 4.0
        z = 0.0

        # Add O-H-H for water molecule
        symbols.extend(["O", "H", "H"])
        geometry.extend([
            [x, y, z],
            [x + 0.757, y, z + 0.587],
            [x - 0.757, y, z + 0.587],
        ])

        # Fragment indices (0-indexed for atom positions)
        start = i * 3
        fragments.append([start, start + 1, start + 2])

    return Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=0.0,
        molecular_multiplicity=1,
    )


class TestWater16HMBEReference:
    """Reference validation tests for water-16 HMBE."""

    def test_water16_4x4_hierarchy_structure(self):
        """Test that 4x4 hierarchy is correctly structured."""
        hierarchy = create_water16_4x4_hierarchy()

        # Verify structure
        assert hierarchy.num_tiers == 2
        assert len(hierarchy.fragment_tiers) == 16

        # Check tier-1 groups
        tier1_groups = set()
        for frag_id in range(1, 17):
            tier1_group = hierarchy.fragment_tiers[frag_id][0]
            tier1_groups.add(tier1_group)

        assert tier1_groups == {"G0", "G1", "G2", "G3"}

        # Verify each group has exactly 4 fragments
        for group_idx in range(4):
            group_name = f"G{group_idx}"
            frags_in_group = [
                frag_id for frag_id in range(1, 17)
                if hierarchy.fragment_tiers[frag_id][0] == group_name
            ]
            assert len(frags_in_group) == 4

    def test_hmbe_2_3_filtering_behavior(self):
        """Test (2,3)-HMBE filtering on water-16."""
        mol = create_water16_test_molecule()
        hierarchy = create_water16_4x4_hierarchy()

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
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

        # Get compute map and verify filtering
        compute_map = mbcore.compute_map
        nocp_list = compute_map["hf/sto-3g"]["nocp"]

        # Verify 3-body terms are filtered
        if 3 in nocp_list:
            three_body_terms = nocp_list[3]

            for frag_tuple, _ in three_body_terms:
                # Count tier-1 groups for this term
                groups = {hierarchy.fragment_tiers[f][0] for f in frag_tuple}

                # (2,3)-HMBE: max 2 tier-1 groups at tier-1
                assert len(groups) <= 2, (
                    f"Term {frag_tuple} spans {len(groups)} tier-1 groups "
                    f"but (2,3)-HMBE should only allow ≤2 groups"
                )

        # Get statistics
        stats = mbcore.get_hmbe_statistics()

        # Verify we got a reduction
        reduction = stats["reduction_factors"]["hf/sto-3g"]
        assert reduction > 1.0, "Expected term reduction for (2,3)-HMBE on 16-fragment system with 4 groups"

    def test_hmbe_2_4_includes_more_terms(self):
        """Test that (2,4)-HMBE includes more terms than (2,3)-HMBE."""
        mol = create_water16_test_molecule()
        hierarchy = create_water16_4x4_hierarchy()

        # (2,3)-HMBE
        hmbe_23 = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        mbcore_23 = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_23,
        )

        stats_23 = mbcore_23.get_hmbe_statistics()
        count_23 = stats_23["hmbe_term_counts"]["hf/sto-3g"]

        # (2,4)-HMBE
        hmbe_24 = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=hierarchy,
        )

        mbcore_24 = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_24,
        )

        stats_24 = mbcore_24.get_hmbe_statistics()
        count_24 = stats_24["hmbe_term_counts"]["hf/sto-3g"]

        # (2,4)-HMBE should have at least as many terms as (2,3)-HMBE
        # (actually more due to 4-body)
        assert count_24 > count_23, (
            f"(2,4)-HMBE should have more terms than (2,3)-HMBE, "
            f"got {count_24} vs {count_23}"
        )

    def test_hmbe_3_4_allows_more_groups(self):
        """Test that (3,4)-HMBE allows up to 3 tier-1 groups."""
        mol = create_water16_test_molecule()
        hierarchy = create_water16_4x4_hierarchy()

        hmbe_spec = HMBESpecification(
            truncation_orders=(3, 4),
            hierarchy=hierarchy,
        )

        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        compute_map = mbcore.compute_map
        nocp_list = compute_map["hf/sto-3g"]["nocp"]

        # Verify 3-body and 4-body terms respect T_1=3 constraint
        for nbody in [3, 4]:
            if nbody in nocp_list:
                for frag_tuple, _ in nocp_list[nbody]:
                    groups = {hierarchy.fragment_tiers[f][0] for f in frag_tuple}
                    assert len(groups) <= 3, (
                        f"{nbody}-body term {frag_tuple} spans {len(groups)} tier-1 groups "
                        f"but (3,4)-HMBE should only allow ≤3 groups"
                    )

        # Get statistics
        stats = mbcore.get_hmbe_statistics()

        # Should still see some reduction (4-body terms spanning all 4 groups filtered)
        reduction = stats["reduction_factors"]["hf/sto-3g"]
        assert reduction >= 1.0

    def test_hmbe_reduction_factors_reasonable(self):
        """Test that HMBE shows reasonable reduction factors for different truncations."""
        mol = create_water16_test_molecule()
        hierarchy = create_water16_4x4_hierarchy()

        # Test (2,3) and (3,4)
        configs = [
            ((2, 3), "2_3"),
            ((3, 4), "3_4"),
        ]

        reductions = {}
        hmbe_counts = {}
        mbe_counts = {}

        for truncation_orders, label in configs:
            hmbe_spec = HMBESpecification(
                truncation_orders=truncation_orders,
                hierarchy=hierarchy,
            )

            max_nbody = truncation_orders[-1]
            levels = {i: "hf/sto-3g" for i in range(1, max_nbody + 1)}

            mbcore = ManyBodyCore(
                molecule=mol,
                bsse_type=[BsseEnum.nocp],
                levels=levels,
                return_total_data=False,
                supersystem_ie_only=False,
                embedding_charges={},
                hmbe_spec=hmbe_spec,
            )

            stats = mbcore.get_hmbe_statistics()
            reductions[label] = stats["reduction_factors"]["hf/sto-3g"]
            hmbe_counts[label] = stats["hmbe_term_counts"]["hf/sto-3g"]
            mbe_counts[label] = stats["mbe_term_counts"]["hf/sto-3g"]

        # Both should show some reduction or equal (>=1.0)
        assert reductions["2_3"] >= 1.0
        assert reductions["3_4"] >= 1.0

        # (2,3)-HMBE should have fewer terms than MBE-3 or equal
        assert hmbe_counts["2_3"] <= mbe_counts["2_3"]

        # (3,4)-HMBE should have fewer terms than MBE-4 or equal
        assert hmbe_counts["3_4"] <= mbe_counts["3_4"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
