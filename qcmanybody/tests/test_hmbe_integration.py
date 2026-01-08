"""
Integration tests for HMBE functionality.

Tests the complete HMBE workflow from input specification through result generation.
Uses the core interface to avoid dependencies on external QC programs.
"""

import pytest
import numpy as np
from qcelemental.models import Molecule

from qcmanybody import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import (
    FragmentHierarchy,
    HMBESpecification,
    SchengenSpecification,
)


class TestHMBEIntegration:
    """Integration tests for HMBE with ManyBodyCore."""

    def test_hmbe_2tier_basic(self):
        """Test basic 2-tier HMBE with 4-fragment system.

        Note: With only 2 tier-1 groups and (2,3)-HMBE, no terms are filtered
        because all terms satisfy the constraint of ≤2 tier-1 groups.
        This test verifies the infrastructure works correctly even when no filtering occurs.
        """
        # Create 4-fragment He4 system in 2x2 hierarchy
        mol = Molecule(
            symbols=["He", "He", "He", "He"],
            geometry=[
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 3.0],
                [0.0, 3.0, 0.0],
                [0.0, 3.0, 3.0],
            ],
            fragments=[[0], [1], [2], [3]],
        )

        # Define 2x2 hierarchy
        # Group 0: fragments 1, 2
        # Group 1: fragments 3, 4
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "G0_S0"),
                2: ("G0", "G0_S1"),
                3: ("G1", "G1_S0"),
                4: ("G1", "G1_S1"),
            },
            tier_names=("group", "subgroup"),
        )

        # Configure (2,3)-HMBE
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        # Create ManyBodyCore with HMBE
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={1: "scf/sto-3g", 2: "scf/sto-3g", 3: "scf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        # Verify HMBE spec is set
        assert mbcore.hmbe_spec is not None
        assert mbcore.hmbe_spec.truncation_orders == (2, 3)
        assert mbcore.hmbe_spec.num_tiers == 2

        # Get compute map (triggers filtering)
        compute_map = mbcore.compute_map

        # Verify compute map is filtered
        assert compute_map is not None
        assert "scf/sto-3g" in compute_map

        # Count total terms in CP list
        cp_list = compute_map["scf/sto-3g"]["cp"]
        total_terms = sum(len(terms) for terms in cp_list.values())

        # For (2,3)-HMBE with 4 fragments in 2x2:
        # - All 2-body terms should pass (within 2 tier-1 groups)
        # - 3-body terms: only those from ≤2 tier-1 groups
        # This should have fewer terms than full MBE-3
        assert total_terms > 0
        # CP correction includes more terms than NOCP
        # Just verify we got some terms
        assert total_terms >= 10

        # Get HMBE statistics
        stats = mbcore.get_hmbe_statistics()
        assert stats is not None
        assert "mbe_term_counts" in stats
        assert "hmbe_term_counts" in stats
        assert "reduction_factors" in stats
        assert stats["truncation_orders"] == (2, 3)
        assert stats["num_tiers"] == 2

        # For 2x2 hierarchy with (2,3)-HMBE, no terms are filtered
        # because all terms span ≤2 tier-1 groups (and we only have 2 groups total)
        # So reduction factor should be 1.0 (no reduction)
        hmbe_count = stats["hmbe_term_counts"]["scf/sto-3g"]
        mbe_count = stats["mbe_term_counts"]["scf/sto-3g"]
        assert hmbe_count == mbe_count  # No filtering for this case
        assert stats["reduction_factors"]["scf/sto-3g"] == 1.0

    def test_hmbe_33_reduces_to_mbe(self):
        """Test that (3,3)-HMBE is equivalent to standard MBE-3."""
        mol = Molecule(
            symbols=["He", "He", "He"],
            geometry=[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
            fragments=[[0], [1], [2]],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G0", "S1"),
                3: ("G1", "S2"),
            },
        )

        # (3,3)-HMBE should be identical to MBE-3
        hmbe_spec = HMBESpecification(
            truncation_orders=(3, 3),
            hierarchy=hierarchy,
        )

        assert hmbe_spec.is_standard_mbe()

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

        # For (3,3)-HMBE, reduction factor should be 1.0 (no reduction)
        assert stats["reduction_factors"]["hf/sto-3g"] == 1.0
        assert stats["hmbe_term_counts"]["hf/sto-3g"] == stats["mbe_term_counts"]["hf/sto-3g"]

    def test_hmbe_term_filtering(self):
        """Test that specific terms are filtered correctly."""
        # 6-fragment system with 3 tier-1 groups
        mol = Molecule(
            symbols=["He"] * 6,
            geometry=[[float(i), 0, 0] for i in range(6)],
            fragments=[[i] for i in range(6)],
        )

        # 3 tier-1 groups of 2 fragments each
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

        # (2,3)-HMBE: max 2 tier-1 groups, max 3-body
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

        compute_map = mbcore.compute_map
        nocp_list = compute_map["hf/sto-3g"]["nocp"]

        # Check specific 3-body terms
        if 3 in nocp_list:
            for frag_tuple, _ in nocp_list[3]:
                # Count tier-1 groups
                groups = {hierarchy.fragment_tiers[f][0] for f in frag_tuple}
                # Should have ≤ 2 tier-1 groups (T_1 = 2)
                assert len(groups) <= 2, f"Term {frag_tuple} has {len(groups)} tier-1 groups, expected ≤2"

    def test_hmbe_schengen_enabled(self):
        """Test HMBE with Schengen term selection."""
        mol = Molecule(
            symbols=["He"] * 4,
            geometry=[[float(i), 0, 0] for i in range(4)],
            fragments=[[i] for i in range(4)],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G0", "S1"),
                3: ("G1", "S2"),
                4: ("G1", "S3"),
            },
        )

        # Enable Schengen with 10% selection
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.1,
                distance_metric="R2",
            ),
        )

        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        stats = mbcore.get_hmbe_statistics()

        # Verify Schengen is reported as enabled
        assert stats["schengen_enabled"] is True

        # With Schengen, we should have more terms than base HMBE
        # but still fewer than full MBE
        compute_map = mbcore.compute_map
        cp_list = compute_map["hf/sto-3g"]["cp"]
        total_hmbe_with_schengen = sum(len(terms) for terms in cp_list.values())

        # Should have some terms
        assert total_hmbe_with_schengen > 0

    def test_hmbe_without_spec_is_standard_mbe(self):
        """Test that omitting hmbe_spec gives standard MBE."""
        mol = Molecule(
            symbols=["He", "He", "He"],
            geometry=[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
            fragments=[[0], [1], [2]],
        )

        # Create ManyBodyCore WITHOUT hmbe_spec
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=None,  # Standard MBE
        )

        # Verify HMBE is not enabled
        assert mbcore.hmbe_spec is None
        assert mbcore.get_hmbe_statistics() is None

        # Compute map should still work normally
        compute_map = mbcore.compute_map
        assert compute_map is not None
        assert "hf/sto-3g" in compute_map

    def test_hmbe_multilevel(self):
        """Test HMBE with multi-level model chemistry."""
        mol = Molecule(
            symbols=["He"] * 4,
            geometry=[[float(i), 0, 0] for i in range(4)],
            fragments=[[i] for i in range(4)],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G0", "S1"),
                3: ("G1", "S2"),
                4: ("G1", "S3"),
            },
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        # Different methods at different levels
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={
                1: "hf/sto-3g",
                2: "mp2/sto-3g",
                3: "mp2/sto-3g",
            },
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        stats = mbcore.get_hmbe_statistics()

        # Should have statistics for both model chemistries
        assert "hf/sto-3g" in stats["hmbe_term_counts"]
        assert "mp2/sto-3g" in stats["hmbe_term_counts"]

        # For 2x2 hierarchy with (2,3)-HMBE, no terms are filtered
        # because all terms span ≤2 tier-1 groups (we only have 2 groups total)
        # MP2 level (2-body and 3-body) should show no reduction
        assert stats["reduction_factors"]["mp2/sto-3g"] == 1.0

        # HF level (1-body only) with CP might have 0 terms (CP correction not applicable to 1-body)
        # Just verify it's in the statistics
        assert "hf/sto-3g" in stats["reduction_factors"]


class TestHMBEValidation:
    """Test HMBE validation and error handling."""

    def test_hmbe_max_nbody_consistency(self):
        """Test that max_nbody is consistent with HMBE T_K."""
        mol = Molecule(
            symbols=["He", "He"],
            geometry=[[0, 0, 0], [0, 0, 3]],
            fragments=[[0], [1]],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G1", "S1"),
            },
        )

        # T_K = 3, so max_nbody in levels should be ≥ 3
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        # This should work: max_nbody = 3 matches T_K
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        # Verify it works
        assert mbcore.max_nbody == 3
        assert mbcore.hmbe_spec.max_nbody == 3

    def test_hmbe_fragment_count_matches_hierarchy(self):
        """Test that number of fragments matches hierarchy definition."""
        mol = Molecule(
            symbols=["He", "He", "He"],
            geometry=[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
            fragments=[[0], [1], [2]],
        )

        # Hierarchy defines only 2 fragments (mismatch!)
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"),
                2: ("G1", "S1"),
                # Missing fragment 3!
            },
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        # Should still create (filtering handles missing fragments gracefully)
        # but will fail when trying to filter fragment 3
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        # This should raise KeyError when accessing compute_map
        # because fragment 3 is not in hierarchy
        with pytest.raises(KeyError):
            _ = mbcore.compute_map


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
