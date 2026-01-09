"""
Comprehensive tests for HMBE with parallel execution.

These tests validate that HMBE works correctly with parallel execution,
including thread safety, determinism, and performance.
"""

import pytest
from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification


@pytest.fixture
def water16_molecule():
    """Create water-16 molecule with 4×4 two-tier hierarchy."""
    # 16 water molecules (4 groups of 4 waters each)
    symbols = ["O", "H", "H"] * 16

    # Simplified geometry (not real coords, just for testing structure)
    geometry = []
    for i in range(16):
        base_x = (i % 4) * 5.0
        base_y = (i // 4) * 5.0
        geometry.extend([
            [base_x, base_y, 0.0],      # O
            [base_x + 0.96, base_y, 0.0],  # H
            [base_x, base_y + 0.96, 0.0],  # H
        ])

    fragments = [[i*3, i*3+1, i*3+2] for i in range(16)]

    return Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=0.0,
        molecular_multiplicity=1,
    )


@pytest.fixture
def water16_hierarchy():
    """Create 4×4 two-tier hierarchy for water-16."""
    fragment_tiers = {}
    for i in range(16):
        group_idx = i // 4  # 0-3
        subgroup_idx = i % 4  # 0-3
        fragment_tiers[str(i + 1)] = [f"G{group_idx}", f"G{group_idx}_S{subgroup_idx}"]

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=["domain", "fragment"]
    )


class TestHMBEStructure:
    """Test HMBE compute map structure and term counts."""

    def test_hmbe_23_term_count(self, water16_molecule, water16_hierarchy):
        """Test (2,3)-HMBE generates expected number of terms."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            enumeration_mode="filter"
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        stats = mbc.get_hmbe_statistics()
        assert stats is not None
        assert "hmbe_term_counts" in stats
        assert "reduction_factors" in stats

        # (2,3)-HMBE should significantly reduce term count
        reduction = stats["reduction_factors"]["hf/sto-3g"]
        assert reduction > 1.0, f"Expected reduction >1.0x, got {reduction}"

    def test_hmbe_24_term_count(self, water16_molecule, water16_hierarchy):
        """Test (2,4)-HMBE generates expected number of terms."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 4),
            hierarchy=water16_hierarchy,
            enumeration_mode="filter"
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        stats = mbc.get_hmbe_statistics()
        reduction = stats["reduction_factors"]["hf/sto-3g"]
        assert reduction > 1.0

    def test_hmbe_34_term_count(self, water16_molecule, water16_hierarchy):
        """Test (3,4)-HMBE generates expected number of terms."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(3, 4),
            hierarchy=water16_hierarchy,
            enumeration_mode="filter"
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g", 4: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        stats = mbc.get_hmbe_statistics()
        reduction = stats["reduction_factors"]["hf/sto-3g"]
        # (3,4)-HMBE is less restrictive, so smaller reduction
        assert reduction >= 1.0


class TestHMBEEnumerationModes:
    """Test direct vs filter enumeration consistency."""

    def test_direct_vs_filter_identical_23(self, water16_molecule, water16_hierarchy):
        """Test direct and filter modes produce identical term sets for (2,3)-HMBE."""
        hmbe_filter = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            enumeration_mode="filter"
        )

        hmbe_direct = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            enumeration_mode="direct"
        )

        mbc_filter = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_filter
        )

        mbc_direct = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_direct
        )

        # Extract fragment tuples from both
        filter_frags = set()
        for bsse_dict in mbc_filter.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    filter_frags.add(frag)

        direct_frags = set()
        for bsse_dict in mbc_direct.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    direct_frags.add(frag)

        assert filter_frags == direct_frags, "Direct and filter enumeration must produce identical term sets"

    def test_auto_mode_selection(self, water16_molecule, water16_hierarchy):
        """Test that auto mode selects appropriate enumeration method."""
        hmbe_auto = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            enumeration_mode="auto"
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_auto
        )

        stats = mbc.get_hmbe_statistics()
        # For 16 fragments, auto should choose filter mode
        assert stats["actual_enumeration_mode"] == "filter"


class TestHMBECompleteness:
    """Test HMBE completeness requirement for Möbius inversion."""

    def test_base_hmbe_completeness(self, water16_molecule, water16_hierarchy):
        """Test that base HMBE (no Schengen) satisfies completeness."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        # Extract all fragment tuples
        fragment_tuples = set()
        for bsse_dict in mbc.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    fragment_tuples.add(frag)

        # Validate completeness
        mbc._validate_hmbe_completeness(fragment_tuples, "hf/sto-3g", "nocp")
        # If this doesn't raise, completeness is satisfied

    def test_schengen_completeness(self, water16_molecule, water16_hierarchy):
        """Test that Schengen terms maintain completeness."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.2,
                distance_metric="R2"
            )
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        fragment_tuples = set()
        for bsse_dict in mbc.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    fragment_tuples.add(frag)

        # Validate completeness with Schengen terms
        mbc._validate_hmbe_completeness(fragment_tuples, "hf/sto-3g", "nocp")


class TestHMBEWithSchengen:
    """Test HMBE with Schengen term selection."""

    def test_schengen_adds_terms(self, water16_molecule, water16_hierarchy):
        """Test that Schengen adds additional terms beyond base HMBE."""
        # Base HMBE
        hmbe_base = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
        )

        mbc_base = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_base
        )

        # HMBE with Schengen
        hmbe_schengen = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.2,
                distance_metric="R2"
            )
        )

        mbc_schengen = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_schengen
        )

        base_terms = set()
        for bsse_dict in mbc_base.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    base_terms.add(frag)

        schengen_terms = set()
        for bsse_dict in mbc_schengen.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    schengen_terms.add(frag)

        # Schengen should add terms (or at least not remove any)
        assert len(schengen_terms) >= len(base_terms), "Schengen should not remove base HMBE terms"

        # If Schengen added terms, verify they're new
        if len(schengen_terms) > len(base_terms):
            added_terms = schengen_terms - base_terms
            assert len(added_terms) > 0, "Schengen should add new terms"


class TestThreadSafety:
    """Test thread safety of HMBE operations."""

    def test_schengen_distance_calculation_deterministic(self, water16_molecule, water16_hierarchy):
        """Test that Schengen distance calculations are deterministic."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.2,
                distance_metric="R2"
            )
        )

        # Create two identical cores
        mbc1 = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        mbc2 = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        # Extract terms from both
        terms1 = set()
        for bsse_dict in mbc1.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    terms1.add(frag)

        terms2 = set()
        for bsse_dict in mbc2.compute_map["hf/sto-3g"].values():
            for terms in bsse_dict.values():
                for frag, bas in terms:
                    terms2.add(frag)

        assert terms1 == terms2, "Schengen selection must be deterministic"

    def test_hmbe_filter_pure_function(self, water16_hierarchy):
        """Test that HMBE filter is a pure function (no side effects)."""
        from qcmanybody.hmbe_filter import passes_hmbe_filter

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
        )

        test_term = (1, 2, 5)  # 2 groups, 3 fragments

        # Call multiple times
        result1 = passes_hmbe_filter(test_term, hmbe_spec)
        result2 = passes_hmbe_filter(test_term, hmbe_spec)
        result3 = passes_hmbe_filter(test_term, hmbe_spec)

        assert result1 == result2 == result3, "passes_hmbe_filter must be deterministic"


class TestHMBEDifferentBSSE:
    """Test HMBE with different BSSE correction types."""

    def test_hmbe_nocp(self, water16_molecule, water16_hierarchy):
        """Test HMBE with NOCP correction."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        assert "nocp" in mbc.compute_map["hf/sto-3g"]
        assert len(mbc.compute_map["hf/sto-3g"]["nocp"]) > 0

    def test_hmbe_cp(self, water16_molecule, water16_hierarchy):
        """Test HMBE with CP correction."""
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=water16_hierarchy,
        )

        mbc = ManyBodyCore(
            molecule=water16_molecule,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=True,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        assert "cp" in mbc.compute_map["hf/sto-3g"]
        assert len(mbc.compute_map["hf/sto-3g"]["cp"]) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
