"""
End-to-end HMBE workflow tests.

Tests the complete workflow from input specification through to result generation.
These tests validate the full integration without requiring actual QC calculations.
"""

import pytest
import json
from pathlib import Path
from qcelemental.models import Molecule

from qcmanybody.core import ManyBodyCore
from qcmanybody import ManyBodyComputer
from qcmanybody.models.v1 import BsseEnum, ManyBodyInput, ManyBodySpecification, ManyBodyKeywords
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification, SchengenSpecification
from qcmanybody.cli.input_parser import parse_input_file
from qcmanybody.cli.converter import convert_to_manybody_input


class TestHMBEEndToEnd:
    """End-to-end workflow tests for HMBE."""

    def test_python_api_workflow(self):
        """Test complete workflow using Python API."""
        # Create a realistic 9-water system with 3x3 hierarchy
        water_positions = [
            # Group 0 (3 waters)
            [[0, 0, 0], [0, 0.757, 0.587], [0, -0.757, 0.587]],
            [[3, 0, 0], [3, 0.757, 0.587], [3, -0.757, 0.587]],
            [[0, 3, 0], [0, 3.757, 0.587], [0, 2.243, 0.587]],
            # Group 1 (3 waters)
            [[6, 0, 0], [6, 0.757, 0.587], [6, -0.757, 0.587]],
            [[9, 0, 0], [9, 0.757, 0.587], [9, -0.757, 0.587]],
            [[6, 3, 0], [6, 3.757, 0.587], [6, 2.243, 0.587]],
            # Group 2 (3 waters)
            [[12, 0, 0], [12, 0.757, 0.587], [12, -0.757, 0.587]],
            [[15, 0, 0], [15, 0.757, 0.587], [15, -0.757, 0.587]],
            [[12, 3, 0], [12, 3.757, 0.587], [12, 2.243, 0.587]],
        ]

        geometry = []
        fragments = []
        for i, water in enumerate(water_positions):
            start_idx = i * 3
            fragments.append([start_idx, start_idx + 1, start_idx + 2])
            geometry.extend(water)

        mol = Molecule(
            symbols=["O", "H", "H"] * 9,
            geometry=geometry,
            fragments=fragments,
        )

        # Define (2,3)-HMBE with 3 tier-1 groups
        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "W0"), 2: ("G0", "W1"), 3: ("G0", "W2"),
                4: ("G1", "W3"), 5: ("G1", "W4"), 6: ("G1", "W5"),
                7: ("G2", "W6"), 8: ("G2", "W7"), 9: ("G2", "W8"),
            },
            tier_names=("domain", "water"),
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.1,
                distance_metric="R2",
            ),
        )

        # Create ManyBodyCore with HMBE
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "mp2/sto-3g", 3: "mp2/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        # Get HMBE statistics
        stats = mbcore.get_hmbe_statistics()

        # Verify statistics are generated
        assert stats is not None
        assert "mbe_term_counts" in stats
        assert "hmbe_term_counts" in stats
        assert "reduction_factors" in stats

        # With 9 fragments and (2,3)-HMBE spanning 3 groups:
        # Standard MBE-3: C(9,1) + C(9,2) + C(9,3) = 9 + 36 + 84 = 129 fragment tuples
        # HMBE (2,3): Only terms from ≤2 tier-1 groups allowed
        # 3-body terms spanning all 3 groups should be filtered

        # Verify we got a reduction
        mp2_reduction = stats["reduction_factors"]["mp2/sto-3g"]
        assert mp2_reduction > 1.0, "HMBE should reduce term count for 9-fragment system with 3 groups"

        # Verify the compute map is generated correctly
        compute_map = mbcore.compute_map
        assert "mp2/sto-3g" in compute_map
        assert "cp" in compute_map["mp2/sto-3g"]

        # Count actual terms
        cp_terms = compute_map["mp2/sto-3g"]["cp"]
        total_terms = sum(len(terms) for terms in cp_terms.values())
        assert total_terms > 0

    def test_cli_workflow(self):
        """Test complete workflow using CLI input file."""
        example_file = Path(__file__).parent.parent.parent / "examples" / "cli" / "hmbe_water6_example.json"

        if not example_file.exists():
            pytest.skip(f"Example file not found: {example_file}")

        # Parse input file
        cli_input = parse_input_file(str(example_file))

        # Convert to ManyBodyInput
        mb_input = convert_to_manybody_input(cli_input, str(example_file))

        # Create ManyBodyComputer from ManyBodyInput (without building tasks)
        mbc = ManyBodyComputer.from_manybodyinput(mb_input, build_tasks=False)

        # Verify HMBE is configured
        assert mbc.qcmb_core.hmbe_spec is not None
        assert mbc.qcmb_core.hmbe_spec.truncation_orders == (2, 3)

        # Get statistics
        stats = mbc.qcmb_core.get_hmbe_statistics()
        assert stats is not None

        # Verify term filtering occurred
        assert "reduction_factors" in stats
        # With 6 fragments in 3 groups and (2,3)-HMBE, we should see some reduction
        # For 3-body terms spanning all 3 groups

        # Get the actual key (single-level uses max_nbody as key)
        reduction_factors = stats["reduction_factors"]
        assert len(reduction_factors) > 0, "Should have at least one reduction factor"

        # Check that we got some reduction (at least one level should show reduction)
        any_reduction = any(rf > 1.0 for rf in reduction_factors.values())
        assert any_reduction, f"HMBE should show reduction for some level, got: {reduction_factors}"

    def test_hmbe_workflow_with_term_verification(self):
        """Test HMBE workflow with detailed term verification."""
        # Create a 6-fragment system with 3 tier-1 groups (2 fragments each)
        mol = Molecule(
            symbols=["He"] * 6,
            geometry=[[float(i * 3), 0, 0] for i in range(6)],
            fragments=[[i] for i in range(6)],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"), 2: ("G0", "S1"),
                3: ("G1", "S2"), 4: ("G1", "S3"),
                5: ("G2", "S4"), 6: ("G2", "S5"),
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

        # Get compute map and verify specific terms
        compute_map = mbcore.compute_map
        nocp_list = compute_map["hf/sto-3g"]["nocp"]

        # Check 3-body terms
        if 3 in nocp_list:
            three_body_terms = nocp_list[3]

            for frag_tuple, _ in three_body_terms:
                # Count tier-1 groups
                groups = {hierarchy.fragment_tiers[f][0] for f in frag_tuple}

                # With (2,3)-HMBE, no 3-body term should span all 3 tier-1 groups
                assert len(groups) <= 2, (
                    f"3-body term {frag_tuple} spans {len(groups)} tier-1 groups "
                    f"({groups}), but (2,3)-HMBE should filter terms with >2 groups"
                )

                # Specifically, (1,3,5) should be FILTERED (spans G0, G1, G2)
                # And (1,2,3) should be KEPT (spans G0, G1 only)
                if frag_tuple == (1, 3, 5):
                    pytest.fail(f"Term {frag_tuple} should have been filtered by HMBE (spans 3 groups)")

            # Verify (1,2,3) IS present (spans only 2 groups: G0, G1)
            frag_tuples = {frag for frag, _ in three_body_terms}
            assert (1, 2, 3) in frag_tuples, "Term (1,2,3) should be present (spans only 2 tier-1 groups)"

    def test_hmbe_schengen_workflow(self):
        """Test HMBE workflow with Schengen term selection."""
        mol = Molecule(
            symbols=["He"] * 6,
            geometry=[
                [0, 0, 0], [3, 0, 0],  # Group 0
                [0, 3, 0], [3, 3, 0],  # Group 1
                [0, 6, 0], [3, 6, 0],  # Group 2
            ],
            fragments=[[i] for i in range(6)],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"), 2: ("G0", "S1"),
                3: ("G1", "S2"), 4: ("G1", "S3"),
                5: ("G2", "S4"), 6: ("G2", "S5"),
            },
        )

        # Base (2,3)-HMBE
        hmbe_base = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        mbcore_base = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_base,
        )

        base_stats = mbcore_base.get_hmbe_statistics()
        base_count = base_stats["hmbe_term_counts"]["hf/sto-3g"]

        # (2,3)-HMBE with Schengen
        hmbe_schengen = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=SchengenSpecification(
                enabled=True,
                selection_fraction=0.2,  # Add back 20% of filtered terms
                distance_metric="R2",
            ),
        )

        mbcore_schengen = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_schengen,
        )

        schengen_stats = mbcore_schengen.get_hmbe_statistics()
        schengen_count = schengen_stats["hmbe_term_counts"]["hf/sto-3g"]

        # Schengen should add back some terms, so count should be higher
        assert schengen_count >= base_count, (
            "HMBE with Schengen should have at least as many terms as base HMBE "
            f"(got {schengen_count} vs {base_count})"
        )

        # Verify Schengen is reported in stats
        assert schengen_stats.get("schengen_enabled") is True

    def test_hmbe_33_standard_mbe_equivalence(self):
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

        # Reduction factor should be 1.0 (no reduction)
        assert stats["reduction_factors"]["hf/sto-3g"] == 1.0
        assert stats["hmbe_term_counts"]["hf/sto-3g"] == stats["mbe_term_counts"]["hf/sto-3g"]


class TestHMBEOutputMetadata:
    """Test HMBE metadata in output results."""

    def test_hmbe_metadata_structure(self):
        """Test that HMBE metadata is properly structured in output."""
        mol = Molecule(
            symbols=["He"] * 4,
            geometry=[[float(i), 0, 0] for i in range(4)],
            fragments=[[i] for i in range(4)],
        )

        hierarchy = FragmentHierarchy(
            num_tiers=2,
            fragment_tiers={
                1: ("G0", "S0"), 2: ("G0", "S1"),
                3: ("G1", "S2"), 4: ("G1", "S3"),
            },
        )

        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
        )

        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "mp2/sto-3g", 3: "mp2/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec,
        )

        stats = mbcore.get_hmbe_statistics()

        # Verify required fields
        assert "mbe_term_counts" in stats
        assert "hmbe_term_counts" in stats
        assert "reduction_factors" in stats
        assert "truncation_orders" in stats
        assert "num_tiers" in stats

        # Verify values
        assert stats["truncation_orders"] == (2, 3)
        assert stats["num_tiers"] == 2

        # Verify per-level data
        assert "hf/sto-3g" in stats["mbe_term_counts"]
        assert "mp2/sto-3g" in stats["mbe_term_counts"]
        assert "hf/sto-3g" in stats["reduction_factors"]
        assert "mp2/sto-3g" in stats["reduction_factors"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
