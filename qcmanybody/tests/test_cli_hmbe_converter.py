"""
Tests for CLI HMBE converter functionality.
"""

import pytest
from qcmanybody.cli.input_parser import parse_input_file
from qcmanybody.cli.converter import convert_to_manybody_input, convert_hmbe_specification
from qcmanybody.cli.schemas.input_schema import HMBESchema, HierarchySchema, SchengenSchema
from qcmanybody.models.hierarchy import HMBESpecification


class TestHMBEConverter:
    """Test HMBE conversion from CLI schema to internal models."""

    def test_convert_hmbe_specification_basic(self):
        """Test basic HMBE specification conversion."""
        # Create CLI HMBE schema
        hmbe_schema = HMBESchema(
            truncation_orders=[2, 3],
            hierarchy=HierarchySchema(
                num_tiers=2,
                fragment_tiers={
                    1: ["G0", "S0"],
                    2: ["G0", "S1"],
                    3: ["G1", "S2"],
                    4: ["G1", "S3"],
                },
                tier_names=["domain", "fragment"],
            ),
        )

        # Convert to internal model
        hmbe_spec = convert_hmbe_specification(hmbe_schema)

        # Verify conversion
        assert isinstance(hmbe_spec, HMBESpecification)
        assert hmbe_spec.truncation_orders == (2, 3)
        assert hmbe_spec.num_tiers == 2
        assert hmbe_spec.hierarchy.num_tiers == 2
        assert hmbe_spec.hierarchy.tier_names == ("domain", "fragment")
        assert hmbe_spec.schengen is None

        # Verify fragment tiers are tuples
        assert hmbe_spec.hierarchy.fragment_tiers[1] == ("G0", "S0")
        assert hmbe_spec.hierarchy.fragment_tiers[2] == ("G0", "S1")

    def test_convert_hmbe_specification_with_schengen(self):
        """Test HMBE specification conversion with Schengen terms."""
        hmbe_schema = HMBESchema(
            truncation_orders=[2, 3],
            hierarchy=HierarchySchema(
                num_tiers=2,
                fragment_tiers={
                    1: ["G0", "S0"],
                    2: ["G0", "S1"],
                    3: ["G1", "S2"],
                },
            ),
            schengen=SchengenSchema(
                enabled=True,
                selection_fraction=0.15,
                distance_metric="R_inv",
            ),
        )

        hmbe_spec = convert_hmbe_specification(hmbe_schema)

        # Verify Schengen conversion
        assert hmbe_spec.schengen is not None
        assert hmbe_spec.schengen.enabled is True
        assert hmbe_spec.schengen.selection_fraction == 0.15
        assert hmbe_spec.schengen.distance_metric == "R_inv"

    def test_convert_hmbe_specification_3tier(self):
        """Test 3-tier HMBE specification conversion."""
        hmbe_schema = HMBESchema(
            truncation_orders=[2, 3, 4],
            hierarchy=HierarchySchema(
                num_tiers=3,
                fragment_tiers={
                    1: ["D0", "R0", "A0"],
                    2: ["D0", "R0", "A1"],
                    3: ["D1", "R1", "A2"],
                },
                tier_names=["domain", "residue", "atom"],
            ),
        )

        hmbe_spec = convert_hmbe_specification(hmbe_schema)

        assert hmbe_spec.truncation_orders == (2, 3, 4)
        assert hmbe_spec.num_tiers == 3
        assert hmbe_spec.hierarchy.tier_names == ("domain", "residue", "atom")
        assert hmbe_spec.hierarchy.fragment_tiers[1] == ("D0", "R0", "A0")

    def test_cli_json_to_manybody_input_with_hmbe(self):
        """Test full conversion from CLI JSON to ManyBodyInput with HMBE."""
        # Use the example JSON file
        import json
        from pathlib import Path

        example_file = Path(__file__).parent.parent.parent / "examples" / "cli" / "hmbe_water6_example.json"

        if not example_file.exists():
            pytest.skip(f"Example file not found: {example_file}")

        # Parse the input file
        cli_input = parse_input_file(str(example_file))

        # Verify CLI input has HMBE
        assert cli_input.manybody is not None
        assert cli_input.manybody.hmbe is not None
        assert cli_input.manybody.hmbe.truncation_orders == [2, 3]

        # Convert to ManyBodyInput
        mb_input = convert_to_manybody_input(cli_input, str(example_file))

        # Verify conversion
        assert mb_input.specification.keywords.hmbe_spec is not None
        hmbe_spec = mb_input.specification.keywords.hmbe_spec

        assert hmbe_spec.truncation_orders == (2, 3)
        assert hmbe_spec.num_tiers == 2
        assert hmbe_spec.schengen is not None
        assert hmbe_spec.schengen.enabled is True
        assert hmbe_spec.schengen.selection_fraction == 0.1
        assert hmbe_spec.schengen.distance_metric == "R2"

        # Verify hierarchy
        assert hmbe_spec.hierarchy.num_tiers == 2
        assert hmbe_spec.hierarchy.fragment_tiers[1] == ("G0", "G0_S0")
        assert hmbe_spec.hierarchy.fragment_tiers[6] == ("G2", "G2_S1")

        # Verify molecule
        assert len(mb_input.molecule.symbols) == 18  # 6 waters * 3 atoms
        assert len(mb_input.molecule.fragments) == 6

    def test_cli_without_hmbe(self):
        """Test that CLI conversion works without HMBE (standard MBE)."""
        from qcelemental.models import Molecule
        from qcmanybody.cli.schemas.input_schema import (
            QCManyBodyInput,
            MoleculeSchema,
            InlineMoleculeSchema,
            CalculationSchema,
            SingleLevelCalculationSchema,
            BsseSchema,
            ManyBodySchema,
        )

        # Create simple CLI input without HMBE
        cli_input = QCManyBodyInput(
            molecule=MoleculeSchema(
                source="inline",
                inline=InlineMoleculeSchema(
                    symbols=["He", "He"],
                    geometry=[[0, 0, 0], [0, 0, 3]],
                    fragments=[[0], [1]],
                ),
            ),
            calculation=CalculationSchema(
                type="single",
                single=SingleLevelCalculationSchema(
                    driver="energy",
                    method="scf",
                    basis="sto-3g",
                    program="psi4",
                ),
            ),
            bsse=BsseSchema(type=["cp"]),
            manybody=ManyBodySchema(max_nbody=2),
        )

        # Convert to ManyBodyInput
        mb_input = convert_to_manybody_input(cli_input)

        # Verify no HMBE
        assert mb_input.specification.keywords.hmbe_spec is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
