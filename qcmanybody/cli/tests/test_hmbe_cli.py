"""
Tests for HMBE CLI integration.

Verifies that HMBE hierarchical molecules can be specified via CLI input
and correctly converted to ManyBodyInput with HMBE mode enabled.
"""

import sys
from pathlib import Path

# Ensure we're using the local qcmanybody package
project_root = Path(__file__).parent.parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
from qcmanybody.cli.converter import convert_to_manybody_input
from qcmanybody.cli.molecule_loader import convert_hmbe_hierarchical
from qcmanybody.cli.schemas.input_schema import (
    BsseTypeEnum,
    BsseSchema,
    CalculationSchema,
    DriverEnum,
    HMBEFragmentSchema,
    HMBEHierarchicalMoleculeSchema,
    ManyBodySchema,
    MoleculeSchema,
    MoleculeSourceEnum,
    MultiLevelCalculationSchema,
    QCManyBodyInput,
    SingleLevelCalculationSchema,
)


class TestHMBESchemaValidation:
    """Test HMBE schema validation."""

    def test_simple_2tier_hierarchy_schema(self):
        """Test valid 2-tier hierarchy schema."""
        hmbe_schema = HMBEHierarchicalMoleculeSchema(
            tiers=2,
            max_primary_per_nmer=2,
            fragments=[
                HMBEFragmentSchema(
                    id="cluster_A",
                    sub_fragments=[
                        HMBEFragmentSchema(id="He_A1", symbols=["He"], geometry=[[0.0, 0.0, 0.0]]),
                        HMBEFragmentSchema(id="He_A2", symbols=["He"], geometry=[[1.0, 0.0, 0.0]]),
                    ],
                ),
                HMBEFragmentSchema(
                    id="cluster_B",
                    sub_fragments=[
                        HMBEFragmentSchema(id="He_B1", symbols=["He"], geometry=[[0.0, 1.0, 0.0]]),
                        HMBEFragmentSchema(id="He_B2", symbols=["He"], geometry=[[1.0, 1.0, 0.0]]),
                    ],
                ),
            ],
            units="angstrom",
        )

        assert hmbe_schema.tiers == 2
        assert hmbe_schema.max_primary_per_nmer == 2
        assert len(hmbe_schema.fragments) == 2
        assert hmbe_schema.units == "angstrom"

    def test_3tier_hierarchy_schema(self):
        """Test valid 3-tier hierarchy schema."""
        hmbe_schema = HMBEHierarchicalMoleculeSchema(
            tiers=3,
            max_primary_per_nmer=2,
            fragments=[
                HMBEFragmentSchema(
                    id="primary",
                    sub_fragments=[
                        HMBEFragmentSchema(
                            id="secondary",
                            sub_fragments=[
                                HMBEFragmentSchema(id="H1", symbols=["H"], geometry=[[0, 0, 0]]),
                                HMBEFragmentSchema(id="H2", symbols=["H"], geometry=[[1, 0, 0]]),
                            ],
                        )
                    ],
                )
            ],
            units="bohr",
        )

        assert hmbe_schema.tiers == 3
        assert hmbe_schema.units == "bohr"

    def test_invalid_tiers(self):
        """Test that tiers < 2 raises validation error."""
        with pytest.raises(ValueError, match="ensure this value is greater than or equal to 2"):
            HMBEHierarchicalMoleculeSchema(
                tiers=1,  # Invalid: must be >= 2
                max_primary_per_nmer=2,
                fragments=[],
            )

    def test_empty_fragments(self):
        """Test that empty fragments list raises validation error."""
        with pytest.raises(ValueError, match="Must specify at least one primary fragment"):
            HMBEHierarchicalMoleculeSchema(
                tiers=2,
                max_primary_per_nmer=2,
                fragments=[],  # Invalid: must have at least one
            )

    def test_geometry_validation(self):
        """Test geometry coordinate validation."""
        # Valid: 3 coordinates per atom
        frag = HMBEFragmentSchema(id="test", symbols=["H", "H"], geometry=[[0, 0, 0], [1, 0, 0]])
        assert len(frag.geometry) == 2

        # Invalid: wrong number of coordinates
        with pytest.raises(ValueError, match="Each coordinate must have exactly 3 values"):
            HMBEFragmentSchema(id="test", symbols=["H"], geometry=[[0, 0]])  # Missing z

    def test_geometry_symbols_mismatch(self):
        """Test that geometry length must match symbols length."""
        with pytest.raises(ValueError, match="Geometry has .* coordinates but .* atoms"):
            HMBEFragmentSchema(
                id="test",
                symbols=["H", "H"],
                geometry=[[0, 0, 0]],  # Only 1 coordinate for 2 atoms
            )


class TestHMBEMoleculeLoading:
    """Test HMBE molecule loading and conversion."""

    @pytest.fixture
    def simple_hmbe_schema(self):
        """Simple 2-tier HMBE molecule schema."""
        return HMBEHierarchicalMoleculeSchema(
            tiers=2,
            max_primary_per_nmer=2,
            fragments=[
                HMBEFragmentSchema(
                    id="cluster_A",
                    sub_fragments=[
                        HMBEFragmentSchema(id="He_A1", symbols=["He"], geometry=[[0.0, 0.0, 0.0]]),
                        HMBEFragmentSchema(id="He_A2", symbols=["He"], geometry=[[1.0, 0.0, 0.0]]),
                    ],
                ),
                HMBEFragmentSchema(
                    id="cluster_B",
                    sub_fragments=[
                        HMBEFragmentSchema(id="He_B1", symbols=["He"], geometry=[[0.0, 1.0, 0.0]]),
                        HMBEFragmentSchema(id="He_B2", symbols=["He"], geometry=[[1.0, 1.0, 0.0]]),
                    ],
                ),
            ],
            units="angstrom",
        )

    def test_convert_hmbe_hierarchical(self, simple_hmbe_schema):
        """Test converting HMBE schema to hierarchy dict."""
        hierarchy_dict = convert_hmbe_hierarchical(simple_hmbe_schema)

        # Check top-level structure
        assert hierarchy_dict["tiers"] == 2
        assert hierarchy_dict["max_primary_per_nmer"] == 2
        assert hierarchy_dict["units"] == "angstrom"
        assert len(hierarchy_dict["fragments"]) == 2

        # Check primary fragment structure
        primary_A = hierarchy_dict["fragments"][0]
        assert primary_A["id"] == "cluster_A"
        assert "sub_fragments" in primary_A
        assert len(primary_A["sub_fragments"]) == 2

        # Check elementary fragment
        elem_A1 = primary_A["sub_fragments"][0]
        assert elem_A1["id"] == "He_A1"
        assert elem_A1["symbols"] == ["He"]
        assert elem_A1["geometry"] == [[0.0, 0.0, 0.0]]

    def test_molecule_schema_with_hmbe(self, simple_hmbe_schema):
        """Test MoleculeSchema with HMBE source."""
        mol_schema = MoleculeSchema(source=MoleculeSourceEnum.hmbe_hierarchical, hmbe=simple_hmbe_schema)

        assert mol_schema.source == MoleculeSourceEnum.hmbe_hierarchical
        assert mol_schema.hmbe is not None

    def test_molecule_schema_hmbe_validation(self):
        """Test that hmbe must be provided when source is hmbe_hierarchical."""
        with pytest.raises(ValueError, match="'hmbe' must be specified"):
            MoleculeSchema(
                source=MoleculeSourceEnum.hmbe_hierarchical,
                hmbe=None,  # Invalid: must provide hmbe when source is hmbe_hierarchical
            )


class TestHMBEConverter:
    """Test converter integration with HMBE."""

    @pytest.fixture
    def hmbe_cli_input(self):
        """Create CLI input with HMBE molecule."""
        return QCManyBodyInput(
            schema_version=1,
            molecule=MoleculeSchema(
                source=MoleculeSourceEnum.hmbe_hierarchical,
                hmbe=HMBEHierarchicalMoleculeSchema(
                    tiers=2,
                    max_primary_per_nmer=1,  # Only within-cluster interactions
                    fragments=[
                        HMBEFragmentSchema(
                            id="cluster_A",
                            sub_fragments=[
                                HMBEFragmentSchema(id="He_A1", symbols=["He"], geometry=[[0.0, 0.0, 0.0]]),
                                HMBEFragmentSchema(id="He_A2", symbols=["He"], geometry=[[1.0, 0.0, 0.0]]),
                            ],
                        ),
                        HMBEFragmentSchema(
                            id="cluster_B",
                            sub_fragments=[
                                HMBEFragmentSchema(id="He_B1", symbols=["He"], geometry=[[0.0, 1.0, 0.0]]),
                                HMBEFragmentSchema(id="He_B2", symbols=["He"], geometry=[[1.0, 1.0, 0.0]]),
                            ],
                        ),
                    ],
                    units="angstrom",
                ),
            ),
            calculation=CalculationSchema(
                type="single",
                single=SingleLevelCalculationSchema(
                    driver=DriverEnum.energy,
                    method="hf",
                    basis="sto-3g",
                    program="psi4",
                ),
            ),
            bsse=BsseSchema(type=[BsseTypeEnum.cp]),
            manybody=ManyBodySchema(max_nbody=2),
        )

    def test_convert_hmbe_to_manybody_input(self, hmbe_cli_input):
        """Test converting HMBE CLI input to ManyBodyInput."""
        mb_input = convert_to_manybody_input(hmbe_cli_input)

        # Check that HMBE hierarchy is stored in keywords
        assert mb_input.specification.keywords.hmbe_hierarchy is not None
        hierarchy = mb_input.specification.keywords.hmbe_hierarchy

        # Check hierarchy structure
        assert hierarchy["tiers"] == 2
        assert hierarchy["max_primary_per_nmer"] == 1
        assert len(hierarchy["fragments"]) == 2

        # Check that geometry was converted to bohr
        assert hierarchy["units"] == "bohr"
        # Original was [[0.0, 0.0, 0.0]] in angstrom
        # Should be converted to bohr (1 angstrom ≈ 1.889726 bohr)
        elem_geom = hierarchy["fragments"][0]["sub_fragments"][0]["geometry"]
        assert elem_geom[0][0] == pytest.approx(0.0)  # Origin stays at origin

        # Check that a placeholder molecule was created
        assert mb_input.molecule is not None
        assert len(mb_input.molecule.symbols) == 1  # Placeholder He atom

        # Check other settings
        assert mb_input.specification.keywords.max_nbody == 2

    def test_hmbe_with_multilevel(self):
        """Test HMBE with multilevel calculation."""
        hmbe_cli_input = QCManyBodyInput(
            schema_version=1,
            molecule=MoleculeSchema(
                source=MoleculeSourceEnum.hmbe_hierarchical,
                hmbe=HMBEHierarchicalMoleculeSchema(
                    tiers=2,
                    max_primary_per_nmer=2,
                    fragments=[
                        HMBEFragmentSchema(
                            id="cluster",
                            sub_fragments=[
                                HMBEFragmentSchema(id="H1", symbols=["H"], geometry=[[0, 0, 0]]),
                                HMBEFragmentSchema(id="H2", symbols=["H"], geometry=[[1, 0, 0]]),
                                HMBEFragmentSchema(id="H3", symbols=["H"], geometry=[[2, 0, 0]]),
                            ],
                        )
                    ],
                    units="bohr",
                ),
            ),
            calculation=CalculationSchema(
                type="multi",
                multi=MultiLevelCalculationSchema(
                    driver=DriverEnum.energy,
                    levels={
                        "1": {"method": "ccsd(t)", "basis": "cc-pvtz", "program": "psi4"},
                        "2": {"method": "mp2", "basis": "cc-pvdz", "program": "psi4"},
                    },
                ),
            ),
            bsse=BsseSchema(type=[BsseTypeEnum.nocp]),
            manybody=ManyBodySchema(max_nbody=2),
        )

        mb_input = convert_to_manybody_input(hmbe_cli_input)

        # Check HMBE hierarchy present
        assert mb_input.specification.keywords.hmbe_hierarchy is not None

        # Check multilevel settings
        assert mb_input.specification.keywords.levels is not None
        assert 1 in mb_input.specification.keywords.levels
        assert 2 in mb_input.specification.keywords.levels


class TestHMBEEndToEnd:
    """Test end-to-end HMBE workflow."""

    def test_hmbe_input_to_preprocessor_ready(self):
        """Test that HMBE CLI input produces preprocessor-ready ManyBodyInput."""
        cli_input = QCManyBodyInput(
            schema_version=1,
            molecule=MoleculeSchema(
                source=MoleculeSourceEnum.hmbe_hierarchical,
                hmbe=HMBEHierarchicalMoleculeSchema(
                    tiers=2,
                    max_primary_per_nmer=2,
                    fragments=[
                        HMBEFragmentSchema(
                            id="water_dimer",
                            sub_fragments=[
                                HMBEFragmentSchema(
                                    id="water1",
                                    symbols=["O", "H", "H"],
                                    geometry=[[0, 0, 0], [1, 0, 0], [-1, 0, 0]],
                                    molecular_charge=0.0,
                                    molecular_multiplicity=1,
                                ),
                                HMBEFragmentSchema(
                                    id="water2",
                                    symbols=["O", "H", "H"],
                                    geometry=[[0, 3, 0], [1, 3, 0], [-1, 3, 0]],
                                    molecular_charge=0.0,
                                    molecular_multiplicity=1,
                                ),
                            ],
                        )
                    ],
                    units="bohr",
                ),
            ),
            calculation=CalculationSchema(
                type="single",
                single=SingleLevelCalculationSchema(
                    driver=DriverEnum.energy,
                    method="hf",
                    basis="sto-3g",
                    program="psi4",
                ),
            ),
            bsse=BsseSchema(type=[BsseTypeEnum.cp]),
            manybody=ManyBodySchema(max_nbody=2),
        )

        mb_input = convert_to_manybody_input(cli_input)

        # Verify HMBE preprocessor can detect this
        from qcmanybody.hmbe.preprocessor import HMBEPreprocessor

        assert HMBEPreprocessor.is_hmbe_input(mb_input)

        # Verify the preprocessor can process it
        # (This will build the flat molecule and generate restricted tuples)
        preprocessed = HMBEPreprocessor.preprocess_hmbe_input(mb_input)

        # Check that preprocessing produced a valid ManyBodyInput
        assert preprocessed.molecule is not None
        # Should have 6 atoms (2 water molecules with 3 atoms each)
        assert len(preprocessed.molecule.symbols) == 6
        # Should have 2 elementary fragments
        assert len(preprocessed.molecule.fragments) == 2

        # Check that HMBE metadata is stored
        assert "hmbe" in preprocessed.extras
        assert preprocessed.extras["hmbe"]["hmbe_mode"] is True

    def test_hmbe_different_units(self):
        """Test HMBE with different unit specifications."""
        # Test with angstrom units (should be converted to bohr)
        cli_input_ang = QCManyBodyInput(
            schema_version=1,
            molecule=MoleculeSchema(
                source=MoleculeSourceEnum.hmbe_hierarchical,
                hmbe=HMBEHierarchicalMoleculeSchema(
                    tiers=2,
                    max_primary_per_nmer=2,
                    fragments=[
                        HMBEFragmentSchema(
                            id="cluster",
                            sub_fragments=[
                                HMBEFragmentSchema(id="He1", symbols=["He"], geometry=[[0.0, 0.0, 0.0]]),
                                HMBEFragmentSchema(id="He2", symbols=["He"], geometry=[[1.0, 0.0, 0.0]]),
                            ],
                        )
                    ],
                    units="angstrom",
                ),
            ),
            calculation=CalculationSchema(
                type="single",
                single=SingleLevelCalculationSchema(
                    driver=DriverEnum.energy,
                    method="hf",
                    basis="sto-3g",
                    program="psi4",
                ),
            ),
            bsse=BsseSchema(type=[BsseTypeEnum.nocp]),
            manybody=ManyBodySchema(max_nbody=2),
        )

        mb_input = convert_to_manybody_input(cli_input_ang)

        # Check that units were converted to bohr
        hierarchy = mb_input.specification.keywords.hmbe_hierarchy
        assert hierarchy["units"] == "bohr"

        # Check that geometry was actually converted (1.0 angstrom ≈ 1.889726 bohr)
        geom = hierarchy["fragments"][0]["sub_fragments"][1]["geometry"]
        assert geom[0][0] == pytest.approx(1.889726, abs=0.001)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
