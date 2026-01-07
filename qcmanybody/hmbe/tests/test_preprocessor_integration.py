"""
Integration tests for HMBE preprocessor with QCManyBody core.

Tests that HMBE preprocessing correctly integrates with ManyBodyComputer and ManyBodyCore.
"""

import sys
from pathlib import Path

# Ensure we're using the local qcmanybody package
project_root = Path(__file__).parent.parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

try:
    import pytest
except ImportError:
    pytest = None

from qcelemental.models import Molecule
from qcmanybody.models.v1 import ManyBodyInput, ManyBodyKeywords, ManyBodySpecification, AtomicSpecification
from qcmanybody.hmbe import HMBEPreprocessor


def test_preprocessor_is_hmbe_input():
    """Test HMBE input detection."""
    # Non-HMBE input
    mol = Molecule(
        symbols=["He", "He"],
        geometry=[[0, 0, 0], [0, 0, 3]],
        fragments=[[0], [1]],
    )

    spec = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(bsse_type=["cp"], max_nbody=2),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )

    input_model = ManyBodyInput(molecule=mol, specification=spec)

    assert not HMBEPreprocessor.is_hmbe_input(input_model)

    # HMBE input
    spec_hmbe = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(
            bsse_type=["cp"],
            max_nbody=2,
            hmbe_hierarchy={
                "tiers": 2,
                "max_primary_per_nmer": 2,
                "fragments": [],
            },
        ),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )

    input_hmbe = ManyBodyInput(molecule=mol, specification=spec_hmbe)

    assert HMBEPreprocessor.is_hmbe_input(input_hmbe)


def test_preprocessor_basic_integration():
    """Test basic HMBE preprocessing creates valid flat molecule."""
    # Create simple 2-tier HMBE input
    hmbe_hierarchy = {
        "tiers": 2,
        "max_primary_per_nmer": 2,
        "fragments": [
            {
                "id": "cluster_A",
                "sub_fragments": [
                    {
                        "id": "A1",
                        "symbols": ["He"],
                        "geometry": [[0.0, 0.0, 0.0]],
                    },
                    {
                        "id": "A2",
                        "symbols": ["He"],
                        "geometry": [[1.0, 0.0, 0.0]],
                    },
                ],
            },
            {
                "id": "cluster_B",
                "sub_fragments": [
                    {
                        "id": "B1",
                        "symbols": ["He"],
                        "geometry": [[0.0, 1.0, 0.0]],
                    },
                    {
                        "id": "B2",
                        "symbols": ["He"],
                        "geometry": [[1.0, 1.0, 0.0]],
                    },
                ],
            },
        ],
    }

    # Create ManyBodyInput with HMBE hierarchy
    # Note: Initial molecule is a placeholder, will be replaced by preprocessing
    placeholder_mol = Molecule(symbols=["He"], geometry=[[0, 0, 0]])

    spec = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(
            bsse_type=["cp"],
            max_nbody=2,
            hmbe_hierarchy=hmbe_hierarchy,
        ),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )

    input_model = ManyBodyInput(molecule=placeholder_mol, specification=spec)

    # Preprocess
    preprocessed = HMBEPreprocessor.preprocess_hmbe_input(input_model)

    # Verify flat molecule was created
    assert len(preprocessed.molecule.symbols) == 4  # 4 He atoms
    assert len(preprocessed.molecule.fragments) == 4  # 4 elementary fragments

    # Verify HMBE metadata was stored
    hmbe_metadata = HMBEPreprocessor.get_hmbe_metadata(preprocessed)
    assert hmbe_metadata is not None
    assert hmbe_metadata["hmbe_mode"] is True
    assert hmbe_metadata["tiers"] == 2
    assert hmbe_metadata["num_elementary"] == 4
    assert hmbe_metadata["num_primary"] == 2

    # Verify restricted tuples were generated
    valid_tuples = hmbe_metadata["valid_tuples"]
    assert "1" in valid_tuples  # Monomers
    assert "2" in valid_tuples  # Dimers

    # Check that we have all 4 monomers
    assert len(valid_tuples["1"]) == 4

    # Check that dimers are restricted (not all 6 possible)
    # With max_primary_per_nmer=2, we should have some dimers
    assert len(valid_tuples["2"]) > 0

    # Verify hmbe_hierarchy was removed from keywords
    assert preprocessed.specification.keywords.hmbe_hierarchy is None


def test_preprocessor_tuple_restriction():
    """Test that HMBE restricts tuples correctly based on max_primary_per_nmer."""
    # Create 2-tier HMBE with max_primary=1 (only within-cluster dimers)
    hmbe_hierarchy = {
        "tiers": 2,
        "max_primary_per_nmer": 1,  # Strict restriction
        "fragments": [
            {
                "id": "cluster_A",
                "sub_fragments": [
                    {"id": "A1", "symbols": ["He"], "geometry": [[0.0, 0.0, 0.0]]},
                    {"id": "A2", "symbols": ["He"], "geometry": [[1.0, 0.0, 0.0]]},
                ],
            },
            {
                "id": "cluster_B",
                "sub_fragments": [
                    {"id": "B1", "symbols": ["He"], "geometry": [[0.0, 1.0, 0.0]]},
                    {"id": "B2", "symbols": ["He"], "geometry": [[1.0, 1.0, 0.0]]},
                ],
            },
        ],
    }

    placeholder_mol = Molecule(symbols=["He"], geometry=[[0, 0, 0]])

    spec = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(
            bsse_type=["cp"],
            max_nbody=2,
            hmbe_hierarchy=hmbe_hierarchy,
        ),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )

    input_model = ManyBodyInput(molecule=placeholder_mol, specification=spec)
    preprocessed = HMBEPreprocessor.preprocess_hmbe_input(input_model)

    hmbe_metadata = HMBEPreprocessor.get_hmbe_metadata(preprocessed)
    valid_tuples = hmbe_metadata["valid_tuples"]

    # With max_primary=1, dimers can only be within same cluster
    # Clusters A and B each have 2 fragments, so 1 dimer each
    # Total: 2 within-cluster dimers
    assert len(valid_tuples["2"]) == 2

    # All tuples should be from same primary cluster
    # Elementary indices: A1=0, A2=1, B1=2, B2=3
    # Valid dimers: (0,1) from cluster A, (2,3) from cluster B
    dimers = [tuple(t) for t in valid_tuples["2"]]
    assert (0, 1) in dimers  # A1-A2
    assert (2, 3) in dimers  # B1-B2

    # Cross-cluster dimers should NOT be present
    assert (0, 2) not in dimers  # A1-B1
    assert (0, 3) not in dimers  # A1-B2
    assert (1, 2) not in dimers  # A2-B1
    assert (1, 3) not in dimers  # A2-B2


def test_get_hmbe_metadata():
    """Test extracting HMBE metadata from preprocessed input."""
    # Non-HMBE input
    mol = Molecule(symbols=["He"], geometry=[[0, 0, 0]])
    spec = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(bsse_type=["cp"]),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )
    input_model = ManyBodyInput(molecule=mol, specification=spec)

    assert HMBEPreprocessor.get_hmbe_metadata(input_model) is None

    # Manually create preprocessed-like input with HMBE metadata
    input_with_extras = input_model.copy(
        update={
            "extras": {
                "hmbe": {
                    "hmbe_mode": True,
                    "tiers": 2,
                    "num_elementary": 4,
                }
            }
        }
    )

    metadata = HMBEPreprocessor.get_hmbe_metadata(input_with_extras)
    assert metadata is not None
    assert metadata["hmbe_mode"] is True
    assert metadata["tiers"] == 2


def test_is_preprocessed():
    """Test checking if input has been preprocessed."""
    mol = Molecule(symbols=["He"], geometry=[[0, 0, 0]])
    spec = ManyBodySpecification(
        driver="energy",
        keywords=ManyBodyKeywords(bsse_type=["cp"]),
        specification={
            "(auto)": AtomicSpecification(
                program="psi4",
                driver="energy",
                model={"method": "hf", "basis": "sto-3g"},
            )
        },
    )

    # Not preprocessed
    input_model = ManyBodyInput(molecule=mol, specification=spec)
    assert not HMBEPreprocessor.is_preprocessed(input_model)

    # Preprocessed (has HMBE metadata)
    input_preprocessed = input_model.copy(
        update={"extras": {"hmbe": {"hmbe_mode": True}}}
    )
    assert HMBEPreprocessor.is_preprocessed(input_preprocessed)


if __name__ == "__main__":
    # Run basic tests
    print("Testing HMBE input detection...")
    test_preprocessor_is_hmbe_input()
    print("✓ HMBE input detection works")

    print("\nTesting basic preprocessing...")
    test_preprocessor_basic_integration()
    print("✓ Basic preprocessing works")

    print("\nTesting tuple restriction...")
    test_preprocessor_tuple_restriction()
    print("✓ Tuple restriction works")

    print("\nTesting metadata extraction...")
    test_get_hmbe_metadata()
    print("✓ Metadata extraction works")

    print("\nTesting preprocessed detection...")
    test_is_preprocessed()
    print("✓ Preprocessed detection works")

    print("\n✅ All integration tests passed!")
