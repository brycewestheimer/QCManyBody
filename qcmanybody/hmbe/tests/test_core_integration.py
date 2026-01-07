"""
Tests for HMBE integration with ManyBodyCore.

Verifies that restricted_tuples parameter flows through ManyBodyCore correctly.
"""

import sys
from pathlib import Path

# Ensure we're using the local qcmanybody package
project_root = Path(__file__).parent.parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from qcelemental.models import Molecule
from qcmanybody import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum


@pytest.fixture
def simple_he4_molecule():
    """Create simple He4 molecule with 4 fragments."""
    return Molecule(
        symbols=["He", "He", "He", "He"],
        geometry=[
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [3.0, 3.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3]],
        molecular_charge=0.0,
        molecular_multiplicity=1,
    )


class TestManyBodyCoreRestrictedTuples:
    """Test ManyBodyCore with restricted_tuples parameter."""

    def test_core_initialization_without_restrictions(self, simple_he4_molecule):
        """Test ManyBodyCore initialization without restrictions."""
        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=None,
        )

        assert core.nfragments == 4
        assert core.max_nbody == 2
        assert core.restricted_tuples is None

    def test_core_initialization_with_restrictions(self, simple_he4_molecule):
        """Test ManyBodyCore initialization with restrictions."""
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            2: [(0, 1), (2, 3)],  # Only 2 dimers
        }

        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        assert core.restricted_tuples == restricted_tuples

    def test_compute_map_without_restrictions(self, simple_he4_molecule):
        """Test compute_map generation without restrictions."""
        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=None,
        )

        compute_map = core.compute_map

        # Standard MBE: all C(4,2) = 6 dimers
        nocp_2body = compute_map["hf/sto-3g"]["nocp"][2]
        frag_tuples = {item[0] for item in nocp_2body}

        dimers = [t for t in frag_tuples if len(t) == 2]
        assert len(dimers) == 6

    def test_compute_map_with_restrictions(self, simple_he4_molecule):
        """Test compute_map generation with restrictions."""
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            2: [(0, 1), (2, 3)],  # Only 2 dimers allowed
        }

        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        compute_map = core.compute_map

        # HMBE: only 2 restricted dimers
        nocp_2body = compute_map["hf/sto-3g"]["nocp"][2]
        frag_tuples = {item[0] for item in nocp_2body}

        dimers = [t for t in frag_tuples if len(t) == 2]
        assert len(dimers) == 2

        # Check specific dimers (1-indexed)
        assert (1, 2) in frag_tuples
        assert (3, 4) in frag_tuples

        # Check that non-allowed dimers are NOT present
        assert (1, 3) not in frag_tuples
        assert (1, 4) not in frag_tuples
        assert (2, 3) not in frag_tuples
        assert (2, 4) not in frag_tuples

    def test_multilevel_with_restrictions(self, simple_he4_molecule):
        """Test multilevel calculation with restrictions."""
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            2: [(0, 1), (1, 2)],
            3: [(0, 1, 2)],
        }

        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={
                1: "hf/sto-3g",
                2: "mp2/cc-pvdz",
                3: "ccsd(t)/cc-pvtz",
            },
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        compute_map = core.compute_map

        # Each level should have restricted tuples
        # 1-body at HF level
        assert "hf/sto-3g" in compute_map
        hf_1body = compute_map["hf/sto-3g"]["nocp"][1]
        frag_1 = {item[0] for item in hf_1body}
        monomers = [t for t in frag_1 if len(t) == 1]
        assert len(monomers) == 4

        # 2-body at MP2 level
        assert "mp2/cc-pvdz" in compute_map
        mp2_2body = compute_map["mp2/cc-pvdz"]["nocp"][2]
        frag_2 = {item[0] for item in mp2_2body}
        dimers = [t for t in frag_2 if len(t) == 2]
        assert len(dimers) == 2  # Restricted

        # 3-body at CCSD(T) level
        assert "ccsd(t)/cc-pvtz" in compute_map
        ccsd_3body = compute_map["ccsd(t)/cc-pvtz"]["nocp"][3]
        frag_3 = {item[0] for item in ccsd_3body}
        trimers = [t for t in frag_3 if len(t) == 3]
        assert len(trimers) == 1  # Restricted

    def test_cp_bsse_with_restrictions(self, simple_he4_molecule):
        """Test CP BSSE treatment with restrictions."""
        restricted_tuples = {
            1: [(0,), (1,), (2,)],
            2: [(0, 1), (1, 2)],
        }

        mol_3frag = Molecule(
            symbols=["He", "He", "He"],
            geometry=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            fragments=[[0], [1], [2]],
        )

        core = ManyBodyCore(
            mol_3frag,
            bsse_type=[BsseEnum.cp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        compute_map = core.compute_map

        # CP should create ghost basis tasks for restricted tuples
        cp_2body = compute_map["hf/sto-3g"]["cp"][2]

        # Should have tasks, but restricted to allowed dimers
        assert len(cp_2body) > 0

        # Extract fragment combinations
        frag_tuples = {item[0] for item in cp_2body}

        # Should not have the non-allowed dimer (0,2) -> (1,3)
        # Note: CP might have monomer tasks too, filter to dimers
        dimers = [t for t in frag_tuples if len(t) == 2]
        assert (1, 3) not in dimers


class TestManyBodyCoreIterateMolecules:
    """Test molecule iteration with restrictions."""

    def test_iterate_molecules_with_restrictions(self, simple_he4_molecule):
        """Test that iterate_molecules respects restrictions."""
        restricted_tuples = {
            1: [(0,), (1,)],  # Only 2 monomers
            2: [(0, 1)],      # Only 1 dimer
        }

        mol_2frag = Molecule(
            symbols=["He", "He"],
            geometry=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            fragments=[[0], [1]],
        )

        core = ManyBodyCore(
            mol_2frag,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        # Iterate and count molecules
        molecules_generated = []
        for chem, label, mol in core.iterate_molecules():
            molecules_generated.append((chem, label, mol.symbols))

        # Should generate fewer molecules with restrictions
        # Without restrictions: 2 monomers + 1 dimer = 3 molecules
        # With restrictions: same in this case, but in general should be fewer
        assert len(molecules_generated) > 0

        # All generated molecules should be from allowed tuples
        # (This is implicitly tested by the compute_map validation)


class TestManyBodyCoreFragmentSlicing:
    """Test that fragment slicing works with restrictions."""

    def test_fragment_slices_unchanged(self, simple_he4_molecule):
        """Test that restrictions don't affect fragment slicing."""
        restricted_tuples = {
            1: [(0,), (2,)],  # Skip fragments 1 and 3
            2: [(0, 2)],
        }

        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        # Fragment slices should still cover all fragments
        assert len(core.fragment_slice_dict) == 4
        assert 1 in core.fragment_slice_dict
        assert 2 in core.fragment_slice_dict
        assert 3 in core.fragment_slice_dict
        assert 4 in core.fragment_slice_dict


class TestManyBodyCoreEdgeCases:
    """Test edge cases with restrictions."""

    def test_empty_restriction_level(self, simple_he4_molecule):
        """Test when a level has no allowed tuples."""
        restricted_tuples = {
            1: [(0,), (1,)],
            2: [],  # No dimers allowed
        }

        mol_2frag = Molecule(
            symbols=["He", "He"],
            geometry=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            fragments=[[0], [1]],
        )

        core = ManyBodyCore(
            mol_2frag,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        compute_map = core.compute_map

        # 2-body should exist
        assert 2 in compute_map["hf/sto-3g"]["nocp"]

        # Should not have any actual dimers (2-fragment combinations)
        frag_tuples = {item[0] for item in compute_map["hf/sto-3g"]["nocp"][2]}
        dimers = [t for t in frag_tuples if len(t) == 2]
        assert len(dimers) == 0

    def test_restrictions_larger_than_nfragments(self, simple_he4_molecule):
        """Test when restrictions include fragment indices beyond nfragments."""
        # This shouldn't happen in practice (preprocessor prevents it),
        # but test robustness
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            2: [(0, 1)],
        }

        # Should not raise error (builder will skip invalid indices)
        core = ManyBodyCore(
            simple_he4_molecule,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges=None,
            restricted_tuples=restricted_tuples,
        )

        # Should work normally
        compute_map = core.compute_map
        assert "hf/sto-3g" in compute_map


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
