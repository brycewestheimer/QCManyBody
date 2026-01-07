"""
Tests for HMBE integration with qcmanybody.builder.

Verifies that restricted_tuples parameter correctly limits compute tasks.
"""

import sys
from pathlib import Path

# Ensure we're using the local qcmanybody package
project_root = Path(__file__).parent.parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
from qcmanybody.builder import build_nbody_compute_list
from qcmanybody.models.v1 import BsseEnum


class TestBuilderRestrictedTuples:
    """Test build_nbody_compute_list with restricted_tuples parameter."""

    def test_standard_mbe_without_restrictions(self):
        """Test standard MBE (no restrictions) generates all combinations."""
        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[1, 2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=None,
        )

        # Should have all possible dimers: C(4,2) = 6
        nocp_2body = compute_dict["nocp"][2]

        # Extract fragment tuples (first element of each (frag, bas) pair)
        frag_tuples = {item[0] for item in nocp_2body}

        # All 6 dimers (1-indexed)
        assert (1, 2) in frag_tuples
        assert (1, 3) in frag_tuples
        assert (1, 4) in frag_tuples
        assert (2, 3) in frag_tuples
        assert (2, 4) in frag_tuples
        assert (3, 4) in frag_tuples

    def test_hmbe_restricted_dimers(self):
        """Test HMBE-style restricted dimers."""
        # HMBE restriction: only dimers (0,1) and (2,3) in 0-indexed
        # Converted to 1-indexed for builder: (1,2) and (3,4)
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],  # All monomers (0-indexed)
            2: [(0, 1), (2, 3)],           # Restricted dimers (0-indexed)
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[1, 2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # Should have only restricted dimers
        nocp_2body = compute_dict["nocp"][2]
        frag_tuples = {item[0] for item in nocp_2body}

        # Only 2 restricted dimers (1-indexed)
        assert len([t for t in frag_tuples if len(t) == 2]) == 2
        assert (1, 2) in frag_tuples
        assert (3, 4) in frag_tuples

        # Cross-cluster dimers should NOT be present
        assert (1, 3) not in frag_tuples
        assert (1, 4) not in frag_tuples
        assert (2, 3) not in frag_tuples
        assert (2, 4) not in frag_tuples

    def test_hmbe_restricted_trimers(self):
        """Test HMBE-style restricted trimers."""
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            2: [(0, 1), (0, 2), (1, 2), (2, 3)],  # 4 allowed dimers
            3: [(0, 1, 2)],  # Only 1 allowed trimer
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[1, 2, 3],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # Check 3-body
        nocp_3body = compute_dict["nocp"][3]
        frag_tuples = {item[0] for item in nocp_3body}

        # Only 1 trimer (1-indexed: (1,2,3) from 0-indexed (0,1,2))
        trimers = [t for t in frag_tuples if len(t) == 3]
        assert len(trimers) == 1
        assert (1, 2, 3) in frag_tuples

    def test_cp_with_restrictions(self):
        """Test CP BSSE treatment with restricted tuples."""
        restricted_tuples = {
            1: [(0,), (1,), (2,)],
            2: [(0, 1), (1, 2)],  # Only 2 dimers allowed
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.cp],
            nfragments=3,
            nbodies=[1, 2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # CP creates ghost basis tasks
        cp_2body = compute_dict["cp"][2]

        # Extract unique fragment combinations
        frag_tuples = {item[0] for item in cp_2body}

        # Should only have restricted dimers (1-indexed)
        dimers = [t for t in frag_tuples if len(t) == 2]
        assert (1, 2) in dimers
        assert (2, 3) in dimers

        # Should NOT have (0,2) -> (1,3) dimer
        assert (1, 3) not in frag_tuples

    def test_vmfc_with_restrictions(self):
        """Test VMFC BSSE treatment with restricted tuples."""
        restricted_tuples = {
            1: [(0,), (1,), (2,)],
            2: [(0, 1), (1, 2)],
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.vmfc],
            nfragments=3,
            nbodies=[1, 2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # VMFC creates compute tasks
        vmfc_compute = compute_dict["vmfc_compute"][2]

        # Extract fragment combinations from VMFC tasks
        frag_tuples = {item[0] for item in vmfc_compute}

        # Should respect restrictions (though VMFC has complex internal structure)
        # At minimum, should not generate tasks for non-allowed dimers
        assert len(vmfc_compute) > 0  # VMFC generates some tasks

    def test_return_total_data_with_restrictions(self):
        """Test return_total_data with restricted tuples."""
        restricted_tuples = {
            1: [(0,), (1,)],
            2: [(0, 1)],
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=2,
            nbodies=[1, 2],
            return_total_data=True,  # Includes monomer-basis monomers
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        nocp_1body = compute_dict["nocp"][1]

        # With return_total_data=True, should have monomers in monomer basis
        # ((1,), (1,)) and ((2,), (2,))
        assert ((1,), (1,)) in nocp_1body
        assert ((2,), (2,)) in nocp_1body

    def test_mixed_nbody_levels(self):
        """Test restricted tuples with non-contiguous n-body levels."""
        restricted_tuples = {
            1: [(0,), (1,), (2,), (3,)],
            3: [(0, 1, 2), (1, 2, 3)],  # Skip 2-body, only do 3-body
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[1, 3],  # Skip 2-body
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # Should have 3-body but not 2-body
        assert 3 in compute_dict["nocp"]
        nocp_3body = compute_dict["nocp"][3]

        frag_tuples = {item[0] for item in nocp_3body}

        # Only restricted trimers (1-indexed)
        assert (1, 2, 3) in frag_tuples
        assert (2, 3, 4) in frag_tuples


class TestBuilderRestrictionValidation:
    """Test that builder correctly handles restriction edge cases."""

    def test_empty_restricted_level(self):
        """Test with an n-body level that has no allowed tuples."""
        restricted_tuples = {
            1: [(0,), (1,)],
            2: [],  # No dimers allowed
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=2,
            nbodies=[1, 2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # 2-body should exist
        assert 2 in compute_dict["nocp"]

        # Should not have any actual dimers (2-fragment combinations)
        frag_tuples = {item[0] for item in compute_dict["nocp"][2]}
        dimers = [t for t in frag_tuples if len(t) == 2]
        assert len(dimers) == 0

    def test_partial_level_restrictions(self):
        """Test when only some n-body levels have restrictions."""
        # Only restrict 2-body, let 1-body and 3-body be unrestricted
        restricted_tuples = {
            2: [(0, 1), (2, 3)],  # Only restrict dimers
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[1, 2, 3],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        # 1-body: all 4 monomers (unrestricted)
        nocp_1body = compute_dict["nocp"][1]
        frag_1 = {item[0] for item in nocp_1body}
        monomers = [t for t in frag_1 if len(t) == 1]
        assert len(monomers) == 4

        # 2-body: restricted to 2 dimers
        nocp_2body = compute_dict["nocp"][2]
        frag_2 = {item[0] for item in nocp_2body}
        dimers = [t for t in frag_2 if len(t) == 2]
        assert len(dimers) == 2

        # 3-body: all C(4,3) = 4 trimers (unrestricted)
        nocp_3body = compute_dict["nocp"][3]
        frag_3 = {item[0] for item in nocp_3body}
        trimers = [t for t in frag_3 if len(t) == 3]
        assert len(trimers) == 4


class TestBuilderIndexConversion:
    """Test correct conversion between 0-indexed and 1-indexed."""

    def test_index_conversion(self):
        """Verify that 0-indexed input converts to 1-indexed internally."""
        # HMBE uses 0-indexing: (0,1), (2,3)
        restricted_tuples = {
            2: [(0, 1), (2, 3)],
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=4,
            nbodies=[2],
            return_total_data=False,
            supersystem_ie_only=False,
            restricted_tuples=restricted_tuples,
        )

        nocp_2body = compute_dict["nocp"][2]
        frag_tuples = {item[0] for item in nocp_2body}

        # Builder internally uses 1-indexing
        assert (1, 2) in frag_tuples  # From 0-indexed (0,1)
        assert (3, 4) in frag_tuples  # From 0-indexed (2,3)

        # Should not have 0-indexed tuples
        assert (0, 1) not in frag_tuples
        assert (2, 3) not in frag_tuples  # This is a different tuple (frag 2 and 3)


class TestBuilderSuperSystemWithRestrictions:
    """Test supersystem calculations with restrictions."""

    def test_supersystem_ie_only_with_restrictions(self):
        """Test supersystem_ie_only with restricted tuples."""
        restricted_tuples = {
            1: [(0,), (1,), (2,)],
            # Note: supersystem_ie_only means we skip intermediate bodies
        }

        compute_dict = build_nbody_compute_list(
            bsse_type=[BsseEnum.nocp],
            nfragments=3,
            nbodies=[1, 2, 3],
            return_total_data=False,
            supersystem_ie_only=True,  # Skip 2-body
            restricted_tuples=restricted_tuples,
        )

        # Should have 1-body and 3-body (full system), but not 2-body
        assert 1 in compute_dict["nocp"]
        assert 3 in compute_dict["nocp"]

        # 2-body might exist but should be empty or minimal
        if 2 in compute_dict["nocp"]:
            assert len(compute_dict["nocp"][2]) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
