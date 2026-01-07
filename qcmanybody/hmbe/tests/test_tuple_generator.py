"""
Comprehensive tests for HMBE tuple generator.

Tests the HMBETupleGenerator class which generates valid N-mer tuples
respecting hierarchical constraints.
"""

import sys
from pathlib import Path

# Ensure we're using the local qcmanybody package
project_root = Path(__file__).parent.parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
from qcmanybody.hmbe import FragmentHierarchy, HMBETupleGenerator


@pytest.fixture
def simple_2tier_hierarchy():
    """Create a simple 2-tier hierarchy for testing.

    Structure:
        Cluster A: elem_0, elem_1
        Cluster B: elem_2, elem_3
    """
    data = {
        "tiers": 2,
        "max_primary_per_nmer": 2,
        "fragments": [
            {
                "id": "cluster_A",
                "sub_fragments": [
                    {"id": "elem_0", "symbols": ["He"], "geometry": [[0, 0, 0]]},
                    {"id": "elem_1", "symbols": ["He"], "geometry": [[1, 0, 0]]},
                ],
            },
            {
                "id": "cluster_B",
                "sub_fragments": [
                    {"id": "elem_2", "symbols": ["He"], "geometry": [[0, 1, 0]]},
                    {"id": "elem_3", "symbols": ["He"], "geometry": [[1, 1, 0]]},
                ],
            },
        ],
    }
    return FragmentHierarchy.from_dict(data)


@pytest.fixture
def three_cluster_hierarchy():
    """Create 3-cluster 2-tier hierarchy for testing.

    Structure:
        Cluster A: elem_0, elem_1
        Cluster B: elem_2
        Cluster C: elem_3, elem_4, elem_5
    """
    data = {
        "tiers": 2,
        "max_primary_per_nmer": 2,
        "fragments": [
            {
                "id": "cluster_A",
                "sub_fragments": [
                    {"id": "elem_0", "symbols": ["He"], "geometry": [[0, 0, 0]]},
                    {"id": "elem_1", "symbols": ["He"], "geometry": [[1, 0, 0]]},
                ],
            },
            {
                "id": "cluster_B",
                "sub_fragments": [
                    {"id": "elem_2", "symbols": ["Ne"], "geometry": [[0, 1, 0]]},
                ],
            },
            {
                "id": "cluster_C",
                "sub_fragments": [
                    {"id": "elem_3", "symbols": ["Ar"], "geometry": [[0, 0, 1]]},
                    {"id": "elem_4", "symbols": ["Ar"], "geometry": [[1, 0, 1]]},
                    {"id": "elem_5", "symbols": ["Ar"], "geometry": [[2, 0, 1]]},
                ],
            },
        ],
    }
    return FragmentHierarchy.from_dict(data)


@pytest.fixture
def three_tier_hierarchy():
    """Create 3-tier hierarchy for testing.

    Structure:
        Primary 0:
            Secondary 0,0: elem_0
            Secondary 0,1: elem_1, elem_2
        Primary 1:
            Secondary 1,0: elem_3
    """
    data = {
        "tiers": 3,
        "max_primary_per_nmer": 2,
        "fragments": [
            {
                "id": "primary_0",
                "sub_fragments": [
                    {
                        "id": "secondary_0_0",
                        "sub_fragments": [
                            {"id": "elem_0", "symbols": ["He"], "geometry": [[0, 0, 0]]},
                        ],
                    },
                    {
                        "id": "secondary_0_1",
                        "sub_fragments": [
                            {"id": "elem_1", "symbols": ["He"], "geometry": [[1, 0, 0]]},
                            {"id": "elem_2", "symbols": ["He"], "geometry": [[2, 0, 0]]},
                        ],
                    },
                ],
            },
            {
                "id": "primary_1",
                "sub_fragments": [
                    {
                        "id": "secondary_1_0",
                        "sub_fragments": [
                            {"id": "elem_3", "symbols": ["He"], "geometry": [[0, 1, 0]]},
                        ],
                    },
                ],
            },
        ],
    }
    return FragmentHierarchy.from_dict(data)


class TestHMBETupleGeneratorBasic:
    """Basic functionality tests."""

    def test_initialization(self, simple_2tier_hierarchy):
        """Test generator initialization."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        assert generator.hierarchy == simple_2tier_hierarchy
        assert generator.max_primary == 2
        assert len(generator._elementary_to_primary) == 4

    def test_elementary_to_primary_mapping(self, simple_2tier_hierarchy):
        """Test that elementary fragment to primary mapping is correct."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        # elem_0 and elem_1 belong to cluster_A (primary index 0)
        assert generator._elementary_to_primary[0] == 0
        assert generator._elementary_to_primary[1] == 0

        # elem_2 and elem_3 belong to cluster_B (primary index 1)
        assert generator._elementary_to_primary[2] == 1
        assert generator._elementary_to_primary[3] == 1

    def test_monomer_generation(self, simple_2tier_hierarchy):
        """Test that all monomers are generated."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)
        tuples = generator.generate_nmer_tuples(max_nbody=1)

        assert 1 in tuples
        monomers = tuples[1]

        # Should have all 4 elementary fragments as monomers
        assert len(monomers) == 4
        assert (0,) in monomers
        assert (1,) in monomers
        assert (2,) in monomers
        assert (3,) in monomers


class TestHMBETupleGeneratorRestrictions:
    """Test tuple restriction based on max_primary_per_nmer."""

    def test_max_primary_1_dimers(self, simple_2tier_hierarchy):
        """With max_primary=1, only within-cluster dimers are allowed."""
        # Modify hierarchy to have max_primary=1
        simple_2tier_hierarchy.max_primary_per_nmer = 1
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=2)
        dimers = tuples[2]

        # Only 2 within-cluster dimers: (0,1) from A, (2,3) from B
        assert len(dimers) == 2
        assert (0, 1) in dimers  # A: elem_0, elem_1
        assert (2, 3) in dimers  # B: elem_2, elem_3

    def test_max_primary_2_dimers(self, simple_2tier_hierarchy):
        """With max_primary=2, cross-cluster dimers are also allowed."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)
        tuples = generator.generate_nmer_tuples(max_nbody=2)
        dimers = tuples[2]

        # All 6 possible dimers: C(4,2) = 6
        assert len(dimers) == 6

        # Within-cluster dimers
        assert (0, 1) in dimers
        assert (2, 3) in dimers

        # Cross-cluster dimers
        assert (0, 2) in dimers
        assert (0, 3) in dimers
        assert (1, 2) in dimers
        assert (1, 3) in dimers

    def test_max_primary_1_trimers(self, three_cluster_hierarchy):
        """With max_primary=1, only within-cluster trimers are allowed."""
        three_cluster_hierarchy.max_primary_per_nmer = 1
        generator = HMBETupleGenerator(three_cluster_hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=3)
        trimers = tuples[3]

        # Only cluster C has ≥3 fragments: elem_3, elem_4, elem_5
        assert len(trimers) == 1
        assert (3, 4, 5) in trimers

    def test_max_primary_2_trimers(self, three_cluster_hierarchy):
        """With max_primary=2, cross-cluster trimers from ≤2 clusters allowed."""
        generator = HMBETupleGenerator(three_cluster_hierarchy)
        tuples = generator.generate_nmer_tuples(max_nbody=3)
        trimers = tuples[3]

        # Cluster A (2 elem) + Cluster B (1 elem): C(2,2)*C(1,1) = 1 trimer
        # Cluster A (2 elem) + Cluster C (3 elem): C(2,2)*C(3,1) + C(2,1)*C(3,2) = 3 + 6 = 9 trimers
        # Cluster B (1 elem) + Cluster C (3 elem): C(1,1)*C(3,2) = 3 trimers
        # Within Cluster C (3 elem): C(3,3) = 1 trimer
        # Total: 1 + 9 + 3 + 1 = 14 trimers
        assert len(trimers) == 14

    def test_max_primary_3_allows_all(self, three_cluster_hierarchy):
        """With max_primary=3, all trimers from all clusters allowed."""
        three_cluster_hierarchy.max_primary_per_nmer = 3
        generator = HMBETupleGenerator(three_cluster_hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=3)
        trimers = tuples[3]

        # All combinations from 6 elementary fragments: C(6,3) = 20
        assert len(trimers) == 20


class TestHMBETupleGeneratorNBody:
    """Test generation for different n-body levels."""

    def test_incremental_nbody_generation(self, simple_2tier_hierarchy):
        """Test generating tuples up to different max_nbody levels."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        # Generate up to 2-body
        tuples_2 = generator.generate_nmer_tuples(max_nbody=2)
        assert 1 in tuples_2
        assert 2 in tuples_2
        assert 3 not in tuples_2

        # Generate up to 3-body
        tuples_3 = generator.generate_nmer_tuples(max_nbody=3)
        assert 1 in tuples_3
        assert 2 in tuples_3
        assert 3 in tuples_3
        assert 4 not in tuples_3

    def test_max_nbody_exceeds_num_elementary(self, simple_2tier_hierarchy):
        """Test that max_nbody is capped at num_elementary."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        # Request max_nbody=10, but only 4 elementary fragments
        tuples = generator.generate_nmer_tuples(max_nbody=10)

        # Should only generate up to 4-body
        assert max(tuples.keys()) == 4

    def test_4body_generation(self, simple_2tier_hierarchy):
        """Test 4-body tuple generation."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)
        tuples = generator.generate_nmer_tuples(max_nbody=4)

        fourbody = tuples[4]

        # Only one 4-body: all fragments
        # With max_primary=2, this is valid (uses 2 clusters)
        assert len(fourbody) == 1
        assert (0, 1, 2, 3) in fourbody


class TestHMBETupleGeneratorThreeTier:
    """Test with 3-tier hierarchies."""

    def test_three_tier_primary_mapping(self, three_tier_hierarchy):
        """Test primary mapping for 3-tier hierarchy."""
        generator = HMBETupleGenerator(three_tier_hierarchy)

        # elem_0, elem_1, elem_2 belong to primary_0
        assert generator._elementary_to_primary[0] == 0
        assert generator._elementary_to_primary[1] == 0
        assert generator._elementary_to_primary[2] == 0

        # elem_3 belongs to primary_1
        assert generator._elementary_to_primary[3] == 1

    def test_three_tier_dimers_max_primary_1(self, three_tier_hierarchy):
        """Test 3-tier dimers with max_primary=1."""
        three_tier_hierarchy.max_primary_per_nmer = 1
        generator = HMBETupleGenerator(three_tier_hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=2)
        dimers = tuples[2]

        # Within primary_0: C(3,2) = 3 dimers
        # Within primary_1: C(1,2) = 0 dimers
        assert len(dimers) == 3
        assert (0, 1) in dimers
        assert (0, 2) in dimers
        assert (1, 2) in dimers

    def test_three_tier_dimers_max_primary_2(self, three_tier_hierarchy):
        """Test 3-tier dimers with max_primary=2."""
        generator = HMBETupleGenerator(three_tier_hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=2)
        dimers = tuples[2]

        # All combinations from 4 elementary fragments: C(4,2) = 6
        assert len(dimers) == 6


class TestHMBETupleGeneratorHelpers:
    """Test helper methods."""

    def test_is_valid_hmbe_tuple(self, simple_2tier_hierarchy):
        """Test tuple validation logic."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        # Within cluster A (primary 0): valid
        assert generator._is_valid_hmbe_tuple((0, 1))

        # Within cluster B (primary 1): valid
        assert generator._is_valid_hmbe_tuple((2, 3))

        # Cross-cluster with max_primary=2: valid
        assert generator._is_valid_hmbe_tuple((0, 2))

        # Now restrict to max_primary=1
        simple_2tier_hierarchy.max_primary_per_nmer = 1
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        # Cross-cluster with max_primary=1: invalid
        assert not generator._is_valid_hmbe_tuple((0, 2))
        assert not generator._is_valid_hmbe_tuple((0, 3))

    def test_get_tuple_info(self, simple_2tier_hierarchy):
        """Test detailed tuple information retrieval."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        info = generator.get_tuple_info((0, 2))

        assert info["nbody"] == 2
        assert info["elementary_indices"] == [0, 2]
        assert info["primary_indices"] == [0, 1]  # From two different clusters
        assert info["num_primary"] == 2
        assert info["is_valid"] is True  # max_primary=2
        assert "elem_0" in info["fragment_ids"]
        assert "elem_2" in info["fragment_ids"]

    def test_count_tuples_by_level(self, simple_2tier_hierarchy):
        """Test counting tuples without generating them."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        counts = generator.count_tuples_by_level(max_nbody=3)

        assert counts[1] == 4  # 4 monomers
        assert counts[2] == 6  # 6 dimers (max_primary=2)
        assert counts[3] == 4  # 4 trimers

    def test_generate_tuples_by_primary(self, three_cluster_hierarchy):
        """Test organizing tuples by primary combination."""
        generator = HMBETupleGenerator(three_cluster_hierarchy)

        organized = generator.generate_tuples_by_primary(max_nbody=2)

        # Should have dimers organized by primary combination
        dimers_by_primary = organized[2]

        # Within cluster A (primary 0)
        assert (0,) in dimers_by_primary
        assert (0, 1) in dimers_by_primary[(0,)]

        # Within cluster C (primary 2)
        assert (2,) in dimers_by_primary
        dimers_from_C = dimers_by_primary[(2,)]
        assert (3, 4) in dimers_from_C
        assert (3, 5) in dimers_from_C
        assert (4, 5) in dimers_from_C

        # Cross-cluster A+B (primary 0,1)
        assert (0, 1) in dimers_by_primary


class TestHMBETupleGeneratorEdgeCases:
    """Test edge cases and error handling."""

    def test_max_nbody_less_than_1(self, simple_2tier_hierarchy):
        """Test that max_nbody < 1 raises error."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)

        with pytest.raises(ValueError, match="max_nbody must be >= 1"):
            generator.generate_nmer_tuples(max_nbody=0)

    def test_single_fragment_cluster(self):
        """Test with clusters containing single fragments."""
        data = {
            "tiers": 2,
            "max_primary_per_nmer": 2,
            "fragments": [
                {
                    "id": "cluster_A",
                    "sub_fragments": [
                        {"id": "elem_0", "symbols": ["He"], "geometry": [[0, 0, 0]]},
                    ],
                },
                {
                    "id": "cluster_B",
                    "sub_fragments": [
                        {"id": "elem_1", "symbols": ["Ne"], "geometry": [[1, 0, 0]]},
                    ],
                },
            ],
        }
        hierarchy = FragmentHierarchy.from_dict(data)
        generator = HMBETupleGenerator(hierarchy)

        tuples = generator.generate_nmer_tuples(max_nbody=2)

        # 2 monomers
        assert len(tuples[1]) == 2

        # 1 cross-cluster dimer (since each cluster has only 1 fragment)
        assert len(tuples[2]) == 1
        assert (0, 1) in tuples[2]

    def test_tuple_sorting(self, simple_2tier_hierarchy):
        """Test that tuples are sorted and unique."""
        generator = HMBETupleGenerator(simple_2tier_hierarchy)
        tuples = generator.generate_nmer_tuples(max_nbody=3)

        for nbody, tuple_list in tuples.items():
            for t in tuple_list:
                # Each tuple should be sorted
                assert list(t) == sorted(t)

            # No duplicates
            assert len(tuple_list) == len(set(tuple_list))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
