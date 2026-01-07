"""
Tests for hierarchical fragment representation.
"""

import numpy as np
import pytest

from qcmanybody.hmbe.hierarchy import FragmentHierarchy, HierarchicalFragment


class TestHierarchicalFragment:
    """Tests for HierarchicalFragment class."""

    def test_elementary_fragment_creation(self):
        """Test creating an elementary fragment with atomic data."""
        frag = HierarchicalFragment(
            fragment_id="He1",
            tier=2,
            parent_path=(0, 0),
            symbols=["He"],
            geometry=np.array([[0.0, 0.0, 0.0]]),
            molecular_charge=0.0,
            molecular_multiplicity=1,
        )

        assert frag.fragment_id == "He1"
        assert frag.tier == 2
        assert frag.parent_path == (0, 0)
        assert frag.is_elementary
        assert not frag.is_composite
        assert frag.natoms == 1
        assert frag.get_elementary_count() == 1

    def test_composite_fragment_creation(self):
        """Test creating a composite fragment with sub-fragments."""
        sub1 = HierarchicalFragment(
            fragment_id="He1",
            tier=2,
            parent_path=(0, 0),
            symbols=["He"],
            geometry=np.array([[0.0, 0.0, 0.0]]),
        )
        sub2 = HierarchicalFragment(
            fragment_id="He2",
            tier=2,
            parent_path=(0, 1),
            symbols=["He"],
            geometry=np.array([[1.0, 0.0, 0.0]]),
        )

        composite = HierarchicalFragment(
            fragment_id="cluster_A",
            tier=1,
            parent_path=(0,),
            sub_fragments=[sub1, sub2],
        )

        assert composite.fragment_id == "cluster_A"
        assert composite.tier == 1
        assert not composite.is_elementary
        assert composite.is_composite
        assert composite.get_elementary_count() == 2
        assert composite.get_total_atoms() == 2

    def test_elementary_fragment_iteration(self):
        """Test iterating over elementary fragments."""
        sub1 = HierarchicalFragment(
            fragment_id="He1", tier=2, parent_path=(0, 0), symbols=["He"], geometry=np.array([[0, 0, 0]])
        )
        sub2 = HierarchicalFragment(
            fragment_id="He2", tier=2, parent_path=(0, 1), symbols=["He"], geometry=np.array([[1, 0, 0]])
        )
        composite = HierarchicalFragment(
            fragment_id="cluster", tier=1, parent_path=(0,), sub_fragments=[sub1, sub2]
        )

        elementary_frags = list(composite.iter_elementary_fragments())
        assert len(elementary_frags) == 2
        assert elementary_frags[0].fragment_id == "He1"
        assert elementary_frags[1].fragment_id == "He2"

    def test_nested_hierarchy_iteration(self):
        """Test iteration over multi-level nested hierarchy."""
        # Create 3-tier structure: primary -> secondary -> tertiary (elementary)
        tert1 = HierarchicalFragment(
            fragment_id="atom1",
            tier=3,
            parent_path=(0, 0, 0),
            symbols=["H"],
            geometry=np.array([[0, 0, 0]]),
        )
        tert2 = HierarchicalFragment(
            fragment_id="atom2",
            tier=3,
            parent_path=(0, 0, 1),
            symbols=["H"],
            geometry=np.array([[1, 0, 0]]),
        )

        sec = HierarchicalFragment(
            fragment_id="secondary", tier=2, parent_path=(0, 0), sub_fragments=[tert1, tert2]
        )

        prim = HierarchicalFragment(fragment_id="primary", tier=1, parent_path=(0,), sub_fragments=[sec])

        elementary_frags = list(prim.iter_elementary_fragments())
        assert len(elementary_frags) == 2
        assert all(f.tier == 3 for f in elementary_frags)

    def test_invalid_elementary_fragment(self):
        """Test that elementary fragments must have symbols and geometry."""
        with pytest.raises(ValueError, match="must have symbols and geometry"):
            HierarchicalFragment(
                fragment_id="bad", tier=1, parent_path=(0,), symbols=None, geometry=None
            )

    def test_invalid_geometry_shape(self):
        """Test that geometry must be (natoms, 3)."""
        with pytest.raises(ValueError, match="geometry must be shape"):
            HierarchicalFragment(
                fragment_id="bad",
                tier=1,
                parent_path=(0,),
                symbols=["He"],
                geometry=np.array([0, 0, 0]),  # Wrong shape (should be [[0,0,0]])
            )

    def test_mismatched_symbols_geometry(self):
        """Test that symbols and geometry lengths must match."""
        with pytest.raises(ValueError, match="has 1 symbols but 2 coordinates"):
            HierarchicalFragment(
                fragment_id="bad",
                tier=1,
                parent_path=(0,),
                symbols=["He"],
                geometry=np.array([[0, 0, 0], [1, 1, 1]]),  # 2 atoms but 1 symbol
            )


class TestFragmentHierarchy:
    """Tests for FragmentHierarchy class."""

    def test_simple_2tier_hierarchy(self):
        """Test creating a simple 2-tier hierarchy."""
        # 2 primary fragments, each with 2 elementary sub-fragments
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="prim_0",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem_00",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["He"],
                            geometry=np.array([[0, 0, 0]]),
                        ),
                        HierarchicalFragment(
                            fragment_id="elem_01",
                            tier=2,
                            parent_path=(0, 1),
                            symbols=["He"],
                            geometry=np.array([[1, 0, 0]]),
                        ),
                    ],
                ),
                HierarchicalFragment(
                    fragment_id="prim_1",
                    tier=1,
                    parent_path=(1,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem_10",
                            tier=2,
                            parent_path=(1, 0),
                            symbols=["He"],
                            geometry=np.array([[0, 1, 0]]),
                        ),
                        HierarchicalFragment(
                            fragment_id="elem_11",
                            tier=2,
                            parent_path=(1, 1),
                            symbols=["He"],
                            geometry=np.array([[1, 1, 0]]),
                        ),
                    ],
                ),
            ],
        )

        assert hierarchy.tiers == 2
        assert hierarchy.num_primary == 2
        assert hierarchy.num_elementary == 4
        assert hierarchy.total_atoms == 4
        assert hierarchy.max_primary_per_nmer == 2

    def test_elementary_fragment_access(self):
        """Test accessing elementary fragments by flat index."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="prim",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem_0",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["He"],
                            geometry=np.array([[0, 0, 0]]),
                        ),
                        HierarchicalFragment(
                            fragment_id="elem_1",
                            tier=2,
                            parent_path=(0, 1),
                            symbols=["He"],
                            geometry=np.array([[1, 0, 0]]),
                        ),
                    ],
                )
            ],
        )

        elem0 = hierarchy.get_elementary_fragment(0)
        assert elem0.fragment_id == "elem_0"
        assert elem0.parent_path == (0, 0)

        elem1 = hierarchy.get_elementary_fragment(1)
        assert elem1.fragment_id == "elem_1"
        assert elem1.parent_path == (0, 1)

    def test_path_based_access(self):
        """Test accessing fragments by hierarchical path."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="prim",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["He"],
                            geometry=np.array([[0, 0, 0]]),
                        )
                    ],
                )
            ],
        )

        frag = hierarchy.get_fragment_by_path((0, 0))
        assert frag.fragment_id == "elem"

        flat_idx = hierarchy.get_flat_index((0, 0))
        assert flat_idx == 0

    def test_primary_index_extraction(self):
        """Test getting primary index from elementary path."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="prim_0",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["He"],
                            geometry=np.array([[0, 0, 0]]),
                        )
                    ],
                ),
                HierarchicalFragment(
                    fragment_id="prim_1",
                    tier=1,
                    parent_path=(1,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="elem",
                            tier=2,
                            parent_path=(1, 0),
                            symbols=["He"],
                            geometry=np.array([[1, 0, 0]]),
                        )
                    ],
                ),
            ],
        )

        assert hierarchy.get_primary_index((0, 0)) == 0
        assert hierarchy.get_primary_index((1, 0)) == 1

    def test_from_dict_simple(self):
        """Test constructing hierarchy from dictionary."""
        data = {
            "tiers": 2,
            "max_primary_per_nmer": 2,
            "fragments": [
                {
                    "id": "primary_1",
                    "sub_fragments": [
                        {"id": "elem_1", "symbols": ["He"], "geometry": [[0, 0, 0]]}
                    ],
                }
            ],
        }

        hierarchy = FragmentHierarchy.from_dict(data)
        assert hierarchy.tiers == 2
        assert hierarchy.num_primary == 1
        assert hierarchy.num_elementary == 1

    def test_from_dict_complex(self):
        """Test constructing complex hierarchy from dictionary."""
        data = {
            "tiers": 3,
            "max_primary_per_nmer": 2,
            "fragments": [
                {
                    "id": "prim_0",
                    "sub_fragments": [
                        {
                            "id": "sec_00",
                            "sub_fragments": [
                                {"id": "elem_000", "symbols": ["H"], "geometry": [[0, 0, 0]]},
                                {"id": "elem_001", "symbols": ["H"], "geometry": [[1, 0, 0]]},
                            ],
                        },
                        {
                            "id": "sec_01",
                            "sub_fragments": [
                                {"id": "elem_010", "symbols": ["O"], "geometry": [[2, 0, 0]]}
                            ],
                        },
                    ],
                },
                {
                    "id": "prim_1",
                    "sub_fragments": [
                        {
                            "id": "sec_10",
                            "sub_fragments": [
                                {"id": "elem_100", "symbols": ["H"], "geometry": [[0, 1, 0]]}
                            ],
                        }
                    ],
                },
            ],
        }

        hierarchy = FragmentHierarchy.from_dict(data)
        assert hierarchy.tiers == 3
        assert hierarchy.num_primary == 2
        assert hierarchy.num_elementary == 4

        # Verify paths
        elem000 = hierarchy.get_fragment_by_path((0, 0, 0))
        assert elem000.fragment_id == "elem_000"
        assert elem000.tier == 3

    def test_invalid_tiers(self):
        """Test that hierarchies must have at least 2 tiers."""
        with pytest.raises(ValueError, match="at least 2 tiers"):
            FragmentHierarchy(tiers=1, max_primary_per_nmer=2, primary_fragments=[])

    def test_invalid_max_primary(self):
        """Test that max_primary_per_nmer must be >= 1."""
        with pytest.raises(ValueError, match="must be >= 1"):
            FragmentHierarchy(tiers=2, max_primary_per_nmer=0, primary_fragments=[])

    def test_empty_primary_fragments(self):
        """Test that at least one primary fragment is required."""
        with pytest.raises(ValueError, match="at least one primary fragment"):
            FragmentHierarchy(tiers=2, max_primary_per_nmer=2, primary_fragments=[])
