"""
Tests for HMBE fragment builder.
"""

import numpy as np
import pytest

from qcmanybody.hmbe.fragment_builder import HMBEFragmentBuilder
from qcmanybody.hmbe.hierarchy import FragmentHierarchy, HierarchicalFragment


class TestHMBEFragmentBuilder:
    """Tests for HMBEFragmentBuilder class."""

    @pytest.fixture
    def simple_2tier_hierarchy(self):
        """Create a simple 2-tier hierarchy for testing."""
        return FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="cluster_A",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="He_A1",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["He"],
                            geometry=np.array([[0.0, 0.0, 0.0]]),
                            molecular_charge=0.0,
                            molecular_multiplicity=1,
                        ),
                        HierarchicalFragment(
                            fragment_id="He_A2",
                            tier=2,
                            parent_path=(0, 1),
                            symbols=["He"],
                            geometry=np.array([[1.0, 0.0, 0.0]]),
                            molecular_charge=0.0,
                            molecular_multiplicity=1,
                        ),
                    ],
                ),
                HierarchicalFragment(
                    fragment_id="cluster_B",
                    tier=1,
                    parent_path=(1,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="He_B1",
                            tier=2,
                            parent_path=(1, 0),
                            symbols=["He"],
                            geometry=np.array([[0.0, 1.0, 0.0]]),
                            molecular_charge=0.0,
                            molecular_multiplicity=1,
                        ),
                        HierarchicalFragment(
                            fragment_id="He_B2",
                            tier=2,
                            parent_path=(1, 1),
                            symbols=["He"],
                            geometry=np.array([[1.0, 1.0, 0.0]]),
                            molecular_charge=0.0,
                            molecular_multiplicity=1,
                        ),
                    ],
                ),
            ],
        )

    def test_build_flat_molecule_basic(self, simple_2tier_hierarchy):
        """Test building a flat molecule from simple hierarchy."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mol = builder.build_flat_molecule(molecule_name="test_system")

        # Check basic properties
        assert mol.name == "test_system"
        assert len(mol.symbols) == 4  # 4 He atoms
        assert len(mol.fragments) == 4  # 4 elementary fragments
        assert all(sym == "He" for sym in mol.symbols)

        # Check fragments are sequential
        assert mol.fragments[0] == [0]
        assert mol.fragments[1] == [1]
        assert mol.fragments[2] == [2]
        assert mol.fragments[3] == [3]

    def test_flat_molecule_geometry(self, simple_2tier_hierarchy):
        """Test that geometry is correctly flattened."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mol = builder.build_flat_molecule()

        # Expected geometry in order: (0,0,0), (1,0,0), (0,1,0), (1,1,0)
        expected_geom = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])

        np.testing.assert_allclose(mol.geometry, expected_geom)

    def test_flat_molecule_charges(self, simple_2tier_hierarchy):
        """Test that fragment charges are preserved."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mol = builder.build_flat_molecule()

        assert len(mol.fragment_charges) == 4
        assert all(charge == 0.0 for charge in mol.fragment_charges)

    def test_flat_molecule_multiplicities(self, simple_2tier_hierarchy):
        """Test that fragment multiplicities are preserved."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mol = builder.build_flat_molecule()

        assert len(mol.fragment_multiplicities) == 4
        assert all(mult == 1 for mult in mol.fragment_multiplicities)

    def test_fragment_mapping(self, simple_2tier_hierarchy):
        """Test fragment mapping from flat index to hierarchical path."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mapping = builder.get_fragment_mapping()

        assert len(mapping) == 4
        assert mapping[0] == (0, 0)  # cluster_A, fragment 0
        assert mapping[1] == (0, 1)  # cluster_A, fragment 1
        assert mapping[2] == (1, 0)  # cluster_B, fragment 0
        assert mapping[3] == (1, 1)  # cluster_B, fragment 1

    def test_primary_mapping(self, simple_2tier_hierarchy):
        """Test mapping from flat index to primary fragment index."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        mapping = builder.get_primary_mapping()

        assert len(mapping) == 4
        assert mapping[0] == 0  # belongs to primary 0
        assert mapping[1] == 0  # belongs to primary 0
        assert mapping[2] == 1  # belongs to primary 1
        assert mapping[3] == 1  # belongs to primary 1

    def test_fragment_info(self, simple_2tier_hierarchy):
        """Test getting detailed fragment information."""
        builder = HMBEFragmentBuilder(simple_2tier_hierarchy)
        info = builder.get_fragment_info()

        assert len(info) == 4

        # Check first fragment
        frag0 = info[0]
        assert frag0["id"] == "He_A1"
        assert frag0["path"] == (0, 0)
        assert frag0["tier"] == 2
        assert frag0["primary_index"] == 0
        assert frag0["natoms"] == 1
        assert frag0["charge"] == 0.0
        assert frag0["multiplicity"] == 1

    def test_3tier_hierarchy(self):
        """Test building flat molecule from 3-tier hierarchy."""
        hierarchy = FragmentHierarchy(
            tiers=3,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="primary",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="secondary",
                            tier=2,
                            parent_path=(0, 0),
                            sub_fragments=[
                                HierarchicalFragment(
                                    fragment_id="H1",
                                    tier=3,
                                    parent_path=(0, 0, 0),
                                    symbols=["H"],
                                    geometry=np.array([[0, 0, 0]]),
                                ),
                                HierarchicalFragment(
                                    fragment_id="H2",
                                    tier=3,
                                    parent_path=(0, 0, 1),
                                    symbols=["H"],
                                    geometry=np.array([[1, 0, 0]]),
                                ),
                            ],
                        )
                    ],
                )
            ],
        )

        builder = HMBEFragmentBuilder(hierarchy)
        mol = builder.build_flat_molecule()

        assert len(mol.symbols) == 2  # 2 H atoms
        assert len(mol.fragments) == 2  # 2 elementary fragments
        assert mol.fragments[0] == [0]
        assert mol.fragments[1] == [1]

        # Check mapping
        mapping = builder.get_fragment_mapping()
        assert mapping[0] == (0, 0, 0)
        assert mapping[1] == (0, 0, 1)

    def test_multi_atom_fragments(self):
        """Test building molecule with multi-atom elementary fragments."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="primary",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="water",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["O", "H", "H"],
                            geometry=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]),
                        )
                    ],
                )
            ],
        )

        builder = HMBEFragmentBuilder(hierarchy)
        mol = builder.build_flat_molecule()

        assert len(mol.symbols) == 3  # O, H, H
        assert len(mol.fragments) == 1  # 1 elementary fragment
        assert mol.fragments[0] == [0, 1, 2]  # Water molecule
        assert mol.symbols == ["O", "H", "H"]

    def test_charged_fragments(self):
        """Test handling of charged fragments."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="primary",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="cation",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["Na"],
                            geometry=np.array([[0, 0, 0]]),
                            molecular_charge=1.0,
                        ),
                        HierarchicalFragment(
                            fragment_id="anion",
                            tier=2,
                            parent_path=(0, 1),
                            symbols=["Cl"],
                            geometry=np.array([[2, 0, 0]]),
                            molecular_charge=-1.0,
                        ),
                    ],
                )
            ],
        )

        builder = HMBEFragmentBuilder(hierarchy)
        mol = builder.build_flat_molecule()

        assert mol.fragment_charges == [1.0, -1.0]
        assert mol.molecular_charge == 0.0  # Total charge

    def test_different_multiplicities(self):
        """Test handling of different spin multiplicities."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="primary",
                    tier=1,
                    parent_path=(0,),
                    sub_fragments=[
                        HierarchicalFragment(
                            fragment_id="singlet",
                            tier=2,
                            parent_path=(0, 0),
                            symbols=["H", "H"],
                            geometry=np.array([[0, 0, 0], [1, 0, 0]]),
                            molecular_multiplicity=1,
                        ),
                        HierarchicalFragment(
                            fragment_id="doublet",
                            tier=2,
                            parent_path=(0, 1),
                            symbols=["H"],
                            geometry=np.array([[2, 0, 0]]),
                            molecular_multiplicity=2,
                        ),
                    ],
                )
            ],
        )

        builder = HMBEFragmentBuilder(hierarchy)
        mol = builder.build_flat_molecule()

        assert mol.fragment_multiplicities == [1, 2]

    def test_invalid_units(self):
        """Test that invalid units raise an error."""
        hierarchy = FragmentHierarchy(
            tiers=2,
            max_primary_per_nmer=2,
            primary_fragments=[
                HierarchicalFragment(
                    fragment_id="primary",
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

        builder = HMBEFragmentBuilder(hierarchy)
        with pytest.raises(ValueError, match="Units must be"):
            builder.build_flat_molecule(units="invalid")
