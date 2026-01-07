"""
Fragment builder for converting hierarchical HMBE structure to flat QCElemental Molecule.

This module handles the conversion from hierarchical fragment representation
to flat elementary fragments suitable for standard ManyBodyComputer.
"""

from __future__ import annotations

import logging
from typing import Dict, Tuple

import numpy as np
from qcelemental.models import Molecule

from .hierarchy import FragmentHierarchy, HierarchicalFragment

logger = logging.getLogger(__name__)

__all__ = ["HMBEFragmentBuilder"]


class HMBEFragmentBuilder:
    """Builds flat fragment list and QCElemental Molecule from hierarchical structure.

    The builder flattens the N-tier hierarchy into a list of elementary fragments
    compatible with standard QCManyBody workflows. It maintains mapping information
    to trace elementary fragments back to their hierarchical origins.

    Parameters
    ----------
    hierarchy : FragmentHierarchy
        The hierarchical fragment structure to convert
    """

    def __init__(self, hierarchy: FragmentHierarchy):
        self.hierarchy = hierarchy

    def build_flat_molecule(
        self,
        molecule_name: str = "hmbe_system",
        units: str = "bohr",
    ) -> Molecule:
        """Construct QCElemental Molecule with flat elementary fragment structure.

        Converts the hierarchical fragment structure into a single Molecule object
        where each fragment is an elementary fragment from the hierarchy.

        Parameters
        ----------
        molecule_name : str, optional
            Name for the molecule (default: "hmbe_system")
        units : str, optional
            Units for geometry, either "bohr" or "angstrom" (default: "bohr")

        Returns
        -------
        Molecule
            QCElemental Molecule where fragments are elementary fragments
            in hierarchy traversal order

        Raises
        ------
        ValueError
            If units are invalid or geometry data is malformed
        """
        if units not in ["bohr", "angstrom"]:
            raise ValueError(f"Units must be 'bohr' or 'angstrom', got '{units}'")

        all_symbols = []
        all_geometry_flat = []
        fragments = []  # List of atom index lists
        fragment_charges = []
        fragment_multiplicities = []

        atom_offset = 0

        # Iterate through elementary fragments in flat index order
        for elem_frag in self.hierarchy.iter_elementary_fragments():
            natoms = elem_frag.natoms

            # Append atomic symbols
            all_symbols.extend(elem_frag.symbols)

            # Append geometry (flatten to 1D)
            geom_array = np.asarray(elem_frag.geometry)
            all_geometry_flat.extend(geom_array.flatten().tolist())

            # Define fragment as list of atom indices
            fragment_indices = list(range(atom_offset, atom_offset + natoms))
            fragments.append(fragment_indices)

            # Fragment properties
            fragment_charges.append(elem_frag.molecular_charge)
            fragment_multiplicities.append(elem_frag.molecular_multiplicity)

            atom_offset += natoms

        # Reshape geometry to (natoms, 3)
        geometry_array = np.array(all_geometry_flat).reshape(-1, 3)

        # Create QCElemental Molecule
        mol = Molecule(
            symbols=all_symbols,
            geometry=geometry_array.tolist(),
            fragments=fragments,
            fragment_charges=fragment_charges,
            fragment_multiplicities=fragment_multiplicities,
            name=molecule_name,
            molecular_charge=sum(fragment_charges),
            molecular_multiplicity=1,  # TODO: compute properly for full system
        )

        logger.info(
            f"Built flat molecule '{molecule_name}': "
            f"{len(all_symbols)} atoms, {len(fragments)} fragments"
        )

        return mol

    def get_fragment_mapping(self) -> Dict[int, Tuple[int, ...]]:
        """Map flat fragment index to hierarchical path.

        Returns
        -------
        Dict[int, Tuple[int, ...]]
            Mapping from flat elementary fragment index to hierarchical path.
            Example: {0: (0, 0), 1: (0, 1), 2: (1, 0)}
            means fragment 0 is primary[0].secondary[0],
                  fragment 1 is primary[0].secondary[1],
                  fragment 2 is primary[1].secondary[0]
        """
        mapping = {}
        for flat_idx, elem_frag in enumerate(self.hierarchy.iter_elementary_fragments()):
            mapping[flat_idx] = elem_frag.parent_path
        return mapping

    def get_primary_mapping(self) -> Dict[int, int]:
        """Map flat fragment index to primary fragment index.

        Returns
        -------
        Dict[int, int]
            Mapping from flat elementary fragment index to primary fragment index.
            Example: {0: 0, 1: 0, 2: 1, 3: 1}
            means fragments 0-1 belong to primary fragment 0,
                  fragments 2-3 belong to primary fragment 1
        """
        mapping = {}
        for flat_idx, elem_frag in enumerate(self.hierarchy.iter_elementary_fragments()):
            primary_idx = self.hierarchy.get_primary_index(elem_frag.parent_path)
            mapping[flat_idx] = primary_idx
        return mapping

    def get_fragment_info(self) -> Dict[int, Dict[str, any]]:
        """Get detailed information about each elementary fragment.

        Returns
        -------
        Dict[int, Dict[str, any]]
            Mapping from flat index to fragment information dictionary with keys:
            - 'id': fragment identifier
            - 'path': hierarchical path
            - 'tier': tier level
            - 'primary_index': primary fragment index
            - 'natoms': number of atoms
            - 'charge': molecular charge
            - 'multiplicity': spin multiplicity
        """
        info = {}
        for flat_idx, elem_frag in enumerate(self.hierarchy.iter_elementary_fragments()):
            info[flat_idx] = {
                "id": elem_frag.fragment_id,
                "path": elem_frag.parent_path,
                "tier": elem_frag.tier,
                "primary_index": self.hierarchy.get_primary_index(elem_frag.parent_path),
                "natoms": elem_frag.natoms,
                "charge": elem_frag.molecular_charge,
                "multiplicity": elem_frag.molecular_multiplicity,
            }
        return info
