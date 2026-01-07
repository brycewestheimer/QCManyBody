"""
Hierarchical fragment representation for HMBE calculations.

This module provides classes for representing and managing hierarchical fragment structures
with arbitrary nesting depth (N-tier hierarchies).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Iterator, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

__all__ = ["HierarchicalFragment", "FragmentHierarchy"]


class HierarchicalFragment:
    """Represents a molecular fragment at any tier in a hierarchical structure.

    Supports arbitrary nesting through recursive `sub_fragments` structure.
    Elementary fragments (leaves) have no sub-fragments and contain atomic data.
    Composite fragments (internal nodes) contain sub-fragments.

    Parameters
    ----------
    fragment_id : str
        Unique identifier for this fragment
    tier : int
        Tier level in hierarchy (1=primary, 2=secondary, 3=tertiary, etc.)
    parent_path : Tuple[int, ...]
        Path from root to this fragment via indices
        Example: (0, 2, 1) means primary[0].secondary[2].tertiary[1]
    symbols : Optional[List[str]]
        Atomic symbols (required for elementary fragments)
    geometry : Optional[np.ndarray]
        Atomic coordinates in shape (natoms, 3) (required for elementary fragments)
    sub_fragments : Optional[List[HierarchicalFragment]]
        Child fragments (for composite fragments)
    molecular_charge : float
        Total molecular charge
    molecular_multiplicity : int
        Spin multiplicity
    """

    def __init__(
        self,
        fragment_id: str,
        tier: int,
        parent_path: Tuple[int, ...],
        symbols: Optional[List[str]] = None,
        geometry: Optional[np.ndarray] = None,
        sub_fragments: Optional[List[HierarchicalFragment]] = None,
        molecular_charge: float = 0.0,
        molecular_multiplicity: int = 1,
    ):
        self.fragment_id = fragment_id
        self.tier = tier
        self.parent_path = parent_path
        self.symbols = symbols
        self.geometry = geometry if geometry is None else np.asarray(geometry)
        self.sub_fragments = sub_fragments if sub_fragments is not None else []
        self.molecular_charge = molecular_charge
        self.molecular_multiplicity = molecular_multiplicity

        # Validate elementary fragments have atomic data
        if self.is_elementary:
            if self.symbols is None or self.geometry is None:
                raise ValueError(
                    f"Elementary fragment '{fragment_id}' must have symbols and geometry. "
                    f"Got symbols={symbols is not None}, geometry={geometry is not None}"
                )
            if len(self.symbols) == 0:
                raise ValueError(f"Elementary fragment '{fragment_id}' must have at least one atom")

            # Validate geometry shape
            geom_array = np.asarray(self.geometry)
            if geom_array.ndim != 2 or geom_array.shape[1] != 3:
                raise ValueError(
                    f"Fragment '{fragment_id}' geometry must be shape (natoms, 3), "
                    f"got {geom_array.shape}"
                )
            if geom_array.shape[0] != len(self.symbols):
                raise ValueError(
                    f"Fragment '{fragment_id}' has {len(self.symbols)} symbols but "
                    f"{geom_array.shape[0]} coordinates"
                )

    @property
    def is_elementary(self) -> bool:
        """True if this is an elementary (leaf) fragment with no sub-fragments."""
        return len(self.sub_fragments) == 0

    @property
    def is_composite(self) -> bool:
        """True if this is a composite (internal node) fragment with sub-fragments."""
        return len(self.sub_fragments) > 0

    @property
    def natoms(self) -> int:
        """Number of atoms in this fragment (elementary only)."""
        if not self.is_elementary:
            raise ValueError(f"Cannot get natoms for composite fragment '{self.fragment_id}'")
        return len(self.symbols)

    def iter_elementary_fragments(self) -> Iterator[HierarchicalFragment]:
        """Recursively iterate over all elementary (leaf) fragments beneath this fragment.

        Yields
        ------
        HierarchicalFragment
            Elementary fragments in depth-first traversal order
        """
        if self.is_elementary:
            yield self
        else:
            for sub_frag in self.sub_fragments:
                yield from sub_frag.iter_elementary_fragments()

    def get_elementary_count(self) -> int:
        """Count total number of elementary fragments beneath this fragment.

        Returns
        -------
        int
            Number of elementary fragments (1 if this is elementary, otherwise sum of children)
        """
        if self.is_elementary:
            return 1
        return sum(sub.get_elementary_count() for sub in self.sub_fragments)

    def get_total_atoms(self) -> int:
        """Count total number of atoms in all elementary fragments beneath this fragment.

        Returns
        -------
        int
            Total atom count
        """
        if self.is_elementary:
            return self.natoms
        return sum(sub.get_total_atoms() for sub in self.sub_fragments)

    def __repr__(self) -> str:
        """String representation for debugging."""
        if self.is_elementary:
            return (
                f"HierarchicalFragment(id='{self.fragment_id}', tier={self.tier}, "
                f"path={self.parent_path}, natoms={self.natoms})"
            )
        else:
            return (
                f"HierarchicalFragment(id='{self.fragment_id}', tier={self.tier}, "
                f"path={self.parent_path}, n_subfragments={len(self.sub_fragments)})"
            )


class FragmentHierarchy:
    """Manages complete N-tier hierarchical fragment structure.

    Provides indexing, querying, and iteration methods for hierarchical fragments.
    Builds flat indices for efficient access to elementary fragments.

    Parameters
    ----------
    tiers : int
        Number of hierarchy tiers (must be >= 2)
    max_primary_per_nmer : int
        Maximum number of primary fragments allowed per N-mer (HMBE constraint)
    primary_fragments : List[HierarchicalFragment]
        Top-level (tier-1) fragments
    """

    def __init__(
        self,
        tiers: int,
        max_primary_per_nmer: int,
        primary_fragments: List[HierarchicalFragment],
    ):
        if tiers < 2:
            raise ValueError(f"HMBE requires at least 2 tiers, got {tiers}")
        if max_primary_per_nmer < 1:
            raise ValueError(f"max_primary_per_nmer must be >= 1, got {max_primary_per_nmer}")
        if len(primary_fragments) == 0:
            raise ValueError("Must have at least one primary fragment")

        self.tiers = tiers
        self.max_primary_per_nmer = max_primary_per_nmer
        self.primary_fragments = primary_fragments

        # Build flat index of all elementary fragments
        self._elementary_fragments: List[HierarchicalFragment] = []
        self._path_to_index: Dict[Tuple[int, ...], int] = {}
        self._build_indices()

        logger.info(
            f"Built FragmentHierarchy: {tiers} tiers, {len(primary_fragments)} primary fragments, "
            f"{self.num_elementary} elementary fragments"
        )

    def _build_indices(self) -> None:
        """Build flat index and path-to-index mapping for elementary fragments."""
        idx = 0
        for prim_frag in self.primary_fragments:
            for elem_frag in prim_frag.iter_elementary_fragments():
                self._elementary_fragments.append(elem_frag)
                self._path_to_index[elem_frag.parent_path] = idx
                idx += 1

    @property
    def num_primary(self) -> int:
        """Number of primary (tier-1) fragments."""
        return len(self.primary_fragments)

    @property
    def num_elementary(self) -> int:
        """Total number of elementary (leaf) fragments."""
        return len(self._elementary_fragments)

    @property
    def total_atoms(self) -> int:
        """Total number of atoms across all elementary fragments."""
        return sum(frag.get_total_atoms() for frag in self.primary_fragments)

    def get_elementary_fragment(self, flat_index: int) -> HierarchicalFragment:
        """Get elementary fragment by flat index.

        Parameters
        ----------
        flat_index : int
            Flat index in range [0, num_elementary)

        Returns
        -------
        HierarchicalFragment
            Elementary fragment at that index
        """
        if not (0 <= flat_index < self.num_elementary):
            raise IndexError(f"Index {flat_index} out of range [0, {self.num_elementary})")
        return self._elementary_fragments[flat_index]

    def get_fragment_by_path(self, path: Tuple[int, ...]) -> HierarchicalFragment:
        """Get elementary fragment by hierarchical path.

        Parameters
        ----------
        path : Tuple[int, ...]
            Hierarchical path (e.g., (0, 1, 2) for primary[0].secondary[1].tertiary[2])

        Returns
        -------
        HierarchicalFragment
            Elementary fragment at that path
        """
        if path not in self._path_to_index:
            raise KeyError(f"Path {path} not found in hierarchy")
        flat_idx = self._path_to_index[path]
        return self._elementary_fragments[flat_idx]

    def get_flat_index(self, path: Tuple[int, ...]) -> int:
        """Get flat index for elementary fragment with given path.

        Parameters
        ----------
        path : Tuple[int, ...]
            Hierarchical path

        Returns
        -------
        int
            Flat index
        """
        if path not in self._path_to_index:
            raise KeyError(f"Path {path} not found in hierarchy")
        return self._path_to_index[path]

    def get_primary_index(self, elementary_path: Tuple[int, ...]) -> int:
        """Get primary fragment index for an elementary fragment.

        Parameters
        ----------
        elementary_path : Tuple[int, ...]
            Path to elementary fragment

        Returns
        -------
        int
            Primary fragment index (first element of path)
        """
        if len(elementary_path) == 0:
            raise ValueError("Path cannot be empty")
        return elementary_path[0]

    def iter_elementary_fragments(self) -> Iterator[HierarchicalFragment]:
        """Iterate over all elementary fragments in flat index order.

        Yields
        ------
        HierarchicalFragment
            Elementary fragments
        """
        yield from self._elementary_fragments

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> FragmentHierarchy:
        """Construct FragmentHierarchy from dictionary (e.g., parsed from JSON/YAML).

        Supports field aliases for backward compatibility:
        - 'tiers' OR 'hmbe_order' (reference format uses string, converted to int)

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with keys:
            - tiers: int (or hmbe_order: int/str)
            - max_primary_per_nmer: int (optional, default 2)
            - fragments: List[Dict] (primary fragments)

        Returns
        -------
        FragmentHierarchy
            Constructed hierarchy

        Example
        -------
        >>> data = {
        ...     "tiers": 2,
        ...     "max_primary_per_nmer": 2,
        ...     "fragments": [
        ...         {
        ...             "id": "primary_1",
        ...             "sub_fragments": [
        ...                 {"id": "elem_1", "symbols": ["He"], "geometry": [[0, 0, 0]]}
        ...             ]
        ...         }
        ...     ]
        ... }
        >>> hierarchy = FragmentHierarchy.from_dict(data)
        """
        # Support both 'tiers' and 'hmbe_order' (reference format)
        tiers = data.get("tiers") or data.get("hmbe_order")
        if tiers is None:
            raise ValueError("Missing required field 'tiers' or 'hmbe_order' in hierarchy data")

        # Convert string to int if needed (reference format uses string)
        if isinstance(tiers, str):
            tiers = int(tiers)

        max_primary = data.get("max_primary_per_nmer", 2)
        fragments_data = data.get("fragments")
        if fragments_data is None:
            raise ValueError("Missing required field 'fragments' in hierarchy data")

        # Parse primary fragments recursively
        primary_fragments = []
        for prim_idx, prim_data in enumerate(fragments_data):
            prim_frag = cls._parse_fragment_recursive(
                frag_data=prim_data, tier=1, parent_path=(prim_idx,), tiers=tiers
            )
            primary_fragments.append(prim_frag)

        return cls(
            tiers=tiers, max_primary_per_nmer=max_primary, primary_fragments=primary_fragments
        )

    @classmethod
    def _parse_fragment_recursive(
        cls, frag_data: Dict[str, Any], tier: int, parent_path: Tuple[int, ...], tiers: int
    ) -> HierarchicalFragment:
        """Recursively parse fragment from dictionary data.

        Supports field aliases for backward compatibility with reference implementation:
        - 'symbols' OR 'elements'
        - 'geometry' OR 'coordinates'
        - 'sub_fragments' OR 'fragments'
        - 'molecular_charge' OR 'charge'
        - 'molecular_multiplicity' OR 'multiplicity'

        Parameters
        ----------
        frag_data : Dict[str, Any]
            Fragment data with keys: id, symbols, geometry, sub_fragments, etc.
        tier : int
            Current tier level
        parent_path : Tuple[int, ...]
            Path to this fragment
        tiers : int
            Total number of tiers in hierarchy

        Returns
        -------
        HierarchicalFragment
            Parsed fragment
        """
        fragment_id = frag_data.get("id")
        if fragment_id is None:
            raise ValueError(f"Fragment at path {parent_path} missing required field 'id'")

        # Get optional fields with alias support (prefer QCElemental-standard names)
        symbols = frag_data.get("symbols") or frag_data.get("elements")
        geometry = frag_data.get("geometry") or frag_data.get("coordinates")
        molecular_charge = frag_data.get("molecular_charge") or frag_data.get("charge", 0.0)
        molecular_multiplicity = frag_data.get("molecular_multiplicity") or frag_data.get("multiplicity", 1)

        # Support both 'sub_fragments' (preferred) and 'fragments' (reference format)
        sub_fragments_data = frag_data.get("sub_fragments") or frag_data.get("fragments", [])

        # Parse sub-fragments if present
        sub_fragments = []
        if sub_fragments_data:
            for sub_idx, sub_data in enumerate(sub_fragments_data):
                sub_path = parent_path + (sub_idx,)
                sub_frag = cls._parse_fragment_recursive(
                    frag_data=sub_data, tier=tier + 1, parent_path=sub_path, tiers=tiers
                )
                sub_fragments.append(sub_frag)

        # Convert geometry to numpy array if present
        geometry_array = None
        if geometry is not None:
            geometry_array = np.array(geometry, dtype=float)

        return HierarchicalFragment(
            fragment_id=fragment_id,
            tier=tier,
            parent_path=parent_path,
            symbols=symbols,
            geometry=geometry_array,
            sub_fragments=sub_fragments,
            molecular_charge=molecular_charge,
            molecular_multiplicity=molecular_multiplicity,
        )

    def __repr__(self) -> str:
        """String representation for debugging."""
        return (
            f"FragmentHierarchy(tiers={self.tiers}, "
            f"num_primary={self.num_primary}, "
            f"num_elementary={self.num_elementary}, "
            f"max_primary_per_nmer={self.max_primary_per_nmer})"
        )
