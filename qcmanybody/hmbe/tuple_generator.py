"""
N-mer tuple generation for HMBE calculations.

This module implements the core HMBE algorithm for generating valid N-mer tuples
that respect hierarchical constraints. The primary constraint is that each N-mer
can contain elementary fragments from at most max_primary_per_nmer different
primary fragments.
"""

from __future__ import annotations

import itertools
import logging
from typing import Dict, List, Set, Tuple

from .hierarchy import FragmentHierarchy

logger = logging.getLogger(__name__)

__all__ = ["HMBETupleGenerator"]


class HMBETupleGenerator:
    """Generates valid N-mer tuples respecting hierarchical constraints.

    The core HMBE constraint is that each N-mer can contain elementary fragments
    from at most `max_primary_per_nmer` different primary fragments. This dramatically
    reduces the number of N-mers compared to standard MBE.

    Parameters
    ----------
    hierarchy : FragmentHierarchy
        The hierarchical fragment structure

    Examples
    --------
    >>> hierarchy = FragmentHierarchy.from_dict({...})
    >>> generator = HMBETupleGenerator(hierarchy)
    >>> tuples = generator.generate_nmer_tuples(max_nbody=3)
    >>> print(tuples[1])  # Monomers
    [(0,), (1,), (2,), (3,)]
    >>> print(tuples[2])  # Dimers (restricted by hierarchy)
    [(0, 1), (0, 2), (1, 2), (2, 3)]
    """

    def __init__(self, hierarchy: FragmentHierarchy):
        self.hierarchy = hierarchy
        self.max_primary = hierarchy.max_primary_per_nmer

        # Pre-compute primary index for each elementary fragment
        self._elementary_to_primary: Dict[int, int] = {}
        for flat_idx in range(self.hierarchy.num_elementary):
            elem_frag = self.hierarchy.get_elementary_fragment(flat_idx)
            primary_idx = self.hierarchy.get_primary_index(elem_frag.parent_path)
            self._elementary_to_primary[flat_idx] = primary_idx

        logger.info(
            f"Initialized HMBETupleGenerator: "
            f"{self.hierarchy.num_elementary} elementary fragments, "
            f"max_primary_per_nmer={self.max_primary}"
        )

    def generate_nmer_tuples(
        self,
        max_nbody: int,
    ) -> Dict[int, List[Tuple[int, ...]]]:
        """Generate all valid N-mer tuples for each n-body level.

        Parameters
        ----------
        max_nbody : int
            Maximum n-body level to generate (e.g., 2 for dimers, 3 for trimers)

        Returns
        -------
        Dict[int, List[Tuple[int, ...]]]
            Mapping from n-body level to list of elementary fragment index tuples.
            Example: {1: [(0,), (1,), (2,)], 2: [(0,1), (0,2), (1,2)]}

        Notes
        -----
        Tuples are sorted and deduplicated to avoid redundant calculations.
        """
        if max_nbody < 1:
            raise ValueError(f"max_nbody must be >= 1, got {max_nbody}")

        if max_nbody > self.hierarchy.num_elementary:
            logger.warning(
                f"max_nbody={max_nbody} exceeds num_elementary={self.hierarchy.num_elementary}, "
                f"capping at {self.hierarchy.num_elementary}"
            )
            max_nbody = self.hierarchy.num_elementary

        tuples_by_level = {}

        for nbody in range(1, max_nbody + 1):
            tuples = self._generate_nbody_tuples(nbody)
            tuples_by_level[nbody] = tuples
            logger.info(f"Generated {len(tuples)} valid {nbody}-body tuples")

        return tuples_by_level

    def _generate_nbody_tuples(self, nbody: int) -> List[Tuple[int, ...]]:
        """Generate tuples for a specific n-body level.

        Parameters
        ----------
        nbody : int
            N-body level (1=monomer, 2=dimer, etc.)

        Returns
        -------
        List[Tuple[int, ...]]
            List of valid tuples as sorted integer tuples
        """
        valid_tuples = []
        num_elem = self.hierarchy.num_elementary

        # Generate all combinations of nbody elementary fragments
        for combo in itertools.combinations(range(num_elem), nbody):
            if self._is_valid_hmbe_tuple(combo):
                valid_tuples.append(combo)

        return valid_tuples

    def _is_valid_hmbe_tuple(self, elementary_indices: Tuple[int, ...]) -> bool:
        """Check if tuple satisfies HMBE constraints.

        The primary constraint is that elementary fragments must come from at most
        `max_primary_per_nmer` different primary fragments.

        Parameters
        ----------
        elementary_indices : Tuple[int, ...]
            Tuple of elementary fragment flat indices

        Returns
        -------
        bool
            True if tuple is valid for HMBE, False otherwise
        """
        # Get primary indices for each elementary fragment
        primary_indices = set(
            self._elementary_to_primary[idx] for idx in elementary_indices
        )

        # Check primary constraint
        if len(primary_indices) > self.max_primary:
            return False

        # Additional tier-specific constraints can be added here
        # For example, could add max secondary fragments per primary
        # For now, only primary constraint is enforced

        return True

    def get_tuple_info(
        self, elementary_tuple: Tuple[int, ...]
    ) -> Dict[str, any]:
        """Get detailed information about an N-mer tuple.

        Parameters
        ----------
        elementary_tuple : Tuple[int, ...]
            Tuple of elementary fragment flat indices

        Returns
        -------
        Dict[str, any]
            Dictionary with keys:
            - 'nbody': n-body level
            - 'elementary_indices': flat indices of elementary fragments
            - 'primary_indices': list of primary indices involved
            - 'num_primary': number of different primary fragments
            - 'is_valid': whether tuple satisfies HMBE constraints
            - 'fragment_ids': list of fragment identifiers
        """
        primary_indices = [
            self._elementary_to_primary[idx] for idx in elementary_tuple
        ]
        unique_primary = set(primary_indices)

        fragment_ids = [
            self.hierarchy.get_elementary_fragment(idx).fragment_id
            for idx in elementary_tuple
        ]

        return {
            "nbody": len(elementary_tuple),
            "elementary_indices": list(elementary_tuple),
            "primary_indices": primary_indices,
            "num_primary": len(unique_primary),
            "is_valid": self._is_valid_hmbe_tuple(elementary_tuple),
            "fragment_ids": fragment_ids,
        }

    def count_tuples_by_level(self, max_nbody: int) -> Dict[int, int]:
        """Count number of valid tuples at each n-body level without generating them.

        Useful for estimating computational cost before running calculations.

        Parameters
        ----------
        max_nbody : int
            Maximum n-body level

        Returns
        -------
        Dict[int, int]
            Mapping from n-body level to tuple count
        """
        counts = {}
        for nbody in range(1, max_nbody + 1):
            tuples = self._generate_nbody_tuples(nbody)
            counts[nbody] = len(tuples)
        return counts

    def generate_tuples_by_primary(
        self, max_nbody: int
    ) -> Dict[int, Dict[Tuple[int, ...], List[Tuple[int, ...]]]]:
        """Generate tuples organized by primary fragment combinations.

        This is useful for analysis and understanding the hierarchical structure
        of the HMBE calculations.

        Parameters
        ----------
        max_nbody : int
            Maximum n-body level

        Returns
        -------
        Dict[int, Dict[Tuple[int, ...], List[Tuple[int, ...]]]]
            Nested dictionary: nbody -> primary_combo -> list of elementary tuples
            Example: {2: {(0,): [(0,1)], (1,): [(2,3)], (0,1): [(0,2), (1,3)]}}
        """
        tuples_by_level = self.generate_nmer_tuples(max_nbody)
        organized = {}

        for nbody, tuples in tuples_by_level.items():
            organized[nbody] = {}

            for elem_tuple in tuples:
                # Get primary combination for this tuple
                primary_combo = tuple(
                    sorted(
                        set(
                            self._elementary_to_primary[idx]
                            for idx in elem_tuple
                        )
                    )
                )

                if primary_combo not in organized[nbody]:
                    organized[nbody][primary_combo] = []

                organized[nbody][primary_combo].append(elem_tuple)

        return organized
