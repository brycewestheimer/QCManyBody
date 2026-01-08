"""
Hierarchical Many-Body Expansion (HMBE) data structures.

This module defines the core data structures for specifying hierarchical
fragment organization and HMBE truncation schemes.
"""

from typing import Any, Dict, FrozenSet, Literal, Optional, Tuple

from pydantic.v1 import BaseModel, Field, validator, root_validator

__all__ = [
    "FragmentHierarchy",
    "SchengenSpecification",
    "HMBESpecification",
]


class FragmentHierarchy(BaseModel):
    """
    Hierarchical organization of fragments into K tiers.

    In HMBE, fragments are organized into a hierarchy of groups at different
    scales. For example, in a protein:
    - Tier 1 (coarsest): Domains
    - Tier 2: Secondary structures (helices, sheets)
    - Tier 3 (finest): Individual residues (elementary fragments)

    Attributes
    ----------
    num_tiers : int
        Number of hierarchical tiers (K). Must be 2 or 3 for current implementation.
        The finest tier (tier-K) contains the elementary fragments.
    fragment_tiers : Dict[int, Tuple[str, ...]]
        Maps 1-indexed fragment ID to tuple of group identifiers at each tier.
        Example for 3-tier: {1: ("Domain1", "Helix1", "Residue1")}
        The tuple length must equal num_tiers for all fragments.
    tier_names : Optional[Tuple[str, ...]]
        Optional human-readable names for each tier (e.g., ("domain", "helix", "residue")).
        Used for documentation and output labeling.

    Examples
    --------
    >>> # 4-fragment system with 2-tier hierarchy (2x2 grouping)
    >>> hierarchy = FragmentHierarchy(
    ...     num_tiers=2,
    ...     fragment_tiers={
    ...         1: ("GroupA", "SubgroupA1"),
    ...         2: ("GroupA", "SubgroupA2"),
    ...         3: ("GroupB", "SubgroupB1"),
    ...         4: ("GroupB", "SubgroupB2"),
    ...     },
    ...     tier_names=("group", "subgroup")
    ... )
    """

    num_tiers: int = Field(..., ge=2, le=3, description="Number of hierarchical tiers (2 or 3)")
    fragment_tiers: Dict[int, Tuple[str, ...]] = Field(
        ..., description="Map from 1-indexed fragment ID to (tier1_group, tier2_group, [tier3_group])"
    )
    tier_names: Optional[Tuple[str, ...]] = Field(
        None, description="Optional names for each tier (e.g., ('domain', 'residue'))"
    )

    class Config:
        extra = "forbid"
        frozen = True  # Immutable after creation for thread safety

    @validator("fragment_tiers")
    @classmethod
    def validate_fragment_tiers(cls, v, values):
        """Ensure all fragments have correct tier depth."""
        num_tiers = values.get("num_tiers")
        if num_tiers is None:
            return v

        for frag_id, tier_path in v.items():
            if len(tier_path) != num_tiers:
                raise ValueError(
                    f"Fragment {frag_id} has {len(tier_path)} tier(s), expected {num_tiers}. "
                    f"All fragments must have exactly {num_tiers} tier identifiers."
                )
        return v

    @validator("tier_names")
    @classmethod
    def validate_tier_names(cls, v, values):
        """Ensure tier_names length matches num_tiers if provided."""
        if v is not None:
            num_tiers = values.get("num_tiers")
            if num_tiers is not None and len(v) != num_tiers:
                raise ValueError(f"tier_names has {len(v)} entries, expected {num_tiers} to match num_tiers")
        return v

    def get_tier_groups(self, tier_idx: int) -> Dict[str, FrozenSet[int]]:
        """
        Get mapping from tier group ID to set of elementary fragments.

        Parameters
        ----------
        tier_idx : int
            0-indexed tier (0 = tier-1, 1 = tier-2, etc.)

        Returns
        -------
        Dict[str, FrozenSet[int]]
            Map from group_id to frozenset of 1-indexed fragment IDs belonging to that group

        Raises
        ------
        ValueError
            If tier_idx is out of range

        Examples
        --------
        >>> groups = hierarchy.get_tier_groups(0)  # Get tier-1 groups
        >>> groups["GroupA"]  # frozenset({1, 2})
        """
        if tier_idx < 0 or tier_idx >= self.num_tiers:
            raise ValueError(f"tier_idx {tier_idx} out of range [0, {self.num_tiers - 1}]")

        groups = {}
        for frag_id, tier_path in self.fragment_tiers.items():
            group_id = tier_path[tier_idx]
            if group_id not in groups:
                groups[group_id] = set()
            groups[group_id].add(frag_id)

        # Convert to frozenset for immutability
        return {k: frozenset(v) for k, v in groups.items()}

    def count_distinct_groups_at_tier(self, fragments: Tuple[int, ...], tier_idx: int) -> int:
        """
        Count distinct tier-t groups containing any of the given fragments.

        This is the core filtering operation for HMBE. An n-body term is kept
        only if the count of distinct groups at each tier is within the truncation
        order for that tier.

        Parameters
        ----------
        fragments : Tuple[int, ...]
            Tuple of 1-indexed fragment IDs
        tier_idx : int
            0-indexed tier (0 = tier-1, 1 = tier-2, etc.)

        Returns
        -------
        int
            Number of distinct groups at the specified tier

        Examples
        --------
        >>> # Fragments 1 and 2 both belong to "GroupA" at tier-1
        >>> hierarchy.count_distinct_groups_at_tier((1, 2), 0)  # Returns 1
        >>> # Fragments 1 and 3 belong to different tier-1 groups
        >>> hierarchy.count_distinct_groups_at_tier((1, 3), 0)  # Returns 2
        """
        distinct_groups = {self.fragment_tiers[frag][tier_idx] for frag in fragments}
        return len(distinct_groups)

    def get_fragments_in_group(self, group_id: str, tier_idx: int) -> FrozenSet[int]:
        """
        Get all fragments belonging to a specific group at a given tier.

        Parameters
        ----------
        group_id : str
            Group identifier at the specified tier
        tier_idx : int
            0-indexed tier (0 = tier-1, 1 = tier-2, etc.)

        Returns
        -------
        FrozenSet[int]
            Frozenset of 1-indexed fragment IDs in the specified group

        Examples
        --------
        >>> # Get all fragments in "GroupA" at tier-1
        >>> frags = hierarchy.get_fragments_in_group("GroupA", 0)
        >>> frags  # frozenset({1, 2, 3, 4})
        """
        tier_groups = self.get_tier_groups(tier_idx)
        return tier_groups.get(group_id, frozenset())


class SchengenSpecification(BaseModel):
    """
    Specification for Schengen term selection in HMBE.

    Schengen terms are n-body interactions excluded by the base HMBE truncation
    but selected for inclusion based on spatial proximity. This can improve
    accuracy at interfaces between hierarchical groups.

    Attributes
    ----------
    enabled : bool
        Whether to add Schengen terms to the base HMBE calculation
    selection_fraction : float
        Fraction of candidate terms to include (0.0 to 1.0).
        E.g., 0.1 means top 10% by proximity metric.
    distance_metric : str
        Metric for ranking fragment proximity. Options:
        - "R2": Sum of squared pairwise distances (recommended, default)
        - "R": Sum of pairwise distances
        - "R_inv": Sum of inverse distances (negative for sorting)
        - "R3_inv": Sum of inverse cube distances
        - "fmo": FMO-style scaling (currently falls back to R2)

    Examples
    --------
    >>> schengen = SchengenSpecification(
    ...     enabled=True,
    ...     selection_fraction=0.1,  # Top 10%
    ...     distance_metric="R2"
    ... )
    """

    enabled: bool = Field(False, description="Whether to add Schengen terms")
    selection_fraction: float = Field(
        0.1, ge=0.0, le=1.0, description="Fraction of candidate terms to include (e.g., 0.1 = top 10%)"
    )
    distance_metric: Literal["R2", "R", "R_inv", "R3_inv", "fmo"] = Field(
        "R2", description="Distance metric for ranking terms. R2 recommended."
    )

    class Config:
        frozen = True


class HMBESpecification(BaseModel):
    """
    Complete specification for Hierarchical Many-Body Expansion calculation.

    HMBE reduces computational cost by applying different truncation orders at
    different hierarchical scales. The truncation orders must be monotonically
    non-decreasing from coarse to fine tiers.

    Attributes
    ----------
    truncation_orders : Tuple[int, ...]
        Truncation orders (T_1, T_2, [T_3]) from coarsest to finest tier.
        Must satisfy T_1 <= T_2 <= [T_3].
        Example: (2, 3) for (2,3)-HMBE means:
        - At tier-1 (coarse): maximum 2-body interactions between groups
        - At tier-2 (fine): maximum 3-body interactions between fragments
    hierarchy : FragmentHierarchy
        Definition of the hierarchical fragment organization
    schengen : Optional[SchengenSpecification]
        Optional Schengen term selection for improved accuracy at interfaces
    enumeration_mode : str
        Method for generating HMBE terms. Options:
        - "filter": Generate all MBE terms, then filter to HMBE subset (default)
        - "direct": Generate only HMBE terms directly (top-down, faster for large systems)
        - "auto": Automatically choose based on system size (<30 fragments: filter, ≥30: direct)

    Examples
    --------
    >>> # (2,3)-HMBE with Schengen terms
    >>> hmbe_spec = HMBESpecification(
    ...     truncation_orders=(2, 3),
    ...     hierarchy=hierarchy,
    ...     schengen=SchengenSpecification(enabled=True, selection_fraction=0.1)
    ... )
    >>> hmbe_spec.max_nbody  # Returns 3
    >>>
    >>> # Use direct enumeration for large system
    >>> hmbe_spec_large = HMBESpecification(
    ...     truncation_orders=(2, 4),
    ...     hierarchy=hierarchy,
    ...     enumeration_mode="direct"
    ... )
    """

    truncation_orders: Tuple[int, ...] = Field(
        ..., description="Truncation orders (T_1, T_2, [T_3]) from coarse to fine"
    )
    hierarchy: FragmentHierarchy = Field(..., description="Fragment hierarchy definition")
    schengen: Optional[SchengenSpecification] = Field(
        None, description="Optional Schengen term specification"
    )
    enumeration_mode: Literal["filter", "direct", "auto"] = Field(
        "auto",
        description="Term enumeration method: 'filter' (generate all MBE then filter), "
        "'direct' (generate only HMBE terms), or 'auto' (choose based on system size)"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @root_validator
    @classmethod
    def validate_spec(cls, values):
        """Validate consistency between truncation_orders and hierarchy."""
        truncation_orders = values.get("truncation_orders")
        hierarchy = values.get("hierarchy")

        if truncation_orders and hierarchy:
            # Check that number of truncation orders matches num_tiers
            if len(truncation_orders) != hierarchy.num_tiers:
                raise ValueError(
                    f"Number of truncation orders ({len(truncation_orders)}) must match "
                    f"number of tiers ({hierarchy.num_tiers})"
                )

            # Check monotonicity: T_1 <= T_2 <= [T_3]
            for i in range(len(truncation_orders) - 1):
                if truncation_orders[i] > truncation_orders[i + 1]:
                    raise ValueError(
                        f"Truncation orders must be non-decreasing (T_1 <= T_2 <= ... <= T_K). "
                        f"Got T_{i + 1}={truncation_orders[i]} > T_{i + 2}={truncation_orders[i + 1]}"
                    )

        return values

    @property
    def max_nbody(self) -> int:
        """Maximum n-body order (T_K, the finest-tier truncation)."""
        return self.truncation_orders[-1]

    @property
    def num_tiers(self) -> int:
        """Number of hierarchical tiers."""
        return self.hierarchy.num_tiers

    def is_standard_mbe(self) -> bool:
        """
        Check if this specification reduces to standard MBE.

        Returns True if all truncation orders are equal, meaning no hierarchical
        truncation is applied (equivalent to standard MBE-n).
        """
        return len(set(self.truncation_orders)) == 1
