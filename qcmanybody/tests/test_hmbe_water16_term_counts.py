"""
Verification and printing of HMBE term counts for a 16-fragment 2-tier (4×4) system.

Covers truncation orders: (2,3), (2,4), (3,4)
Prints both MBE and HMBE counts and verifies mode equivalence.
"""

import math
import pytest
from typing import Dict, Set, Tuple
from qcelemental.models import Molecule

from qcmanybody import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.hmbe_enumerate import estimate_mbe_term_count


def create_uniform_2tier_hierarchy(num_fragments: int, frags_per_group: int) -> FragmentHierarchy:
    """Create a uniform 2-tier hierarchy with equal-sized groups."""
    fragment_tiers = {}
    for i in range(num_fragments):
        frag_id = i + 1
        tier1_group = f"G{i // frags_per_group}"
        tier2_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_name)

    return FragmentHierarchy(num_tiers=2, fragment_tiers=fragment_tiers, tier_names=("group", "fragment"))


def extract_fragment_tuples_from_compute_map(compute_map: Dict) -> Set[Tuple[int, ...]]:
    """Extract all unique fragment tuples from the consolidated 'all' compute set."""
    all_frags = set()
    for mc_dict in compute_map.values():
        bsse_all = mc_dict.get("all", {})
        for nbody_set in bsse_all.values():
            for frag_tuple, _ in nbody_set:
                all_frags.add(frag_tuple)
    return all_frags


# Combinatorial helpers

def C(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def expected_hmbe_count_water16_2tier(truncation_orders: Tuple[int, int]) -> int:
    """Closed-form expected HMBE term count for 16 fragments in 4×4 hierarchy (2-tier).

    Assumptions: 4 tier-1 groups, each with 4 fragments. HMBE constraints:
    - Max n-body = T_K
    - Max number of tier-1 groups involved in any tuple = T_1

    Returns total across 1..T_K bodies.
    """
    T_1, T_K = truncation_orders
    n = 16
    gsize = 4
    G = 4  # number of groups

    # 1-body and 2-body are always full sets for T_1 >= 2
    total = C(n, 1) + (C(n, 2) if T_K >= 2 else 0)

    if T_K >= 3:
        if T_1 == 2:
            # Triples spanning at most 2 groups = one-group + two-group contributions
            one_group = G * C(gsize, 3)
            two_groups = C(G, 2) * (C(2 * gsize, 3) - 2 * C(gsize, 3))
            total += one_group + two_groups
        else:  # T_1 >= 3 → every 3-body triple is allowed
            total += C(n, 3)

    if T_K >= 4:
        if T_1 == 2:
            # Quadruples spanning at most 2 groups = one-group + two-group contributions
            one_group = G * C(gsize, 4)
            two_groups = C(G, 2) * (C(2 * gsize, 4) - 2 * C(gsize, 4))
            total += one_group + two_groups
        elif T_1 == 3:
            # All 4-body combos except those using all 4 groups (exactly one per group)
            total += C(n, 4) - (gsize ** 4)
        else:  # T_1 >= 4 → all 4-body are allowed
            total += C(n, 4)

    return total


class TestWater16TwoTierCounts:
    @pytest.mark.parametrize("truncation_orders", [
        (2, 3),
        (2, 4),
        (3, 4),
    ])
    def test_print_and_verify_counts(self, truncation_orders):
        # Build 16-fragment toy molecule (He atoms per fragment to keep minimal)
        num_fragments = 16
        mol = Molecule(
            symbols=["He"] * num_fragments,
            geometry=[[float(i % 8), float(i // 8), 0.0] for i in range(num_fragments)],
            fragments=[[i] for i in range(num_fragments)],
        )

        hierarchy = create_uniform_2tier_hierarchy(num_fragments=16, frags_per_group=4)

        # Model chemistry levels up to T_K
        T_K = truncation_orders[1]
        levels = {i: "hf/sto-3g" for i in range(1, T_K + 1)}

        # Filter mode
        hmbe_spec_filter = HMBESpecification(
            truncation_orders=truncation_orders,
            hierarchy=hierarchy,
            enumeration_mode="filter",
        )
        mbcore_filter = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_filter,
        )

        # Direct mode
        hmbe_spec_direct = HMBESpecification(
            truncation_orders=truncation_orders,
            hierarchy=hierarchy,
            enumeration_mode="direct",
        )
        mbcore_direct = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels=levels,
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec_direct,
        )

        # Extract tuples and counts
        filter_tuples = extract_fragment_tuples_from_compute_map(mbcore_filter.compute_map)
        direct_tuples = extract_fragment_tuples_from_compute_map(mbcore_direct.compute_map)

        hmbe_count = len(filter_tuples)
        mbe_count = estimate_mbe_term_count(num_fragments=16, max_nbody=T_K)
        expected_hmbe = expected_hmbe_count_water16_2tier(truncation_orders)

        # Print to terminal
        print(f"\nWater16 (4×4), trunc={truncation_orders}: HMBE={hmbe_count}, MBE={mbe_count}")

        # Verifications
        assert filter_tuples == direct_tuples, "Filter and direct modes must match"
        assert hmbe_count == expected_hmbe, (
            f"Expected {expected_hmbe} terms for {truncation_orders}, got {hmbe_count}"
        )
        assert hmbe_count <= mbe_count, "HMBE should not exceed full MBE term count"


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 
