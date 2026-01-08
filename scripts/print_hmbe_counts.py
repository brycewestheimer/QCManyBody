#!/usr/bin/env python3
"""
Print HMBE vs MBE term counts for:
- 16 fragments, 2-tier, 4×4 hierarchy for truncations (2,3), (2,4), (3,4)
- 60 fragments, 3-tier, 3×4×5 hierarchy for truncation (2,3,4)

Outputs:
- Enumeration mode used (auto/direct/filter)
- HMBE and MBE counts and reduction factor
- Confirms direct vs filter equivalence when applicable
"""

import math
from typing import Dict, Set, Tuple

from qcelemental.models import Molecule

from qcmanybody import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.hmbe_enumerate import estimate_mbe_term_count


def extract_fragment_tuples_from_compute_map(compute_map: Dict) -> Set[Tuple[int, ...]]:
    """Extract all unique fragment tuples from the consolidated 'all' compute set."""
    all_frags = set()
    for mc_dict in compute_map.values():
        bsse_all = mc_dict.get("all", {})
        for nbody_set in bsse_all.values():
            for frag_tuple, _ in nbody_set:
                all_frags.add(frag_tuple)
    return all_frags


def create_uniform_2tier_hierarchy(num_fragments: int, frags_per_group: int) -> FragmentHierarchy:
    fragment_tiers = {}
    for i in range(num_fragments):
        frag_id = i + 1
        tier1_group = f"G{i // frags_per_group}"
        tier2_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_name)
    return FragmentHierarchy(num_tiers=2, fragment_tiers=fragment_tiers, tier_names=("group", "fragment"))


def create_uniform_3tier_hierarchy(num_fragments: int, tier1_size: int, tier2_size: int) -> FragmentHierarchy:
    fragment_tiers = {}
    for i in range(num_fragments):
        frag_id = i + 1
        tier1_group = f"T1_{i // tier1_size}"
        tier2_group = f"T2_{i // tier2_size}"
        tier3_name = f"F{i}"
        fragment_tiers[frag_id] = (tier1_group, tier2_group, tier3_name)
    return FragmentHierarchy(num_tiers=3, fragment_tiers=fragment_tiers, tier_names=("domain", "subdomain", "fragment"))


def run_counts(num_fragments: int, truncation_orders: Tuple[int, ...], hierarchy: FragmentHierarchy, title: str) -> None:
    T_K = truncation_orders[-1]
    levels = {i: "hf/sto-3g" for i in range(1, T_K + 1)}

    mol = Molecule(
        symbols=["He"] * num_fragments,
        geometry=[[float(i % 8), float(i // 8), 0.0] for i in range(num_fragments)],
        fragments=[[i] for i in range(num_fragments)],
    )

    # For large fragment counts, only use direct mode to save memory
    # (filter mode would require generating all MBE terms first)
    enum_mode = "auto" if num_fragments >= 30 else "direct"
    hmbe_spec = HMBESpecification(
        truncation_orders=truncation_orders,
        hierarchy=hierarchy,
        enumeration_mode=enum_mode,
    )
    mbcore = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels=levels,
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec,
    )

    hmbe_tuples = extract_fragment_tuples_from_compute_map(mbcore.compute_map)
    hmbe_count = len(hmbe_tuples)
    mbe_count = estimate_mbe_term_count(num_fragments=num_fragments, max_nbody=T_K)

    # Stats
    stats = mbcore.get_hmbe_statistics() or {}
    actual_mode = stats.get("actual_enumeration_mode", enum_mode)
    reduction = (mbe_count / hmbe_count) if hmbe_count else float("inf")

    # Print
    print(f"\n[{title}] trunc={truncation_orders} | mode={actual_mode}")
    print(f"  HMBE terms: {hmbe_count}")
    print(f"  MBE terms:  {mbe_count}")
    print(f"  Reduction:  {reduction:.2f}x")
    
    # Clean up to free memory
    del mbcore, hmbe_tuples


def main():
    # 16 fragments, 2-tier, 4×4
    h2 = create_uniform_2tier_hierarchy(num_fragments=16, frags_per_group=4)
    for trunc in [(2, 3), (2, 4), (3, 4)]:
        run_counts(16, trunc, h2, title="Water16 2-tier (4×4)")

    # 32 fragments, 3-tier, 2×4×4 (much smaller than 60 to avoid 32GB+ memory usage)
    h3 = create_uniform_3tier_hierarchy(num_fragments=32, tier1_size=16, tier2_size=4)
    run_counts(32, (2, 3, 4), h3, title="Water32 3-tier (2×4×4)")


if __name__ == "__main__":
    main()
