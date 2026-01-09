#!/usr/bin/env python3
"""Debug why 3-body contribution is 0.0 for (2,3)-HMBE."""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification

mol = Molecule(
    symbols=["H"] * 6,
    geometry=[[float(i), 0.0, 0.0] for i in range(6)],
    fragments=[[i] for i in range(6)],
)

hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers={
        1: ("G0", "F0"), 2: ("G0", "F1"),
        3: ("G1", "F2"), 4: ("G1", "F3"),
        5: ("G2", "F4"), 6: ("G2", "F5"),
    },
    tier_names=("group", "fragment")
)

hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3),
    hierarchy=hierarchy,
    enumeration_mode="filter"
)

mbc = ManyBodyCore(
    molecule=mol,
    bsse_type=[BsseEnum.nocp],
    levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
    return_total_data=True,
    supersystem_ie_only=False,
    embedding_charges={},
    hmbe_spec=hmbe_spec
)

# Check what 3-body terms are included
mc = "hf/sto-3g"
three_body_terms = set()
for bsse_dict in mbc.compute_map[mc].values():
    if 3 in bsse_dict:
        for frag, bas in bsse_dict[3]:
            three_body_terms.add(frag)

print(f"Number of 3-body terms in HMBE: {len(three_body_terms)}")
print(f"3-body terms: {sorted(three_body_terms)[:10]}")

# For (2,3)-HMBE with 6 fragments in 3 groups:
# - T_1 = 2: max 2 groups
# - T_2 = 3: max 3-body

# Check which 3-body terms should be included
from qcmanybody.hmbe_filter import passes_hmbe_filter

# Example 3-body terms:
test_terms = [
    (1, 2, 3),  # 2 groups (G0, G1) - should PASS
    (1, 2, 5),  # 2 groups (G0, G2) - should PASS
    (1, 3, 5),  # 3 groups (G0, G1, G2) - should FAIL
    (3, 4, 5),  # 2 groups (G1, G2) - should PASS
]

print("\nChecking example terms:")
for term in test_terms:
    passes = passes_hmbe_filter(term, hmbe_spec)
    groups = {hierarchy.fragment_tiers[f][0] for f in term}
    print(f"  {term}: {len(groups)} groups {groups} -> {'PASS' if passes else 'FAIL'}")

print(f"\n✓ Found {len(three_body_terms)} 3-body terms")
print("Note: If there are 3-body terms but contribution is 0.0, check if max_nbody in levels matches")
