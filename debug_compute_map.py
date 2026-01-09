#!/usr/bin/env python3
"""Debug compute_map structure to see what's actually in it."""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
import json

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

print("Compute map structure:")
mc = "hf/sto-3g"
compute_dict = mbc.compute_map[mc]

for bsse_key in compute_dict:
    print(f"\n  BSSE key: '{bsse_key}'")
    if isinstance(compute_dict[bsse_key], dict):
        for nbody_level, terms in compute_dict[bsse_key].items():
            print(f"    nbody_level {nbody_level}: {len(terms)} terms")
            if len(terms) > 0:
                # Show first 3 terms
                for i, (frag, bas) in enumerate(sorted(terms)[:3]):
                    print(f"      {i+1}. frag={frag}, bas={bas}, len(frag)={len(frag)}")
