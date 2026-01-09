#!/usr/bin/env python3
"""Debug why 3-body contribution is 0.0 for (2,3)-HMBE - improved version."""

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

# Extract ALL fragment tuples and organize by size
mc = "hf/sto-3g"
all_frags = set()
for bsse_dict in mbc.compute_map[mc].values():
    for nb_level in bsse_dict:
        for frag, bas in bsse_dict[nb_level]:
            all_frags.add(frag)

# Organize by fragment count
by_size = {i: [] for i in range(1, 4)}
for frag in all_frags:
    n = len(frag)
    if n <= 3:
        by_size[n].append(frag)

print("Fragment terms by size:")
for n in range(1, 4):
    print(f"  {n}-body: {len(by_size[n])} terms")
    if n == 3 and len(by_size[n]) > 0:
        print(f"    Examples: {sorted(by_size[n])[:5]}")

# Now generate mock energies and check assembly
from qcmanybody.utils import labeler

def get_mock_energy(frag_tuple):
    """Energy model with weak attractive interactions."""
    n = len(frag_tuple)
    return -1.0 * n - 0.02 * n * (n - 1) / 2

component_results = {}
for bsse_dict in mbc.compute_map[mc].values():
    for terms in bsse_dict.values():
        for frag, bas in terms:
            lbl = labeler(mc, frag, bas)
            component_results[lbl] = get_mock_energy(frag)

print(f"\nGenerated {len(component_results)} energy values")

# Assemble
results = mbc._assemble_nbody_components("energy", component_results)
body_dict = results["energy_body_dict"]

print("\nEnergy body dict:")
for nb in sorted(body_dict.keys()):
    print(f"  {nb}-body total: {body_dict[nb]:.6f}")

# Compute contributions
if len(body_dict) >= 2:
    contrib_2b = body_dict[2] - body_dict[1]
    print(f"\n2-body contribution: {contrib_2b:.6f}")

if len(body_dict) >= 3:
    contrib_3b = body_dict[3] - body_dict[2]
    print(f"3-body contribution: {contrib_3b:.6f}")

    if len(by_size[3]) == 0:
        print("\n⚠ Note: 0 three-body terms in HMBE filter!")
        print("This explains why 3-body contribution is 0.0")
    else:
        print(f"\n✓ {len(by_size[3])} three-body terms present")
