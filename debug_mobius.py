#!/usr/bin/env python3
"""Debug Möbius inversion to see why 3-body contribution is 0.0."""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.utils import labeler, copy_value

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

# Generate mock energies
mc = "hf/sto-3g"

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

# Now manually perform Möbius inversion with debugging
print("Manual Möbius inversion:")

# Organize energies by fragment tuple
energy_by_frag = {}
for bsse_dict in mbc.compute_map[mc].values():
    if "nocp" in bsse_dict:
        for nbody, terms in bsse_dict["nocp"].items():
            for frag, bas in terms:
                lbl = labeler(mc, frag, bas)
                if lbl in component_results:
                    energy_by_frag[frag] = component_results[lbl]

print(f"\nTotal fragments with energies: {len(energy_by_frag)}")

# Perform Möbius inversion
contrib_by_frag = {}
contribution_by_level = {1: 0.0, 2: 0.0, 3: 0.0}

for n in range(1, 4):
    frags_of_size = sorted([f for f in energy_by_frag if len(f) == n])
    print(f"\nProcessing {len(frags_of_size)} {n}-body terms:")

    for frag in frags_of_size[:3]:  # Show first 3
        subtotal = copy_value(energy_by_frag[frag])
        subs_applied = []

        # Subtract contributions from all proper subsets
        for sub_frag, sub_contrib in contrib_by_frag.items():
            if len(sub_frag) >= n:
                continue
            if set(sub_frag).issubset(frag):
                subtotal -= sub_contrib
                subs_applied.append((sub_frag, sub_contrib))

        contrib_by_frag[frag] = subtotal
        contribution_by_level[n] += subtotal

        print(f"  {frag}: E_total={energy_by_frag[frag]:.6f}, "
              f"subs={len(subs_applied)}, E_contrib={subtotal:.6f}")

    # Sum remaining terms
    for frag in frags_of_size[3:]:
        subtotal = copy_value(energy_by_frag[frag])
        for sub_frag, sub_contrib in contrib_by_frag.items():
            if len(sub_frag) >= n:
                continue
            if set(sub_frag).issubset(frag):
                subtotal -= sub_contrib
        contrib_by_frag[frag] = subtotal
        contribution_by_level[n] += subtotal

    print(f"  Total {n}-body contribution: {contribution_by_level[n]:.6f}")

# Build cumulative body dict
print("\nCumulative energies:")
running = 0.0
for n in range(1, 4):
    running += contribution_by_level[n]
    print(f"  Through {n}-body: {running:.6f} (contrib: {contribution_by_level[n]:.6f})")

print("\n✓ Analysis complete")
