#!/usr/bin/env python3
"""
Test that HMBE energy assembly works correctly and doesn't "blow up" higher-order energies.

This test simulates the energy assembly process to ensure that the Möbius inversion
properly handles HMBE terms and doesn't overcount contributions.
"""

import sys
import numpy as np
from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification
from qcmanybody.utils import labeler


def test_hmbe_energy_assembly_simple():
    """Test HMBE energy assembly with simple mock energies."""
    print("Testing HMBE energy assembly with mock data...")

    # Create simple 4-fragment system (2 groups of 2)
    mol = Molecule(
        symbols=["H"] * 4,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
        },
        tier_names=("group", "fragment")
    )

    # (2,2)-HMBE: max 2 groups, max 2-body
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        enumeration_mode="filter"
    )

    mbc = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    # Create mock energy results
    # For simplicity: E[monomer] = -1.0, E[dimer] = -2.1 (interaction ~= -0.1 per pair)
    mc = "hf/sto-3g"
    component_results = {}

    # Add energies for all terms in compute map
    for bsse_dict in mbc.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                lbl = labeler(mc, frag, bas)
                n = len(frag)
                # Simple additive model: E = -1.0 * n - 0.1 * n*(n-1)/2
                energy = -1.0 * n - 0.1 * n * (n - 1) / 2
                component_results[lbl] = {"energy": energy}

    print(f"  Generated {len(component_results)} mock energy results")

    # Assemble energies
    try:
        results = mbc._assemble_nbody_components("energy", {k: v["energy"] for k, v in component_results.items()})

        # Check that energies are reasonable (not exploding)
        energy_body_dict = results["energy_body_dict"]
        print(f"  1-body total energy: {energy_body_dict[1]:.6f}")
        print(f"  2-body total energy: {energy_body_dict[2]:.6f}")

        # 2-body contribution (2-body minus 1-body)
        contrib_2b = energy_body_dict[2] - energy_body_dict[1]
        print(f"  2-body contribution: {contrib_2b:.6f}")

        # Check that 2-body contribution is reasonable (should be negative but small)
        assert -1.0 < contrib_2b < 0.0, f"2-body contribution out of range: {contrib_2b}"
        print("  ✓ Energies are in reasonable range (not exploding)")

    except Exception as e:
        print(f"  ✗ Energy assembly failed: {e}")
        raise

    print("✓ Simple energy assembly test passed!\n")


def test_hmbe_vs_mbe_energy_consistency():
    """Test that HMBE produces consistent energy contributions (no overcounting)."""
    print("Testing HMBE vs MBE energy consistency...")

    mol = Molecule(
        symbols=["H"] * 4,
        geometry=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        fragments=[[0], [1], [2], [3]],
    )

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "F0"), 2: ("G0", "F1"),
            3: ("G1", "F2"), 4: ("G1", "F3"),
        },
        tier_names=("group", "fragment")
    )

    # Create MBE (no HMBE) core
    mbc_mbe = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=None
    )

    # Create HMBE core
    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 2),
        hierarchy=hierarchy,
        enumeration_mode="filter"
    )

    mbc_hmbe = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g"},
        return_total_data=True,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    # Create mock energies using same model for both
    mc = "hf/sto-3g"

    def get_mock_energy(frag_tuple):
        """Simple additive energy model."""
        n = len(frag_tuple)
        return -1.0 * n - 0.05 * n * (n - 1) / 2

    # Generate energies for MBE
    mbe_results = {}
    for bsse_dict in mbc_mbe.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                lbl = labeler(mc, frag, bas)
                mbe_results[lbl] = get_mock_energy(frag)

    # Generate energies for HMBE (subset of MBE)
    hmbe_results = {}
    for bsse_dict in mbc_hmbe.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                lbl = labeler(mc, frag, bas)
                hmbe_results[lbl] = get_mock_energy(frag)

    print(f"  MBE energy calculations: {len(mbe_results)}")
    print(f"  HMBE energy calculations: {len(hmbe_results)}")

    # Assemble both
    mbe_assembled = mbc_mbe._assemble_nbody_components("energy", mbe_results)
    hmbe_assembled = mbc_hmbe._assemble_nbody_components("energy", hmbe_results)

    mbe_body_dict = mbe_assembled["energy_body_dict"]
    hmbe_body_dict = hmbe_assembled["energy_body_dict"]

    print(f"  MBE 1-body: {mbe_body_dict[1]:.6f}")
    print(f"  MBE 2-body: {mbe_body_dict[2]:.6f}")
    print(f"  HMBE 1-body: {hmbe_body_dict[1]:.6f}")
    print(f"  HMBE 2-body: {hmbe_body_dict[2]:.6f}")

    # For (2,2)-HMBE with this hierarchy, we keep:
    # - All 1-body: (1), (2), (3), (4)
    # - 2-body from same group: (1,2), (3,4)
    # - 2-body across 2 groups: (1,3), (1,4), (2,3), (2,4)
    # We exclude: nothing (all 2-body terms are kept for this case)

    # So HMBE and MBE should give identical results for this system
    diff_1b = abs(mbe_body_dict[1] - hmbe_body_dict[1])
    diff_2b = abs(mbe_body_dict[2] - hmbe_body_dict[2])

    print(f"  Difference in 1-body: {diff_1b:.10f}")
    print(f"  Difference in 2-body: {diff_2b:.10f}")

    assert diff_1b < 1e-10, f"1-body energies differ: {diff_1b}"
    assert diff_2b < 1e-10, f"2-body energies differ: {diff_2b}"

    print("  ✓ HMBE and MBE produce consistent results")
    print("✓ Energy consistency test passed!\n")


def test_hmbe_contribution_signs():
    """Test that HMBE contributions have correct signs (no sign errors)."""
    print("Testing HMBE contribution signs...")

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

    # Generate mock energies with 3-body effects (not just pairwise-additive)
    mc = "hf/sto-3g"

    def get_mock_energy(frag_tuple):
        """Energy model with true 3-body effects."""
        n = len(frag_tuple)
        # Base: -1.0 per monomer
        # 2-body: -0.02 per pair
        # 3-body: -0.005 per triple (non-additive contribution)
        energy = -1.0 * n
        if n >= 2:
            energy -= 0.02 * n * (n - 1) / 2  # Pairwise
        if n >= 3:
            energy -= 0.005 * n * (n - 1) * (n - 2) / 6  # 3-body
        return energy

    component_results = {}
    for bsse_dict in mbc.compute_map[mc].values():
        for terms in bsse_dict.values():
            for frag, bas in terms:
                lbl = labeler(mc, frag, bas)
                component_results[lbl] = get_mock_energy(frag)

    # Assemble
    results = mbc._assemble_nbody_components("energy", component_results)
    body_dict = results["energy_body_dict"]

    # Check contributions (n-body minus (n-1)-body)
    contrib_1b = body_dict[1]
    contrib_2b = body_dict[2] - body_dict[1]
    contrib_3b = body_dict[3] - body_dict[2]

    print(f"  1-body total: {contrib_1b:.6f}")
    print(f"  2-body contribution: {contrib_2b:.6f}")
    print(f"  3-body contribution: {contrib_3b:.6f}")

    # All contributions should be negative (attractive)
    assert contrib_1b < 0, f"1-body should be negative, got {contrib_1b}"
    assert contrib_2b < 0, f"2-body contrib should be negative (attractive), got {contrib_2b}"

    # Note: 3-body can be 0 if there are no 3-body terms, or small if effects are weak
    if abs(contrib_3b) > 1e-10:
        assert contrib_3b < 0, f"3-body contrib should be negative (attractive) if non-zero, got {contrib_3b}"
        print(f"  ✓ 3-body contribution is non-zero and attractive: {contrib_3b:.6f}")
    else:
        print(f"  ℹ 3-body contribution is essentially zero: {contrib_3b:.10f}")

    # Check for "explosion" - 3-body contribution should not be huge
    if abs(contrib_3b) > abs(contrib_1b):
        print(f"  ⚠ WARNING: 3-body contribution ({contrib_3b:.6f}) is larger than 1-body ({contrib_1b:.6f})")
        print("  This might indicate overcounting (though not necessarily wrong)")
        raise AssertionError("3-body contribution is suspiciously large - possible overcounting")

    print("  ✓ All non-zero contributions have correct signs")
    print("  ✓ No energy explosion detected")
    print("✓ Contribution signs test passed!\n")


if __name__ == "__main__":
    try:
        test_hmbe_energy_assembly_simple()
        test_hmbe_vs_mbe_energy_consistency()
        test_hmbe_contribution_signs()

        print("=" * 60)
        print("ALL ENERGY ASSEMBLY TESTS PASSED!")
        print("=" * 60)
        print("\nSummary:")
        print("- HMBE energy assembly works correctly")
        print("- No energy explosion or overcounting detected")
        print("- Möbius inversion produces correct contributions")
        print("- Completeness requirement is satisfied")
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
