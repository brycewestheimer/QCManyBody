"""
Basic HMBE Example: 2-Tier Hierarchical Many-Body Expansion

This example demonstrates how to use the Hierarchical Many-Body Expansion (HMBE)
to reduce computational cost for large systems by organizing fragments into
hierarchical tiers.

HMBE can reduce the number of calculations by 100-1000x compared to standard MBE
while maintaining good accuracy.
"""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import (
    FragmentHierarchy,
    HMBESpecification,
    SchengenSpecification,
)


def main():
    """Run a basic (2,3)-HMBE calculation on a 6-water system."""

    # ========================================================================
    # 1. Define the molecular system with fragments
    # ========================================================================
    # 6 water molecules arranged in a 3x2 hierarchy
    # 3 tier-1 groups (domains), each containing 2 water molecules

    mol = Molecule(
        symbols=["O", "H", "H"] * 6,  # 6 water molecules
        geometry=[
            # Water 1 (Group 0)
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.587],
            [0.0, -0.757, 0.587],
            # Water 2 (Group 0)
            [3.0, 0.0, 0.0],
            [3.0, 0.757, 0.587],
            [3.0, -0.757, 0.587],
            # Water 3 (Group 1)
            [0.0, 3.0, 0.0],
            [0.0, 3.757, 0.587],
            [0.0, 2.243, 0.587],
            # Water 4 (Group 1)
            [3.0, 3.0, 0.0],
            [3.0, 3.757, 0.587],
            [3.0, 2.243, 0.587],
            # Water 5 (Group 2)
            [0.0, 6.0, 0.0],
            [0.0, 6.757, 0.587],
            [0.0, 5.243, 0.587],
            # Water 6 (Group 2)
            [3.0, 6.0, 0.0],
            [3.0, 6.757, 0.587],
            [3.0, 5.243, 0.587],
        ],
        fragments=[
            [0, 1, 2],  # Fragment 1: Water 1
            [3, 4, 5],  # Fragment 2: Water 2
            [6, 7, 8],  # Fragment 3: Water 3
            [9, 10, 11],  # Fragment 4: Water 4
            [12, 13, 14],  # Fragment 5: Water 5
            [15, 16, 17],  # Fragment 6: Water 6
        ],
    )

    # ========================================================================
    # 2. Define the fragment hierarchy
    # ========================================================================
    # Organize 6 fragments into 3 tier-1 groups, each with 2 fragments
    # Group structure: G0: [1,2], G1: [3,4], G2: [5,6]

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers={
            1: ("G0", "G0_S0"),  # Group 0, Subgroup 0
            2: ("G0", "G0_S1"),  # Group 0, Subgroup 1
            3: ("G1", "G1_S0"),  # Group 1, Subgroup 0
            4: ("G1", "G1_S1"),  # Group 1, Subgroup 1
            5: ("G2", "G2_S0"),  # Group 2, Subgroup 0
            6: ("G2", "G2_S1"),  # Group 2, Subgroup 1
        },
        tier_names=("domain", "water"),
    )

    # ========================================================================
    # 3. Define HMBE specification
    # ========================================================================
    # (2,3)-HMBE: T_1 = 2 (max 2 domains), T_2 = 3 (max 3-body)
    # This filters out 3-body terms that span all 3 domains
    # With 6 waters: Standard MBE-3 = 20 3-body terms
    #                (2,3)-HMBE ≈ 12 3-body terms (40% reduction)

    hmbe_spec = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.1,  # Add back 10% of filtered terms (closest by distance)
            distance_metric="R2",  # Sum of squared distances
        ),
    )

    # ========================================================================
    # 4. Run the HMBE calculation
    # ========================================================================
    print("Setting up HMBE calculation...")
    print(f"  Hierarchy: {hierarchy.num_tiers} tiers")
    print(f"  Truncation orders: {hmbe_spec.truncation_orders}")
    print(f"  Schengen enabled: {hmbe_spec.schengen.enabled}")
    print()

    # Create ManyBodyCore with HMBE
    # Note: For actual calculations, use ManyBodyComputer.from_manybodyinput()
    # This example uses ManyBodyCore directly to demonstrate HMBE setup
    mbcore = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.cp],
        levels={
            1: "hf/sto-3g",  # 1-body at HF/STO-3G
            2: "mp2/sto-3g",  # 2-body at MP2/STO-3G
            3: "mp2/sto-3g",  # 3-body at MP2/STO-3G
        },
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec,
    )

    # ========================================================================
    # 5. Print HMBE statistics before running calculations
    # ========================================================================
    stats = mbcore.get_hmbe_statistics()

    print("HMBE Statistics:")
    print("-" * 60)
    for mc_level in stats["mbe_term_counts"]:
        mbe_count = stats["mbe_term_counts"][mc_level]
        hmbe_count = stats["hmbe_term_counts"][mc_level]
        reduction = stats["reduction_factors"][mc_level]

        print(f"Model chemistry: {mc_level}")
        print(f"  Standard MBE terms:  {mbe_count}")
        print(f"  HMBE terms:          {hmbe_count}")
        if reduction > 0:
            print(f"  Reduction factor:    {reduction:.2f}x")
            print(f"  Savings:             {100 * (1 - 1/reduction):.1f}%")
        else:
            print(f"  Reduction factor:    N/A (no terms at this level)")
        print()

    # ========================================================================
    # 6. (Optional) Run the actual calculations
    # ========================================================================
    # NOTE: This requires QCEngine and a QC program (Psi4, NWChem, etc.)
    # Uncomment to actually run the calculations:
    #
    # print("Running HMBE calculation...")
    # print("(This requires QCEngine and a quantum chemistry program)")
    # print()
    #
    # result = mbc.compute()
    #
    # print("Results:")
    # print("-" * 60)
    # print(f"CP-corrected interaction energy: {result.return_result} Hartree")
    # if hasattr(result.properties, "hmbe_metadata"):
    #     print(f"HMBE metadata: {result.properties.hmbe_metadata}")

    print("Example complete!")
    print()
    print("Key takeaways:")
    print("  - HMBE reduces the number of calculations while maintaining accuracy")
    print("  - Organize fragments into hierarchical tiers (domains → residues → atoms)")
    print("  - Use Schengen terms to add back important interface interactions")
    print("  - (2,3)-HMBE typically gives ~40-60% reduction for 2-tier systems")


if __name__ == "__main__":
    main()
