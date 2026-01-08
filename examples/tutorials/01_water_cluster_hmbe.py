"""
Tutorial 1: Water Cluster with (2,3)-HMBE

This tutorial demonstrates HMBE on a 12-water system organized into 3 spatial regions.
We use a 2-tier hierarchy (regions → waters) with (2,3)-HMBE truncation.

Key concepts:
- Creating a hierarchical molecule structure
- Defining fragment tiers
- Using HMBE to reduce computational cost
- Analyzing reduction factors

System: 12 water molecules
Hierarchy: 3 regions × 4 waters each (2-tier)
Truncation: (2,3)-HMBE (max 2 regions, max 3 fragments per term)
Expected reduction: ~1.6x vs standard MBE-3
"""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification


def create_water_cluster_molecule(num_regions=3, waters_per_region=4):
    """Create a simple water cluster arranged in spatial regions.

    Args:
        num_regions: Number of spatial regions (tier-1 groups)
        waters_per_region: Waters per region (fragments per tier-1 group)

    Returns:
        Molecule with waters organized into regions
    """
    total_waters = num_regions * waters_per_region

    # Create waters on a grid
    waters = []
    for region in range(num_regions):
        for local_idx in range(waters_per_region):
            x = region * 6.0  # Regions separated by 6 Angstrom
            y = local_idx * 3.0  # Waters within region separated by 3 Angstrom
            # Standard H2O geometry
            waters.append([
                ("O", [x, y, 0.0]),
                ("H", [x + 0.757, y, 0.587]),
                ("H", [x - 0.757, y, 0.587])
            ])

    # Flatten to QCElemental format
    symbols = []
    geometry = []
    fragments = []

    for i, water in enumerate(waters):
        start_idx = i * 3
        fragments.append([start_idx, start_idx + 1, start_idx + 2])
        for atom_symbol, coord in water:
            symbols.append(atom_symbol)
            geometry.append(coord)

    mol = Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments,
        molecular_charge=0.0,
        molecular_multiplicity=1,
    )

    print(f"Created molecule with {total_waters} water fragments")
    print(f"  Organized into {num_regions} regions")
    print(f"  {waters_per_region} waters per region")
    print()

    return mol


def create_water_cluster_hierarchy(num_regions=3, waters_per_region=4):
    """Create 2-tier hierarchy for water cluster.

    Tier 1: Spatial regions (coarse grouping)
    Tier 2: Individual water molecules (fine fragments)

    Args:
        num_regions: Number of tier-1 groups
        waters_per_region: Fragments per tier-1 group

    Returns:
        FragmentHierarchy for the water cluster
    """
    total_waters = num_regions * waters_per_region

    # Build fragment_tiers mapping
    fragment_tiers = {}
    for i in range(total_waters):
        frag_id = i + 1  # 1-indexed
        region_id = i // waters_per_region
        tier1_group = f"R{region_id}"  # R0, R1, R2
        tier2_name = f"W{i}"  # W0, W1, ..., W11

        fragment_tiers[frag_id] = (tier1_group, tier2_name)

    hierarchy = FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("region", "water")
    )

    print(f"Created 2-tier hierarchy:")
    print(f"  Tier 1 (coarse): {num_regions} regions")
    print(f"  Tier 2 (fine):   {total_waters} water molecules")
    print()

    return hierarchy


def run_hmbe_calculation(mol, hierarchy, truncation_orders=(2, 3)):
    """Run HMBE calculation and analyze results.

    Args:
        mol: Molecule with fragments
        hierarchy: Fragment hierarchy
        truncation_orders: HMBE truncation orders

    Returns:
        ManyBodyCore instance with results
    """
    print(f"Running {truncation_orders}-HMBE calculation...")
    print()

    # Create HMBE specification
    hmbe_spec = HMBESpecification(
        truncation_orders=truncation_orders,
        hierarchy=hierarchy
    )

    # Create ManyBodyCore (without running actual QC - just build compute map)
    mbcore = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.cp],
        levels={
            1: "hf/sto-3g",  # Use cheap method for demo
            2: "hf/sto-3g",
            3: "hf/sto-3g"
        },
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_spec
    )

    return mbcore


def analyze_results(mbcore):
    """Analyze and print HMBE statistics.

    Args:
        mbcore: ManyBodyCore instance
    """
    stats = mbcore.get_hmbe_statistics()

    print("="*70)
    print("HMBE STATISTICS")
    print("="*70)
    print()

    mc_level = list(stats['mbe_term_counts'].keys())[0]

    mbe_count = stats['mbe_term_counts'][mc_level]
    hmbe_count = stats['hmbe_term_counts'][mc_level]
    reduction = stats['reduction_factors'][mc_level]

    print(f"Model Chemistry: {mc_level}")
    print(f"Truncation:      {stats['truncation_orders']}")
    print(f"Tiers:           {stats['num_tiers']}")
    print()

    print("Term Counts:")
    print(f"  Standard MBE-{stats['truncation_orders'][-1]}:  {mbe_count:5d} calculations")
    print(f"  HMBE-{stats['truncation_orders']}:    {hmbe_count:5d} calculations")
    print(f"  Reduction factor:  {reduction:.2f}x")
    print()

    print("Computational Savings:")
    savings_pct = (1 - 1/reduction) * 100
    print(f"  You're running {savings_pct:.1f}% fewer calculations!")
    print()

    # Show breakdown by n-body order
    compute_map = mbcore.compute_map[mc_level]["cp"]

    print("Breakdown by n-body order:")
    for nbody in sorted(compute_map.keys()):
        count = len(compute_map[nbody])
        print(f"  {nbody}-body terms: {count:4d}")
    print()

    print("="*70)
    print()


def main():
    """Run complete tutorial workflow."""

    print()
    print("#" * 70)
    print("# Tutorial 1: Water Cluster with (2,3)-HMBE")
    print("#" * 70)
    print()

    # Step 1: Create molecule
    print("STEP 1: Creating molecule")
    print("-" * 70)
    mol = create_water_cluster_molecule(num_regions=3, waters_per_region=4)

    # Step 2: Create hierarchy
    print("STEP 2: Creating hierarchy")
    print("-" * 70)
    hierarchy = create_water_cluster_hierarchy(num_regions=3, waters_per_region=4)

    # Step 3: Run HMBE calculation
    print("STEP 3: Setting up HMBE calculation")
    print("-" * 70)
    mbcore = run_hmbe_calculation(mol, hierarchy, truncation_orders=(2, 3))

    # Step 4: Analyze results
    print("STEP 4: Analyzing results")
    print("-" * 70)
    analyze_results(mbcore)

    # Comparison: Try different truncation orders
    print("COMPARISON: Different truncation orders")
    print("="*70)
    print()

    for trunc in [(2, 3), (3, 3), (2, 4)]:
        if trunc[-1] > 3:
            # Need to update levels for 4-body
            continue

        print(f"Testing {trunc}-HMBE:")
        mbcore = run_hmbe_calculation(mol, hierarchy, truncation_orders=trunc)
        stats = mbcore.get_hmbe_statistics()
        mc_level = list(stats['mbe_term_counts'].keys())[0]
        reduction = stats['reduction_factors'][mc_level]
        hmbe_count = stats['hmbe_term_counts'][mc_level]
        mbe_count = stats['mbe_term_counts'][mc_level]
        print(f"  {hmbe_count} terms (vs {mbe_count} MBE) = {reduction:.2f}x reduction")
        print()

    print()
    print("Tutorial complete!")
    print()
    print("KEY TAKEAWAYS:")
    print("  1. HMBE organizes fragments into hierarchical tiers")
    print("  2. Truncation orders control filtering aggressiveness")
    print("  3. Always check reduction factor - should be >1.2x to be worthwhile")
    print("  4. For this 12-water system, (2,3)-HMBE gives ~1.6x speedup")
    print()
    print("NEXT STEPS:")
    print("  - Try with more regions (4-5) to see higher reduction")
    print("  - Experiment with Schengen terms (Tutorial 3)")
    print("  - Apply to your own molecular system")
    print()


if __name__ == "__main__":
    main()
